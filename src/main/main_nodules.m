clear; close all; clc;
clear all; clear functions;

%% CONFIGURAZIONE GENERALE

BASE_DIR  = fullfile(pwd, 'data', 'LIDC-IDRI');
DEBUG_DIR = 'debug_images';

MAX_SLICES_PER_PATIENT = 1;   % analizza al massimo 1 slice per paziente (la slice equatoriale)
MIN_EDGEMAP_PTS        = 3;   % numero minimo di punti del contorno per accettare un'annotazione XML

if ~exist(DEBUG_DIR, 'dir')
    mkdir(DEBUG_DIR);
end

patient_entries = dir(fullfile(BASE_DIR, 'LIDC-IDRI-*'));
patient_entries = patient_entries([patient_entries.isdir]);

%% INIZIALIZZAZIONE DELLE VARIABILI DI RISULTATO

all_dice   = [];    % vettore dei Dice score per tutti i casi
all_prec   = [];    % vettore della Precision per tutti i casi
all_sens   = [];    % vettore della Sensitivity per tutti i casi
all_labels = {};    % etichette testuali (paziente + slice) per ogni caso

% ── Dice stratificato per posizione topografica ──
dice_central  = [];   % noduli centrali (all'interno del parenchima)
dice_juxta    = [];   % noduli juxta-pleurali (adesi alla pleura)
dice_isolated = [];   % noduli isolati (circondati da aria)

% ── Dice stratificato per dimensione ──
dice_small    = [];   % noduli piccoli (<10 mm)
dice_medium   = [];   % noduli medi (10–20 mm)
dice_large    = [];   % noduli grandi (>20 mm)

% ── Dice stratificato per carattere densitometrico ──
dice_solid    = [];   % noduli solidi / part-solid
dice_ggo      = [];   % noduli ground-glass

% ── Precision stratificata (stesse categorie) ──
prec_central  = [];
prec_juxta    = [];
prec_isolated = [];
prec_small    = [];
prec_medium   = [];
prec_large    = [];
prec_solid    = [];
prec_ggo      = [];

% ── Sensitivity stratificata (stesse categorie) ──
sens_central  = [];
sens_juxta    = [];
sens_isolated = [];
sens_small    = [];
sens_medium   = [];
sens_large    = [];
sens_solid    = [];
sens_ggo      = [];

fprintf('Found %d patients.\n\n', length(patient_entries));

%% CICLO PRINCIPALE DI VALUTAZIONE

for p = 1:length(patient_entries)   % itera su ogni paziente nel dataset

    %% IDENTIFICAZIONE DEL PAZIENTE
    patient_id   = patient_entries(p).name;               % nome della cartella del paziente
    patient_path = fullfile(BASE_DIR, patient_id);        % percorso completo della cartella paziente

    % Cerca il file XML con le annotazioni e la cartella DICOM nella struttura paziente
    [xml_file, dcm_dir] = find_series(patient_path);
    if isempty(xml_file)   % se non c'è un XML, il paziente non ha annotazioni
        fprintf('[%s] No XML found — skipping.\n', patient_id);
        continue;
    end

    %% PARSING DELLE ANNOTAZIONI XML
    fprintf('[%s] Parsing XML...\n', patient_id);
    % Legge il file XML e restituisce la struttura gt_db con le annotazioni dei radiologi
    gt_db = parse_lidc_xml(xml_file, MIN_EDGEMAP_PTS);

    if isempty(gt_db)   % nessuna ROI valida trovata nell'XML
        fprintf('[%s] No valid nodule ROIs found — skipping.\n', patient_id);
        continue;
    end

    %% FILTRO DI CONSENSO STRETTO 4/4
    % Mantiene solo le slice annotate da tutti e 4 i radiologi
    mask_4 = arrayfun(@(x) length(x.reader_masks) == 4, gt_db);  % è un vettore logico: true per le slice con 4 annotazioni
   
    if any(mask_4)
        gt_db = gt_db(mask_4);   % filtra: mantieni solo le slice con consenso totale
    else
        fprintf('[%s] No 4/4 consensus slices — skipping.\n', patient_id);
        continue;
    end

    fprintf('[%s] %d valid nodule slices in XML.\n', patient_id, length(gt_db));

    %% CARICAMENTO DEI FILE DICOM
    dcm_files = dir(fullfile(dcm_dir, '*.dcm'));   % lista di tutti i file .dcm nella cartella serie
    if isempty(dcm_files)
        fprintf('[%s] No DICOM files found — skipping.\n', patient_id);
        continue;
    end

    fprintf('[%s] Loading volume (%d DICOMs)...\n', patient_id, length(dcm_files));
    % Carica il volume 3D in HU, le posizioni z delle slice, pixel spacing e spessore di slice
    [vol_hu, z_positions, pixel_spacing, slice_thickness] = ...
        load_patient_volume(dcm_dir, dcm_files);

    if isempty(vol_hu)   % se il caricamento è fallito
        fprintf('[%s] Failed to load volume — skipping.\n', patient_id);
        continue;
    end

    [rows, cols, ~] = size(vol_hu);   % estrae dimensioni del volume (righe, colonne, n_slice)
    fprintf('[%s] Volume: %dx%dx%d  px=%.3fmm\n', ...
        patient_id, rows, cols, size(vol_hu,3), pixel_spacing);

    %% SELEZIONE DELLA SLICE EQUATORIALE
    % Sceglie al massimo MAX_SLICES_PER_PATIENT slice rappresentativa per il nodulo
    selected = select_central_slices(gt_db, MAX_SLICES_PER_PATIENT);

    %% CICLO SULLE SLICE SELEZIONATE
    for s = 1:length(selected)

        %% ESTRAZIONE DELLA GROUND TRUTH CORRENTE
        gt_z_pos     = selected(s).z_pos;          % posizione z della slice annotata
        reader_masks = selected(s).reader_masks;   % cell array con le 4 maschere dei radiologi

        % Trova l'indice della slice DICOM più vicina alla z annotata
        [~, z_idx] = min(abs(z_positions - gt_z_pos));
        hu_slice   = vol_hu(:,:,z_idx);   % estrae la slice 2D in HU

        %% CONTESTO 2.5D: SLICE ADIACENTI A ≈5 mm
        n_vol           = size(vol_hu, 3);                          % numero totale di slice nel volume
        z_offset_slices = max(1, round(5.0 / slice_thickness));     % quante slice corrispondono a 5 mm
        z_above         = z_idx - z_offset_slices;                  % indice della slice superiore
        z_below         = z_idx + z_offset_slices;                  % indice della slice inferiore

        if z_above >= 1 && z_below <= n_vol   % controlla che le slice adiacenti siano nel volume
            hu_above = vol_hu(:,:, z_above);   % slice superiore (usata per filtro anti-vasi)
            hu_below = vol_hu(:,:, z_below);   % slice inferiore (usata per filtro anti-vasi)
        else
            hu_above = [];   % slice fuori volume: si disabilita il filtro 2.5D
            hu_below = [];
        end

        %% NORMALIZZAZIONE DIMENSIONALE DELLE MASCHERE GT
        % Le annotazioni XML usano coordinate 512×512; la slice può avere dimensioni diverse
        for r = 1:length(reader_masks)
            rm = reader_masks{r};
            if size(rm,1) ~= rows || size(rm,2) ~= cols
                % Ridimensiona la maschera al formato della slice DICOM con interpolazione nearest
                reader_masks{r} = imresize(rm, [rows cols], 'nearest') > 0;
            end
        end

        %% FILTRO DI FLEISCHNER (ESCLUDE NODULI < 6 mm)
        % Il filtro è applicato prima sull'unione delle maschere GT per ottenere il diametro

        gt_union_f = false(rows, cols);
        for r_f = 1:length(reader_masks)
            gt_union_f = gt_union_f | reader_masks{r_f};   % unione delle annotazioni dei 4 radiologi
        end

        cc_gtf = bwconncomp(gt_union_f, 8);   % componenti connesse dell'unione
        if cc_gtf.NumObjects == 0
            continue;   % nessuna regione annotata: salta
        end

        % Calcola l'area del componente più grande e stima il diametro equivalente (cerchio)
        comp_areas_f = cellfun(@numel, cc_gtf.PixelIdxList);
        gt_diam_f = 2 * sqrt(double(max(comp_areas_f)) / pi) * pixel_spacing;

        if gt_diam_f < 6.0   % nodulo sotto la soglia di Fleischner (< 6 mm)
            fprintf('  [FLEISCHNER SKIP] GT_diam=%.1fmm < 6mm — nodulo sub-soglia escluso.\n', gt_diam_f);
            continue;   % salta questo caso: non clinicamente rilevante
        end

        %% RIMOZIONE DEI COMPONENTI SOTTO-SOGLIA DAI SINGOLI LETTORI
        % Applica lo stesso filtro di Fleischner a ogni singola maschera di radiology reader

        for r_f2 = 1:length(reader_masks)
            rm = reader_masks{r_f2};
            cc_rm = bwconncomp(rm, 8);       % componenti connesse nella maschera del singolo reader
            rm_filt = false(rows, cols);      % nuova maschera filtrata

            for kk = 1:cc_rm.NumObjects
                % Diametro equivalente del k-esimo componente del reader
                d_kk = 2 * sqrt(double(numel(cc_rm.PixelIdxList{kk})) / pi) * pixel_spacing;
                if d_kk >= 6.0
                    rm_filt(cc_rm.PixelIdxList{kk}) = true;   % mantieni solo componenti ≥ 6 mm
                end
            end

            reader_masks{r_f2} = rm_filt;   % sostituisce la maschera originale con la versione filtrata
        end

        %% ELIMINAZIONE DELLE MASCHERE VUOTE
        % Dopo Fleischner, alcuni reader potrebbero avere maschere completamente vuote
        non_empty = cellfun(@(m) any(m(:)), reader_masks);   % vettore logico: true se la maschera non è vuota
        reader_masks = reader_masks(non_empty);               % mantiene solo le maschere non vuote

        if isempty(reader_masks)
            continue;   % se non rimangono maschere: salta
        end

        %% SECONDO CONTROLLO DI CONSENSO 4/4
        if length(reader_masks) < 4
            fprintf('  [SKIP 4/4] Only %d readers remain after Fleischner filter.\n', length(reader_masks));
            continue;   % non si ha più il consenso totale: salta
        end

        %% COSTRUZIONE DELLE MASCHERE COMBINATE PER IL BEST-MATCH
        % Per gestire la variabilità inter-osservatore, si aggiungono le
        % varie combinazioni

        original_masks = reader_masks;         % salva le 4 maschere originali
        n_readers = length(original_masks);    

        if n_readers > 1
            % Majority vote: somma le maschere e accetta i pixel marcati da ≥ n/2 lettori
            sum_mask = zeros(rows, cols);
            for r = 1:n_readers
                sum_mask = sum_mask + double(original_masks{r});   % accumula il voto di ogni reader
            end
            consensus_mask = sum_mask >= (n_readers / 2);    % soglia a metà dei lettori (≥2/4)
            reader_masks{end+1} = consensus_mask;             % aggiunge la consensus mask alla pool

            % Unioni pairwise: ogni coppia di lettori
            for r1 = 1:n_readers-1
                for r2 = r1+1:n_readers
                    pair_mask = original_masks{r1} | original_masks{r2};   % unione dei due lettori
                    reader_masks{end+1} = pair_mask;                        % aggiunge alla pool
                end
            end
        end

        %% SEGMENTAZIONE AUTOMATICA
        % Chiama l'algoritmo principale con la slice corrente e il contesto 2.5D
        pred_mask = find_nodules_advanced(hu_slice, hu_above, hu_below, pixel_spacing);

        %% VALUTAZIONE BEST-MATCH
        % Il Dice finale è il migliore tra la predizione e tutte le GT candidate nella pool
        best_dice   = 0;
        best_prec   = 0;
        best_sens   = 0;
        best_reader = 1;   % indice della GT che ha realizzato il best match

        for r = 1:length(reader_masks)
            [d, p_, s_] = compute_metrics(reader_masks{r}, pred_mask);   % calcola Dice, Prec, Sens
            if d > best_dice       % se questo reader dà un Dice migliore
                best_dice   = d;
                best_prec   = p_;
                best_sens   = s_;
                best_reader = r;   % aggiorna il best match
            end
        end

        % Aggiunge le metriche ai vettori globali
        all_dice(end+1)   = best_dice;
        all_prec(end+1)   = best_prec;
        all_sens(end+1)   = best_sens;
        all_labels{end+1} = sprintf('%s_s%d', patient_id, s);   % etichetta testuale

        %% COSTRUZIONE DELLA GT DI VISUALIZZAZIONE
        gt_union_disp = false(rows, cols);
        for r_d = 1:length(original_masks)
            gt_union_disp = gt_union_disp | original_masks{r_d};   % unione delle 4 maschere originali
        end

        cc_disp       = bwconncomp(gt_union_disp, 8);
        areas_disp    = sort(cellfun(@numel, cc_disp.PixelIdxList), 'descend');   % aree ordinate
        gt_union_area = sum(gt_union_disp(:));   % area totale (pixel) della GT unita

        % Costruisce una stringa con i diametri di ogni componente GT
        diam_parts = {};
        for nd = 1:length(areas_disp)
            d_nd = 2 * sqrt(double(areas_disp(nd)) / pi) * pixel_spacing;
            if d_nd < 6.0
                diam_parts{end+1} = sprintf('[SKIP:%.1fmm]', d_nd);   % componente sotto Fleischner
            else
                diam_parts{end+1} = sprintf('%.1fmm', d_nd);          % componente valido
            end
        end
        gt_nodules_str = strjoin(diam_parts, ' + ');   % unisce i diametri in una stringa

        % Stampa il riepilogo della slice corrente
        fprintf(['  slice %d/%d  z=%.1f (vol_idx=%d)  Best-Dice = %.4f  ' ...
                 '(reader %d/%d)  GT_area=%dpx  GT=%s\n'], ...
            s, length(selected), gt_z_pos, z_idx, best_dice, ...
            best_reader, length(reader_masks), gt_union_area, gt_nodules_str);

        gt_best = reader_masks{best_reader};   % maschera GT che ha realizzato il best Dice

        %% ANALISI STRATIFICATA
        % Classifica il caso nelle categorie di analisi basandosi sulla GT del best match

        cc_gt = bwconncomp(gt_best, 8);
        if cc_gt.NumObjects > 0

            % Proprietà del componente GT più grande (in area)
            stats_gt = regionprops(cc_gt, hu_slice, 'Area', 'Centroid', 'MeanIntensity');
            [~, mx] = max([stats_gt.Area]);   % indice del componente più grande
            gt_s = stats_gt(mx);

            gt_area_px = gt_s.Area;
            gt_diam_mm = 2 * sqrt(double(gt_area_px) / pi) * pixel_spacing;   % diametro equivalente in mm

            % Estrae la maschera polmonare per calcolare la distanza dal bordo
            [lung_mask_eval, ~] = extract_lung_mask_2d(hu_slice, rows, cols);
            lung_dist_eval = bwdist(~lung_mask_eval);   % mappa distanze dal bordo polmonare

            % Coordinate del centroide della GT (clampate ai bordi dell'immagine)
            cy_gt = max(1, min(rows, round(gt_s.Centroid(2))));
            cx_gt = max(1, min(cols, round(gt_s.Centroid(1))));
            dist_gt = lung_dist_eval(cy_gt, cx_gt);   % distanza del centroide dalla pleura (pixel)

            % Costruisce la maschera binaria del componente GT principale
            gt_blob = false(rows, cols);
            gt_blob(cc_gt.PixelIdxList{mx}) = true;

            % Anello contestuale attorno alla GT (stesso principio dell'halo in find_nodules_advanced)
            gt_ring = imdilate(gt_blob, strel('disk', max(2, round(5.0/pixel_spacing)))) & ...
                      ~imdilate(gt_blob, strel('disk', max(1, round(3.0/pixel_spacing))));
            gt_ring = gt_ring & ~gt_blob;   % esclude il blob stesso dall'anello

            if any(gt_ring(:))
                gt_ring_HU = mean(hu_slice(gt_ring));   % HU medio dell'anello attorno alla GT
            else
                gt_ring_HU = -550;
            end

            % Halo score della GT (stesso calcolo di find_nodules_advanced)
            gt_halo = max(0.01, min(1.0, (-gt_ring_HU - 100.0)/600.0));

            %% STRATIFICAZIONE TOPOGRAFICA
            if dist_gt < 15              % centroide vicino al bordo polmonare → juxta
                dice_juxta(end+1) = best_dice;
                prec_juxta(end+1) = best_prec;
                sens_juxta(end+1) = best_sens;
            elseif gt_halo > 0.85        % molto circondato da aria → isolato
                dice_isolated(end+1) = best_dice;
                prec_isolated(end+1) = best_prec;
                sens_isolated(end+1) = best_sens;
            else                         % né juxta né isolato → centrale
                dice_central(end+1) = best_dice;
                prec_central(end+1) = best_prec;
                sens_central(end+1) = best_sens;
            end

            %% STRATIFICAZIONE DIMENSIONALE
            if gt_diam_mm < 10           % piccolo (<10 mm)
                dice_small(end+1) = best_dice;
                prec_small(end+1) = best_prec;
                sens_small(end+1) = best_sens;
            elseif gt_diam_mm < 20       % medio (10–20 mm)
                dice_medium(end+1) = best_dice;
                prec_medium(end+1) = best_prec;
                sens_medium(end+1) = best_sens;
            else                         % grande (>20 mm)
                dice_large(end+1) = best_dice;
                prec_large(end+1) = best_prec;
                sens_large(end+1) = best_sens;
            end

            %% STRATIFICAZIONE DENSITOMETRICA
            if gt_s.MeanIntensity < -100   % HU basso → ground-glass
                dice_ggo(end+1) = best_dice;
                prec_ggo(end+1) = best_prec;
                sens_ggo(end+1) = best_sens;
            else                           % HU più alto → solido / part-solid
                dice_solid(end+1) = best_dice;
                prec_solid(end+1) = best_prec;
                sens_solid(end+1) = best_sens;
            end
        end

        %% SALVATAGGIO DELLE VISUALIZZAZIONI DI DEBUG
        % Genera e salva un'immagine con 3 pannelli: slice originale, GT, predizione
        save_debug_image(DEBUG_DIR, patient_id, s, hu_slice, gt_best, pred_mask, best_dice);
    end
end

%% RIEPILOGO GLOBALE

fprintf('\n');
fprintf('══════════════════════════════════════\n');
fprintf('  RESULTS\n');
fprintf('══════════════════════════════════════\n');

if isempty(all_dice)
    fprintf('  No slices evaluated.\n');   % nessun caso analizzato
else
    fprintf('  Slices evaluated : %d\n',  length(all_dice));
    fprintf('  Mean  Dice       : %.4f\n', mean(all_dice));
    fprintf('  Median Dice      : %.4f\n', median(all_dice));
    fprintf('  Mean  Precision  : %.4f\n', mean(all_prec));
    fprintf('  Mean  Sensitivity: %.4f\n', mean(all_sens));
    fprintf('  Min   Dice       : %.4f\n', min(all_dice));
    fprintf('  Max   Dice       : %.4f\n', max(all_dice));
    fprintf('  Dice > 0.60      : %d/%d slices\n', sum(all_dice > 0.60), length(all_dice));
    fprintf('  Zeros            : %d/%d slices\n', sum(all_dice == 0),   length(all_dice));

    valid_idx = (all_dice > 0);   % indici dei casi dove l'algoritmo ha segmentato qualcosa
    if any(valid_idx)
        fprintf('\n  Mean metrics on detected cases only (Dice > 0):\n');
        fprintf('  Mean  Dice       : %.4f\n', mean(all_dice(valid_idx)));
        fprintf('  Mean  Precision  : %.4f\n', mean(all_prec(valid_idx)));
        fprintf('  Mean  Sensitivity: %.4f\n', mean(all_sens(valid_idx)));
    end

    fprintf('\n  Per-slice results:\n');
    for i = 1:length(all_dice)
        flag = '';
        if all_dice(i) < 0.30
            flag = '  <- LOW';   % segnala i casi con performance molto bassa
        end
        fprintf('    %-30s  Dice=%.4f  Prec=%.4f  Sens=%.4f%s\n', ...
            all_labels{i}, all_dice(i), all_prec(i), all_sens(i), flag);
    end
end
fprintf('══════════════════════════════════════\n');

%% RIEPILOGO STRATIFICATO

fprintf('\n');
fprintf('══════════════════════════════════════\n');
fprintf('  STRATIFIED RESULTS\n');
fprintf('══════════════════════════════════════\n');

% Stampa le statistiche per ogni sottogruppo
print_group_stats('CENTRAL',            dice_central,  prec_central,  sens_central);
print_group_stats('JUXTA',              dice_juxta,    prec_juxta,    sens_juxta);
print_group_stats('ISOLATED',           dice_isolated, prec_isolated, sens_isolated);

print_group_stats('SMALL <10mm',        dice_small,    prec_small,    sens_small);
print_group_stats('MEDIUM 10-20mm',     dice_medium,   prec_medium,   sens_medium);
print_group_stats('LARGE >20mm',        dice_large,    prec_large,    sens_large);

print_group_stats('SOLID / PART-SOLID', dice_solid,    prec_solid,    sens_solid);
print_group_stats('GGO-LIKE',           dice_ggo,      prec_ggo,      sens_ggo);

fprintf('══════════════════════════════════════\n');


%% ── FUNZIONI LOCALI ──────────────────────────────────────────────────────────

function [vol_hu, z_positions, pixel_spacing, slice_thickness] = ...
    load_patient_volume(dcm_dir, dcm_files)
% Carica la serie DICOM, ordina le slice per posizione z e converte in HU.

vol_hu = [];              % volume 3D output (vuoto in caso di errore)
z_positions = [];         % posizioni z (mm) di ogni slice
pixel_spacing = 0.70;     % default: 0.70 mm/pixel
slice_thickness = 2.0;    % default: 2.0 mm di spessore

n = length(dcm_files);    % numero di file DICOM trovati
if n == 0
    return;   % nessun file: esce
end

% Struttura temporanea per memorizzare le slice prima dell'ordinamento
slices = struct('z',{},'img',{},'slope',{},'intercept',{},'px_sp',{});

for i = 1:n
    fpath = fullfile(dcm_dir, dcm_files(i).name);   % percorso del file DICOM i-esimo
    try
        info    = dicominfo(fpath, 'UseDictionaryVR', true);   % legge l'header DICOM
        img_raw = double(dicomread(info));                      % legge i pixel grezzi
    catch
        continue;   % se il file è corrotto o non leggibile: salta
    end

    z = NaN;   % posizione z inizialmente non definita
    if isfield(info,'ImagePositionPatient') && numel(info.ImagePositionPatient)>=3
        z = double(info.ImagePositionPatient(3));   % coordinata z dall'header DICOM (campo standard)
    elseif isfield(info,'SliceLocation')
        z = double(info.SliceLocation);             % fallback: usa SliceLocation se presente
    end
    if isnan(z)
        continue;   % posizione z non recuperabile: salta questa slice
    end

    slope = 1;       % pendenza per la conversione a HU (default = 1)
    intercept = 0;   % intercetta per la conversione a HU (default = 0)
    if isfield(info,'RescaleSlope')
        slope = double(info.RescaleSlope);
    end
    if isfield(info,'RescaleIntercept')
        intercept = double(info.RescaleIntercept);
    end

    px_sp = 0.70;   % spaziatura pixel default
    if isfield(info,'PixelSpacing') && ~isempty(info.PixelSpacing)
        px_sp = double(info.PixelSpacing(1));   % legge il pixel spacing reale dall'header
    end

    k = length(slices) + 1;     % indice della prossima slice da inserire
    slices(k).z         = z;
    slices(k).img       = img_raw;
    slices(k).slope     = slope;
    slices(k).intercept = intercept;
    slices(k).px_sp     = px_sp;
end

if isempty(slices)
    return;   % nessuna slice valida caricata
end

% Ordina le slice in base alla posizione z (asse cranio-caudale)
[z_sorted, sort_idx] = sort([slices.z]);
slices = slices(sort_idx);

pixel_spacing = slices(1).px_sp;   % usa il pixel spacing della prima slice ordinata

if length(z_sorted) > 1
    dz = abs(diff(z_sorted));                  % differenze consecutive tra posizioni z
    slice_thickness = max(0.5, median(dz));    % spessore mediano (robusto agli outlier)
else
    slice_thickness = 2.0;   % se c'è una sola slice: usa il default
end

[rows, cols] = size(slices(1).img);   % dimensioni spaziali (tutte le slice devono essere uguali)
n_slices = length(slices);

vol_hu      = zeros(rows, cols, n_slices, 'double');   % alloca il volume 3D
z_positions = z_sorted;                                 % salva le posizioni z ordinate

for i = 1:n_slices
    sl    = slices(i);
    img_i = sl.img;
    if size(img_i,1) ~= rows || size(img_i,2) ~= cols
        img_i = imresize(img_i, [rows cols], 'nearest');   % riallinea se le dimensioni differiscono
    end
    vol_hu(:,:,i) = double(img_i) .* sl.slope + sl.intercept;   % conversione a HU: HU = pixel*slope + intercept
end
end


function [xml_file, dcm_dir] = find_series(patient_path)
% Cerca nella struttura ad albero del paziente il file XML e la cartella DICOM.

xml_file = '';   % inizializza come vuoto
dcm_dir  = '';

% Lista le sottocartelle di primo livello (date di studio)
date_entries = dir(patient_path);
date_entries = date_entries([date_entries.isdir] & ...
    ~strcmp({date_entries.name},'.') & ~strcmp({date_entries.name},'..'));
% Filtra: solo cartelle, esclude '.' e '..'

for d = 1:length(date_entries)
    date_path = fullfile(patient_path, date_entries(d).name);   % percorso della data di studio

    % Lista le sottocartelle di secondo livello (serie DICOM)
    series_entries = dir(date_path);
    series_entries = series_entries([series_entries.isdir] & ...
        ~strcmp({series_entries.name},'.') & ~strcmp({series_entries.name},'..'));

    for s = 1:length(series_entries)
        series_path = fullfile(date_path, series_entries(s).name);   % percorso di una serie
        xml_list = dir(fullfile(series_path,'*.xml'));                % cerca file .xml in questa serie
        if ~isempty(xml_list)
            xml_file = fullfile(series_path, xml_list(1).name);   % prende il primo XML trovato
            dcm_dir  = series_path;                               % la cartella della serie è anche quella dei DICOM
            return;   % trovato: esce subito dalla funzione
        end
    end
end
end


function gt_db = parse_lidc_xml(xml_file, min_edge_pts)
% Legge il file XML LIDC-IDRI e costruisce la struttura gt_db con le annotazioni.

% Struttura vuota con i campi attesi
gt_db = struct('sop_uid',{},'reader_masks',{},'gt_union',{},'z_pos',{},'is_central',{});

try
    doc = xmlread(xml_file);   % apre e parsa il file XML con il parser Java integrato di MATLAB
catch
    warning('parse_lidc_xml: cannot read %s', xml_file);
    return;   % se il file XML è corrotto: restituisce struttura vuota
end

IMG_SIZE = 512;   % dimensione standard delle immagini LIDC-IDRI

sessions   = doc.getElementsByTagName('readingSession');   % lista delle sessioni di lettura (una per reader)
n_sessions = sessions.getLength();                         % numero di sessioni (di solito 4)

% Una pool per ogni sessione: mappa da SOP_UID → lista di ROI di quel reader
sess_pools = cell(1, n_sessions);
for si = 0:n_sessions-1
    sess_pools{si+1} = containers.Map('KeyType','char','ValueType','any');
end

z_pool       = containers.Map('KeyType','char','ValueType','double');   % SOP_UID → posizione z
central_pool = containers.Map('KeyType','char','ValueType','double');   % SOP_UID → n_edge del nodulo centrale

for si = 0:n_sessions-1   % ciclo sulle sessioni (0-indexed per la Java API)
    session = sessions.item(si);                                    % sessione di lettura i-esima
    nodules = session.getElementsByTagName('unblindedReadNodule'); % lista dei noduli annotati
    pool    = sess_pools{si+1};

    for ni = 0:nodules.getLength()-1   % ciclo su ogni nodulo annotato nella sessione
        nodule = nodules.item(ni);
        rois   = nodule.getElementsByTagName('roi');   % ROI del nodulo (una per slice)
        nodule_items = {};

        for ri = 0:rois.getLength()-1   % ciclo su ogni ROI (slice) del nodulo
            roi = rois.item(ri);

            sop_uid = xml_text(roi, 'imageSOP_UID');   % identificatore univoco della slice DICOM
            if isempty(sop_uid)
                continue;
            end

            z_str = xml_text(roi, 'imageZposition');   % posizione z in mm
            if isempty(z_str)
                continue;
            end
            z_pos = str2double(z_str);   % converte da stringa a numero

            inc_str   = xml_text(roi, 'inclusion');           % 'TRUE' se la ROI è da includere, 'FALSE' da escludere
            inclusion = strcmpi(strtrim(inc_str), 'TRUE');    % converte in booleano

            edges  = roi.getElementsByTagName('edgeMap');   % punti del contorno della ROI
            n_edge = edges.getLength();
            if n_edge < min_edge_pts   % contorno con troppo pochi punti: non affidabile
                continue;
            end

            xv = zeros(1,n_edge);   % vettore delle coordinate x dei punti del contorno
            yv = zeros(1,n_edge);   % vettore delle coordinate y
            ok = true;              % flag di validità

            for ei = 0:n_edge-1    % ciclo su ogni punto del contorno
                em = edges.item(ei);
                xn = em.getElementsByTagName('xCoord');
                yn = em.getElementsByTagName('yCoord');
                if xn.getLength()==0 || yn.getLength()==0
                    ok = false;   % coordinata mancante: contorno invalido
                    break;
                end
                xv(ei+1) = str2double(char(xn.item(0).getTextContent()));   % legge x
                yv(ei+1) = str2double(char(yn.item(0).getTextContent()));   % legge y
            end
            if ~ok
                continue;   % contorno non valido: salta questa ROI
            end

            % Memorizza le coordinate del poligono e il flag di inclusione
            entry.x = xv;
            entry.y = yv;
            entry.inclusion = inclusion;

            % Inizializza la pool per questa slice se non esiste ancora
            if ~isKey(pool, sop_uid)
                pool(sop_uid) = {};
                if ~isKey(z_pool, sop_uid)
                    z_pool(sop_uid) = z_pos;   % registra la posizione z della slice
                end
            end

            % Aggiunge la ROI alla pool della slice
            tmp = pool(sop_uid);
            tmp{end+1} = entry;
            pool(sop_uid) = tmp;

            % Aggiorna anche il dizionario nodule_items per trovare il nodulo "centrale"
            item.sop_uid = sop_uid;
            item.n_edge  = n_edge;
            nodule_items{end+1} = item;
        end

        if isempty(nodule_items)
            continue;   % nessuna ROI valida per questo nodulo
        end

        % Trova la slice del nodulo con più punti di contorno (≈ slice equatoriale)
        best_n = 0;
        best_sop = '';
        for k = 1:length(nodule_items)
            if nodule_items{k}.n_edge > best_n
                best_n = nodule_items{k}.n_edge;
                best_sop = nodule_items{k}.sop_uid;   % aggiorna la slice con il contorno più ricco
            end
        end

        if ~isempty(best_sop)
            % Registra la slice con più punti come "centrale" (slice equatoriale)
            if ~isKey(central_pool,best_sop) || best_n > central_pool(best_sop)
                central_pool(best_sop) = best_n;
            end
        end

        sess_pools{si+1} = pool;   % salva la pool aggiornata
    end
end

% Raccoglie tutti gli SOP_UID unici (unione tra tutti i reader)
all_uids = {};
for si = 1:n_sessions
    ks = keys(sess_pools{si});
    all_uids = union(all_uids, ks);   % aggiunge le chiavi della sessione i-esima senza duplicati
end

% Per ogni slice annotata, costruisce la struttura gt_db
for i = 1:length(all_uids)
    uid = all_uids{i};

    reader_masks = {};   % lista delle maschere binarie (una per reader che ha annotato questa slice)

    for si = 1:n_sessions
        pool = sess_pools{si};
        if ~isKey(pool, uid)
            continue;   % questo reader non ha annotato questa slice
        end

        rois = pool(uid);                        % lista delle ROI di questo reader per questa slice
        mask = false(IMG_SIZE, IMG_SIZE);        % maschera binaria vuota 512×512

        for j = 1:length(rois)
            r  = rois{j};
            xc = min(max(r.x,1),IMG_SIZE);       % clamp delle coordinate x nel range [1, 512]
            yc = min(max(r.y,1),IMG_SIZE);       % clamp delle coordinate y nel range [1, 512]
            try
                filled = poly2mask(xc, yc, IMG_SIZE, IMG_SIZE);   % rasterizza il poligono in maschera binaria
            catch
                continue;   % poly2mask può fallire con poligoni degeneri
            end

            if r.inclusion
                mask = mask | filled;     % ROI da includere: aggiunge i pixel
            else
                mask = mask & ~filled;    % ROI da escludere: rimuove i pixel (annotazione negativa)
            end
        end

        if any(mask(:))
            reader_masks{end+1} = mask;   % aggiunge la maschera di questo reader
        end
    end

    if isempty(reader_masks)
        continue;   % nessun reader ha prodotto una maschera non vuota: salta
    end

    % Calcola l'unione delle maschere di tutti i reader per questa slice
    gt_union = false(IMG_SIZE, IMG_SIZE);
    for r = 1:length(reader_masks)
        gt_union = gt_union | reader_masks{r};
    end
    if ~any(gt_union(:))
        continue;   % unione vuota: salta
    end

    if ~isKey(z_pool, uid)
        continue;   % posizione z non disponibile per questa slice: salta
    end

    % Aggiunge una nuova voce alla struttura gt_db
    k = length(gt_db) + 1;
    gt_db(k).sop_uid      = uid;                        % identificatore DICOM della slice
    gt_db(k).reader_masks = reader_masks;               % maschere dei singoli reader
    gt_db(k).gt_union     = gt_union;                   % unione di tutte le maschere
    gt_db(k).z_pos        = z_pool(uid);                % posizione z in mm
    gt_db(k).is_central   = isKey(central_pool, uid);  % true se è la slice equatoriale del nodulo
end
end


function txt = xml_text(node, tag)
% Restituisce il contenuto testuale del primo elemento XML con il tag specificato.

txt = '';
list = node.getElementsByTagName(tag);   % cerca tutti gli elementi con quel tag nel nodo
if list.getLength()==0
    return;   % tag non trovato: restituisce stringa vuota
end
txt = strtrim(char(list.item(0).getTextContent()));   % prende il testo del primo risultato e toglie spazi
end


function selected = select_central_slices(gt_db, max_slices)
% Seleziona le slice più rappresentative, dando priorità a quelle centrali/equatoriali
% e alle slice con area annotata maggiore.

n = length(gt_db);
if n == 0
    selected = gt_db;   % nessuna slice: restituisce struttura vuota
    return;
end

% Calcola l'area annotata (unione di tutte le maschere) per ogni slice
areas = zeros(1,n);
for i = 1:n
    areas(i) = sum(gt_db(i).gt_union(:));   % numero di pixel annotati
end

is_c  = [gt_db.is_central];    % vettore logico: true per le slice centrali/equatoriali
c_idx = find(is_c);            % indici delle slice centrali
r_idx = find(~is_c);           % indici delle slice non centrali

% Ordina le slice centrali per area decrescente (prima le più grandi)
[~, oc] = sort(areas(c_idx), 'descend');
c_idx = c_idx(oc);

% Ordina le slice non centrali per area decrescente
[~, orr] = sort(areas(r_idx), 'descend');
r_idx = r_idx(orr);

ordered  = [c_idx, r_idx];                             % unisce: prima le centrali poi le altre
keep     = ordered(1:min(max_slices, length(ordered))); % prende al massimo max_slices slice
selected = gt_db(keep);                                 % estrae le slice selezionate
end


function [dice, prec, sens] = compute_metrics(gt, pred_mask)
% Calcola Dice, Precision e Sensitivity tra la maschera GT e la predizione.

gt   = logical(gt(:));        % vettore colonna booleano della GT
pred = logical(pred_mask(:)); % vettore colonna booleano della predizione

tp = sum(gt & pred);    % True Positive: pixel annotati e predetti correttamente
fp = sum(~gt & pred);   % False Positive: pixel predetti ma non annotati
fn = sum(gt & ~pred);   % False Negative: pixel annotati ma non predetti

dice = 2 * tp / (2 * tp + fp + fn + eps);   % Dice = 2TP / (2TP + FP + FN), eps evita divisione per zero
prec = tp / (tp + fp + eps);                % Precision = TP / (TP + FP)
sens = tp / (tp + fn + eps);                % Sensitivity (Recall) = TP / (TP + FN)
end


function save_debug_image(debug_dir, patient_id, slice_idx, img_hu, gt_best, pred_mask, dice)
% Salva una figura diagnostica a 3 pannelli: slice originale, GT, predizione.

fig      = figure('Visible','off','Position',[0 0 1500 500]);   % figura invisibile 1500×500 pixel
img_disp = mat2gray(img_hu, [-1000, 400]);   % normalizza HU nel range [-1000, 400] → [0, 1] per la visualizzazione

subplot(1,3,1);        % pannello 1: slice originale
imshow(img_disp,[]);
title('Original Slice','FontSize',11);

subplot(1,3,2);        % pannello 2: slice con contorno GT in verde
imshow(img_disp,[]);
hold on;
overlay_contour(gt_best,'g');   % disegna il contorno della GT
title('Best-Match GT (green)','FontSize',11);

subplot(1,3,3);        % pannello 3: slice con predizione in rosso e GT in verde
imshow(img_disp,[]);
hold on;
overlay_contour(pred_mask,'r');   % disegna il contorno della predizione
overlay_contour(gt_best,'g');     % sovrappone anche la GT per confronto
title(sprintf('Prediction (red)  Dice = %.3f', dice),'FontSize',11);

sgtitle(sprintf('%s  |  slice %d', strrep(patient_id,'_','-'), slice_idx), ...
    'FontSize',12,'FontWeight','bold');   % titolo generale della figura

% Salva come PNG nella cartella di debug con nome che include paziente, slice e Dice
fname = fullfile(debug_dir, sprintf('%s_s%02d_dice%.3f.png', patient_id, slice_idx, dice));
exportgraphics(fig, fname, 'Resolution', 100);   % esporta a 100 DPI
close(fig);   % chiude la figura per liberare memoria
end


function overlay_contour(mask, color)
% Disegna il contorno di una maschera binaria sull'asse corrente.

if ~any(mask(:))
    return;   % maschera vuota: niente da disegnare
end

bounds = bwboundaries(mask,'noholes');   % estrae i poligoni dei bordi (ignora i buchi interni)
for b = 1:length(bounds)
    % Disegna il contorno: colonna (x) e riga (y) dei punti del bordo
    plot(bounds{b}(:,2), bounds{b}(:,1), color, 'LineWidth', 1.5);
end
end


function [lung_mask, lobes] = extract_lung_mask_2d(hu_slice, rows, cols)
% (Stessa funzione presente anche in find_nodules_advanced; si applica qui per la valutazione)

lung_mask = false(rows, cols);
lobes     = {};

air = hu_slice < -320;   % soglia aria: tutto sotto -320 HU
cc  = bwconncomp(air, 4);
if cc.NumObjects == 0
    return;
end

stats = regionprops(cc, 'Area');
[~, si] = sort([stats.Area], 'descend');   % ordina per area decrescente

n_lobes = 0;
for k = 1:length(si)
    comp = false(rows, cols);
    comp(cc.PixelIdxList{si(k)}) = true;

    % Scarta le regioni che toccano il bordo (sfondo esterno o lettino TC)
    if any(comp(1,:)) || any(comp(end,:)) || any(comp(:,1)) || any(comp(:,end))
        continue;
    end

    cf = imfill(comp, 'holes');    % riempie bronchi e vasi per avere una regione compatta
    lung_mask    = lung_mask | cf;
    lobes{end+1} = cf;
    n_lobes = n_lobes + 1;

    if n_lobes >= 2   % si accettano al massimo 2 lobi (destro e sinistro)
        break;
    end
end
end


function print_group_stats(name, dice_vals, prec_vals, sens_vals)
% Stampa le statistiche riassuntive (n, mean Dice, median Dice, mean Prec, mean Sens, >0.60, zeros)
% per un sottogruppo della stratificazione.

if isempty(dice_vals)
    fprintf('  %-18s : no cases\n', name);   % nessun caso nel sottogruppo
    return;
end

fprintf(['  %-18s : n=%3d | meanD=%.4f | medD=%.4f | ' ...
         'meanP=%.4f | meanS=%.4f | >0.60=%2d | zeros=%2d\n'], ...
    name, length(dice_vals), mean(dice_vals), median(dice_vals), ...
    mean(prec_vals), mean(sens_vals), ...
    sum(dice_vals > 0.60), ...   % casi con Dice superiore al 60%
    sum(dice_vals == 0));        % casi dove l'algoritmo non ha segmentato nulla
end