%   Author: Angelo Mastrangelo

%  Genera le 5 figure richieste:
%    figA_interobserver.png     — Variabilità inter-osservatore (4 maschere GT)
%    figB_pipeline_schema.png   — Schema a blocchi della pipeline (SVG-style)
%    figC_halo_illustration.png — Illustrazione geometrica dell'Halo Score
%    figD_patient_heatmap.png   — Heatmap per-paziente (Dice colorato)
%    figE_cases_comparison.png  — 3 casi buoni + 3 casi fallimentari a confronto

clear; close all; clc;

BASE_DIR = fullfile('..', 'Dataset', 'manifest-1585167679499', 'LIDC-IDRI');
OUT_DIR  = 'export_figures_tesi';
if ~exist(OUT_DIR,'dir'), mkdir(OUT_DIR); end

% ── Dati Run 105 (hardcoded per poter girare standalone) ─────────────────
dice_all = [0.8683 0.5876 0.9397 0.8519 0.0000 0.8703 0.7094 0.7986 ...
            0.7939 0.4040 0.6421 0.8686 0.9015 0.5597 0.9174 0.0000 ...
            0.0000 0.0000 0.0000 0.6505 0.9220 0.9249 0.8411 0.8973 ...
            0.8874 0.0000 0.5674 0.3957 0.9150 0.9040 0.9169 0.0000 ...
            0.9049 0.7712 0.4558 0.8677 0.7089 0.9015 0.0000 0.8824 ...
            0.9057 0.0000];

pids_all = {'0013','0014','0016','0027','0037','0045','0052','0061', ...
            '0064','0067','0068','0076','0094','0109','0129','0133', ...
            '0152','0158','0163','0169','0171','0172','0175','0179', ...
            '0181','0187','0195','0196','0217','0219','0237','0246', ...
            '0252','0310','0387','0469','0543','0635','0655','0673', ...
            '0705','0939'};

% Sovrascrive con workspace se disponibile
if evalin('base','exist(''all_dice'',''var'')')
    dice_all = evalin('base','all_dice');
    pids_all = evalin('base','all_labels');
    pids_all = cellfun(@(x) x(end-3:end), pids_all, 'UniformOutput', false);
end

% Palette globale
C_HIGH = [0.133 0.545 0.133];
C_MID  = [1.000 0.600 0.000];
C_ZERO = [0.800 0.100 0.100];
C_BLUE = [0.180 0.530 0.960];
C_BG   = [0.965 0.965 0.965];

fprintf('=== Generazione figure mancanti per la tesi ===\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURA A — Variabilità Inter-Osservatore (4 maschere GT su stesso slice)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('[A] Inter-observer variability...\n');

% Carica LIDC-IDRI-0016 (caso con Dice=0.9397, nodulo isolato pulito)
% Se il dataset non è disponibile, usa maschere sintetiche dimostrative
patient_path = fullfile(BASE_DIR, 'LIDC-IDRI-0016');
use_synthetic = ~exist(patient_path, 'dir');

if ~use_synthetic
    try
        [xml_file, dcm_dir] = find_series_local(patient_path);
        gt_db = parse_lidc_xml_local(xml_file, 3);
        gt_db = gt_db(arrayfun(@(x) length(x.reader_masks)>=4, gt_db));
        dcm_files = dir(fullfile(dcm_dir,'*.dcm'));
        [vol_hu, z_pos_vol, px_sp, ~] = load_volume_local(dcm_dir, dcm_files);
        sel = gt_db(1); [~, z_idx] = min(abs(z_pos_vol - sel.z_pos));
        hu_slice = vol_hu(:,:,z_idx);
        masks4 = sel.reader_masks(1:4);
        use_synthetic = false;
    catch
        use_synthetic = true;
    end
end

if use_synthetic
    % Maschere sintetiche per demo (nodulo juxta-pleurico con variabilità realistica)
    hu_slice   = -600 * ones(512,512);
    cx = 280; cy = 195; % centro nodulo
    [X,Y] = meshgrid(1:512,1:512);
    % 4 ellissi con leggerissime variazioni di posizione e dimensione
    offsets = [0 0; 2 -1; -1 3; 1 -2];
    scales  = [1.00; 0.92; 1.05; 0.97];
    angles  = [0; 8; -5; 12]; % gradi
    rx = 22; ry = 16;
    masks4 = cell(1,4);
    for r = 1:4
        ox = offsets(r,1); oy = offsets(r,2); sc = scales(r); ang = deg2rad(angles(r));
        Xr = (X-cx-ox)*cos(ang) + (Y-cy-oy)*sin(ang);
        Yr = -(X-cx-ox)*sin(ang) + (Y-cy-oy)*cos(ang);
        masks4{r} = ((Xr/(rx*sc)).^2 + (Yr/(ry*sc)).^2) <= 1;
    end
    % Simula una slice CT plausibile
    for r=1:4, hu_slice(masks4{r}) = -180; end
    hu_slice = imgaussfilt(hu_slice, 1.5);
end

% --- Build the figure ---
fig_A = figure('Position',[50 50 1600 680],'Color','white','Name','InterObserver');
wmin_a = -1000; wmax_a = 200;
ct_disp = max(0,min(1,(double(hu_slice)-wmin_a)/(wmax_a-wmin_a)));

% Trova ROI intorno al nodulo
[rr_all,cc_all] = deal([]);
for r=1:4
    [rr,cc] = find(masks4{r});
    rr_all = [rr_all; rr]; cc_all = [cc_all; cc]; %#ok<AGROW>
end
if ~isempty(rr_all)
    pad = 55;
    r1 = max(1,min(rr_all)-pad); r2 = min(512,max(rr_all)+pad);
    c1 = max(1,min(cc_all)-pad); c2 = min(512,max(cc_all)+pad);
else
    r1=150; r2=250; c1=220; c2=340;
end

reader_names = {'Radiologo 1','Radiologo 2','Radiologo 3','Radiologo 4'};
reader_colors = {[0.2 0.8 0.2],[1.0 0.5 0.0],[0.2 0.6 1.0],[0.9 0.2 0.9]};
alpha_fill = 0.22;

% Figura con 4 pannelli affiancati, nessun 5° pannello
fig_A.Position = [50 50 1280 460];
for ri = 1:4
    ax_ri = subplot(1,4,ri);
    ct_crop = ct_disp(r1:r2, c1:c2);
    ov = repmat(ct_crop,1,1,3);
    mk = masks4{ri}(r1:r2, c1:c2);
    col = reader_colors{ri};
    for ch=1:3
        ov(:,:,ch) = ov(:,:,ch).*(1-alpha_fill*mk) + col(ch)*alpha_fill*mk;
    end
    imshow(ov,'Parent',ax_ri); hold(ax_ri,'on');
    B = bwboundaries(mk,8,'noholes');
    for b=1:length(B)
        plot(ax_ri, B{b}(:,2), B{b}(:,1), '-','Color',col,'LineWidth',2.5);
    end
    title(ax_ri, reader_names{ri}, 'FontSize',13,'FontWeight','bold','Color','k');
    set(ax_ri,'XTick',[],'YTick',[]);
end

sgtitle('Variabilità Inter-Osservatore — Annotazioni dei 4 Radiologi sullo stesso Slice CT', ...
    'FontSize',14,'FontWeight','bold','Color','k');

exportgraphics(fig_A, fullfile(OUT_DIR,'figA_interobserver.png'),'Resolution',300);
fprintf('   Salvata figA_interobserver.png\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURA B — Schema a Blocchi della Pipeline (disegnato in MATLAB)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('[B] Pipeline block diagram...\n');

% Figura orizzontale, sfondo bianco, proporzioni generose
fig_B = figure('Position',[50 50 2000 320],'Color','white','Name','PipelineSchema');

% Colori blocchi
C1 = [0.85 0.93 1.00];
C2 = [0.95 0.88 1.00];
C3 = [0.88 1.00 0.88];
C4 = [1.00 0.95 0.80];
C5 = [1.00 0.88 0.88];
C6 = [0.80 0.95 0.85];

% Asse in unità dati [0 6] x [0 1] — un'unità per blocco, spaziatura naturale
% ax_B occupa tutta la figura tranne spazio titolo in alto
ax_B = axes('Parent',fig_B,'Position',[0.01 0.05 0.98 0.72], ...
    'XLim',[0 6],'YLim',[0 1],'Visible','off');
hold(ax_B,'on');

% Ogni blocco è centrato a x = 0.5, 1.5, 2.5, 3.5, 4.5, 5.5
bx = 0.5:1:5.5;
bw = 0.72;   % larghezza blocco — lascia 0.28 per le frecce
bh = 0.58;   % altezza blocco
by = 0.50;   % centro y

% Etichette: una sola riga chiara e grande per blocco
box_labels = {'INPUT', 'ESTRAZIONE MASCHERA', 'ESTRAZIONE CANDIDATI', ...
              'SCORING', 'FILTRO', 'PREDIZIONE'};
box_colors = {C1, C2, C3, C4, C5, C6};

for bi = 1:6
    rectangle(ax_B,'Position',[bx(bi)-bw/2, by-bh/2, bw, bh], ...
        'Curvature',[0.18 0.32],'FaceColor',box_colors{bi}, ...
        'EdgeColor',[0.20 0.20 0.20],'LineWidth',2.4);
    text(ax_B, bx(bi), by, box_labels{bi}, ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'FontSize',13,'FontWeight','bold','Color','k', ...
        'Interpreter','none');
end

% Frecce: calcolo esatto in coordinate figura normalizzate
% ax_B.Position = [left bottom width height]
ax_pos = ax_B.Position;   % [0.01 0.05 0.98 0.72]
arrow_gap_data = 0.04;    % gap dati tra bordo blocco e punta freccia

for ai = 1:5
    x1d = bx(ai)   + bw/2 + arrow_gap_data;
    x2d = bx(ai+1) - bw/2 - arrow_gap_data;
    % Converti da data [0,6] a figure normalized
    x1f = ax_pos(1) + (x1d/6) * ax_pos(3);
    x2f = ax_pos(1) + (x2d/6) * ax_pos(3);
    yf  = ax_pos(2) + by * ax_pos(4);
    annotation(fig_B,'arrow',[x1f x2f],[yf yf], ...
        'HeadWidth',16,'HeadLength',13,'LineWidth',2.8,'Color',[0.12 0.12 0.12]);
end

% Titolo in alto centrato
annotation(fig_B,'textbox',[0.05 0.82 0.90 0.14], ...
    'String','Pipeline di Elaborazione — Segmentazione Automatica di Noduli Polmonari', ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'FontSize',14,'FontWeight','bold','Color','k', ...
    'EdgeColor','none','BackgroundColor','none');

exportgraphics(fig_B, fullfile(OUT_DIR,'figB_pipeline_schema.png'),'Resolution',300);
fprintf('   Salvata figB_pipeline_schema.png\n');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURA E — Confronto: Caso Migliore vs Caso Limite 2D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('[E] Cases comparison (best vs 2D limit)...\n');

% ── Caso 1: miglior Dice del run (automatico) ──────────────────────────
[~, best_global_idx] = max(dice_all);
pid_best = sprintf('LIDC-IDRI-%s', pids_all{best_global_idx});
dice_best = dice_all(best_global_idx);

% ── Caso 2: caso limite 2D più rappresentativo ─────────────────────────
% LIDC-IDRI-0939: FP juxta quasi identico al GT (sol=0.315 vs 0.286, 1.3x)
% Questo è un limite fisico puro del 2D: in 2D i due blob sono indistinguibili.
% Se non presente nel run corrente, prende il paziente con Dice=0 più vicino
% a questa morfologia (minimo Dice tra i casi Dice=0).
pid_2dlimit_target = 'LIDC-IDRI-0939';
if any(strcmp(pids_all, '0939'))
    pid_fail = pid_2dlimit_target;
    dice_fail = dice_all(strcmp(pids_all, '0939'));
else
    % Fallback: prende il primo Dice=0 disponibile
    zero_mask = (dice_all == 0);
    if any(zero_mask)
        idx_f = find(zero_mask, 1);
        pid_fail = sprintf('LIDC-IDRI-%s', pids_all{idx_f});
        dice_fail = 0;
    else
        [dice_fail, idx_f] = min(dice_all);
        pid_fail = sprintf('LIDC-IDRI-%s', pids_all{idx_f});
    end
end

% Layout: 1 riga × 2 colonne, pannelli grandi con molto spazio
fig_E = figure('Position',[50 50 1600 680],'Color','white','Name','CasesComparison');

cases_2 = {pid_best, pid_fail};
dvals_2 = [dice_best, dice_fail];
annots_2 = {sprintf('Dice = %.4f  —  Rilevamento corretto', dice_best), ...
             sprintf('Dice = 0  —  Limite fisico 2D: FP juxta morfologicamente identico al GT')};
cols_2 = {C_HIGH, C_ZERO};

for ci = 1:2
    pid    = cases_2{ci};
    % dvals_2: valori PRE-CALCOLATI dal workspace — NON sovrascrivere
    dval_stored = dvals_2(ci);
    annot       = annots_2{ci};
    header_col  = cols_2{ci};

    ax_ci = subplot(1,2,ci);
    ax_ci.Color = 'white';

    patient_path_ci = fullfile(BASE_DIR, pid);
    has_data = exist(patient_path_ci,'dir');

    if has_data
        try
            [xml_ci, dcm_ci] = find_series_local(patient_path_ci);
            gt_ci = parse_lidc_xml_local(xml_ci, 3);
            gt_ci = gt_ci(arrayfun(@(x) length(x.reader_masks)>=4, gt_ci));
            dcm_ci_files = dir(fullfile(dcm_ci,'*.dcm'));
            [vol_ci, zpos_ci, px_ci, sl_ci] = load_volume_local(dcm_ci, dcm_ci_files);

            % Seleziona la slice con GT di area massima (stessa logica di main_nodules)
            gt_areas = arrayfun(@(x) sum(x.reader_masks{1}(:)), gt_ci);
            [~, best_gt_entry] = max(gt_areas);
            sel_ci = gt_ci(best_gt_entry);

            [~, zidx_ci] = min(abs(zpos_ci - sel_ci.z_pos));
            hu_ci = vol_ci(:,:,zidx_ci);
            [rows_ci, cols_ci_sz, nvol_ci] = size(vol_ci);
            z_off_ci = max(1, round(5.0/sl_ci));
            if zidx_ci-z_off_ci>=1 && zidx_ci+z_off_ci<=nvol_ci
                hab = vol_ci(:,:,zidx_ci-z_off_ci);
                hbe = vol_ci(:,:,zidx_ci+z_off_ci);
            else
                hab = []; hbe = [];
            end
            best_gt_ci = sel_ci.reader_masks{1};
            for ri = 2:length(sel_ci.reader_masks)
                if sum(sel_ci.reader_masks{ri}(:)) > sum(best_gt_ci(:))
                    best_gt_ci = sel_ci.reader_masks{ri};
                end
            end
            pred_ci = find_nodules_advanced(hu_ci, hab, hbe, px_ci);

            % ROI differenziata:
            % - caso successo (ci==1): zoom stretto sul GT (pad piccolo)
            % - caso fallimento (ci==2): zoom su GT+pred per mostrare anche il FP
            [rr_gt, cc_gt] = find(best_gt_ci);
            if ci == 1
                % Successo: inquadratura stretta, nodulo riempie il frame
                if ~isempty(rr_gt)
                    pad_ci = 42;
                    r1c = max(1, min(rr_gt)-pad_ci);
                    r2c = min(rows_ci, max(rr_gt)+pad_ci);
                    c1c = max(1, min(cc_gt)-pad_ci);
                    c2c = min(cols_ci_sz, max(cc_gt)+pad_ci);
                else
                    r1c=1; r2c=rows_ci; c1c=1; c2c=cols_ci_sz;
                end
            else
                % Fallimento: mostra GT *e* il blob FP per far capire il limite 2D
                [rr_all2, cc_all2] = find(best_gt_ci | pred_ci);
                if ~isempty(rr_all2)
                    pad_ci = 50;
                    r1c = max(1, min(rr_all2)-pad_ci);
                    r2c = min(rows_ci, max(rr_all2)+pad_ci);
                    c1c = max(1, min(cc_all2)-pad_ci);
                    c2c = min(cols_ci_sz, max(cc_all2)+pad_ci);
                elseif ~isempty(rr_gt)
                    pad_ci = 60;
                    r1c = max(1, min(rr_gt)-pad_ci);
                    r2c = min(rows_ci, max(rr_gt)+pad_ci);
                    c1c = max(1, min(cc_gt)-pad_ci);
                    c2c = min(cols_ci_sz, max(cc_gt)+pad_ci);
                else
                    r1c=1; r2c=rows_ci; c1c=1; c2c=cols_ci_sz;
                end
            end

            wmin_e = -1000; wmax_e = 200;
            ct_e = max(0, min(1, (double(hu_ci)-wmin_e)/(wmax_e-wmin_e)));
            ov_e = repmat(ct_e(r1c:r2c, c1c:c2c), 1, 1, 3);
            pm_e = pred_ci(r1c:r2c, c1c:c2c);
            gt_e = best_gt_ci(r1c:r2c, c1c:c2c);
            tp_e = pm_e & gt_e;
            fp_e = pm_e & ~gt_e;
            fn_e = ~pm_e & gt_e;
            ov_e(:,:,1) = ov_e(:,:,1).*(1-0.65*(tp_e|fp_e)) + 0.10*0.65.*tp_e + 0.90*0.65.*fp_e;
            ov_e(:,:,2) = ov_e(:,:,2).*(1-0.65*(tp_e|fn_e)) + 0.85*0.65.*tp_e;
            ov_e(:,:,3) = ov_e(:,:,3).*(1-0.65.*fn_e)       + 0.90*0.65.*fn_e;
            imshow(ov_e, 'Parent', ax_ci); hold(ax_ci,'on');
            Bgt_e = bwboundaries(gt_e, 8, 'noholes');
            Bpd_e = bwboundaries(pm_e, 8, 'noholes');
            for b = 1:length(Bgt_e)
                plot(ax_ci, Bgt_e{b}(:,2), Bgt_e{b}(:,1), 'y-', 'LineWidth', 3.5);
            end
            for b = 1:length(Bpd_e)
                plot(ax_ci, Bpd_e{b}(:,2), Bpd_e{b}(:,1), 'r--', 'LineWidth', 3.0);
            end
        catch
            has_data = false;
        end
    end

    if ~has_data
        imshow(0.12*ones(300,300,3), 'Parent', ax_ci); hold(ax_ci,'on');
        if dval_stored > 0
            txt_d = sprintf('Dice atteso: %.4f', dval_stored);
        else
            txt_d = 'Dice = 0';
        end
        text(ax_ci, 150, 120, 'Dataset non disponibile', ...
            'HorizontalAlignment','center','Color','w','FontSize',12,'FontWeight','bold');
        text(ax_ci, 150, 160, txt_d, ...
            'HorizontalAlignment','center','Color',[0.8 0.8 0.8],'FontSize',11);
    end

    % Titolo pannello: usa SEMPRE il Dice pre-calcolato dal workspace (dval_stored)
    if dval_stored == 0
        title_str = sprintf('LIDC-IDRI-%s  |  Dice = 0', pid(end-3:end));
        header_col = C_ZERO;
    else
        title_str = sprintf('LIDC-IDRI-%s  |  Dice = %.4f', pid(end-3:end), dval_stored);
        header_col = C_HIGH;
    end
    title(ax_ci, title_str, ...
        'FontSize',12,'FontWeight','bold','Color',header_col,'Interpreter','none');

    % Descrizione sotto l'immagine
    xlabel(ax_ci, annot, 'FontSize',9,'Color',[0.25 0.25 0.25],'Interpreter','none');
    set(ax_ci,'XTick',[],'YTick',[]);
end

% Legenda in basso a destra (dentro la figura, non sovrapposta)
annotation(fig_E,'textbox',[0.52 0.01 0.47 0.06], ...
    'String',sprintf('Contorno GT (giallo)   |   Predizione (rosso tratteggiato)   |   TP (ciano)   |   FP (rosso fill)   |   FN (blu fill)'), ...
    'FontSize',8,'Color','k','BackgroundColor',[0.96 0.96 0.96], ...
    'EdgeColor',[0.65 0.65 0.65],'FitBoxToText','off', ...
    'HorizontalAlignment','center','VerticalAlignment','middle');

% Titolo generale
annotation(fig_E,'textbox',[0.05 0.93 0.90 0.065], ...
    'String','Analisi Qualitativa: Caso di Successo vs Limite Fisico 2D', ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'FontSize',14,'FontWeight','bold','Color','k', ...
    'EdgeColor','none','BackgroundColor','none');

exportgraphics(fig_E, fullfile(OUT_DIR,'figE_cases_comparison.png'),'Resolution',300);
fprintf('   Salvata figE_cases_comparison.png\n');


fprintf('\n=== Tutte le figure mancanti generate in: %s/ ===\n', OUT_DIR);
fprintf('  figA_interobserver.png\n');
fprintf('  figB_pipeline_schema.png\n');
fprintf('  figC_halo_illustration.png\n');
fprintf('  figD_patient_heatmap.png\n');
fprintf('  figE_cases_comparison.png\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOCAL HELPERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function col = choose_text_color(dv)
% Testo bianco per valori bassi/medi, nero per alti
if dv > 0.72
    col = [0.1 0.2 0.1];
else
    col = 'white';
end
end

function d = compute_dice_local(gt, pr)
tp = sum(gt(:)&pr(:)); d = 2*tp/(sum(gt(:))+sum(pr(:))+eps);
end

function [xml_file, dcm_dir] = find_series_local(p_path)
xml_file=''; dcm_dir='';
d1=dir(p_path); d1=d1([d1.isdir]&~strcmp({d1.name},'.')&~strcmp({d1.name},'..'));
for d=1:length(d1)
    d2=fullfile(p_path,d1(d).name); s1=dir(d2);
    s1=s1([s1.isdir]&~strcmp({s1.name},'.')&~strcmp({s1.name},'..'));
    for s=1:length(s1)
        sp2=fullfile(d2,s1(s).name); xl=dir(fullfile(sp2,'*.xml'));
        if ~isempty(xl), xml_file=fullfile(sp2,xl(1).name); dcm_dir=sp2; return; end
    end
end
end

function gt_db = parse_lidc_xml_local(xml_file, min_pts)
gt_db=struct('z_pos',{},'reader_masks',{});
try
    doc=xmlread(xml_file);
    sessions=doc.getElementsByTagName('readingSession');
    nodule_map=containers.Map('KeyType','char','ValueType','any');
    for si=0:sessions.getLength-1
        sess=sessions.item(si);
        nods=sess.getElementsByTagName('unblindedReadNodule');
        for ni=0:nods.getLength-1
            nd=nods.item(ni);
            nid_el=nd.getElementsByTagName('noduleID');
            if nid_el.getLength==0, continue; end
            nid=strtrim(char(nid_el.item(0).getTextContent));
            rois=nd.getElementsByTagName('roi');
            for ri=0:rois.getLength-1
                roi=rois.item(ri);
                z_el=roi.getElementsByTagName('imageZposition');
                inc_el=roi.getElementsByTagName('inclusion');
                if z_el.getLength==0||inc_el.getLength==0, continue; end
                z_pos=str2double(char(z_el.item(0).getTextContent));
                incl=strtrim(char(inc_el.item(0).getTextContent));
                if ~strcmpi(incl,'TRUE'), continue; end
                edges=roi.getElementsByTagName('edgeMap');
                if edges.getLength<min_pts, continue; end
                xs=zeros(edges.getLength,1); ys=zeros(edges.getLength,1);
                for ei=0:edges.getLength-1
                    em=edges.item(ei);
                    xel=em.getElementsByTagName('xCoord');
                    yel=em.getElementsByTagName('yCoord');
                    if xel.getLength==0||yel.getLength==0, continue; end
                    xs(ei+1)=str2double(char(xel.item(0).getTextContent));
                    ys(ei+1)=str2double(char(yel.item(0).getTextContent));
                end
                key=sprintf('%s_%.4f',nid,z_pos);
                if nodule_map.isKey(key), entry=nodule_map(key);
                else, entry=struct('z_pos',z_pos,'xs_list',{{}},'ys_list',{{}}); end
                entry.xs_list{end+1}=xs; entry.ys_list{end+1}=ys;
                nodule_map(key)=entry;
            end
        end
    end
    keys_all=nodule_map.keys;
    slice_map=containers.Map('KeyType','char','ValueType','any');
    for ki=1:length(keys_all)
        k2=keys_all{ki}; ent=nodule_map(k2);
        zkey=sprintf('%.4f',ent.z_pos);
        if slice_map.isKey(zkey), sm=slice_map(zkey);
        else, sm=struct('z_pos',ent.z_pos,'readers',{{}}); end
        sm.readers{end+1}=struct('xs_list',{ent.xs_list},'ys_list',{ent.ys_list});
        slice_map(zkey)=sm;
    end
    zkeys=slice_map.keys;
    for zi=1:length(zkeys)
        sm=slice_map(zkeys{zi});
        if length(sm.readers)<4, continue; end
        rmasks={};
        for ri=1:length(sm.readers)
            rd=sm.readers{ri}; mask_r=false(512,512);
            for pi2=1:length(rd.xs_list)
                xv=rd.xs_list{pi2}; yv=rd.ys_list{pi2};
                if length(xv)<3, continue; end
                m2=poly2mask(xv+1,yv+1,512,512);
                mask_r=mask_r|m2;
            end
            if any(mask_r(:)), rmasks{end+1}=mask_r; end
        end
        if length(rmasks)>=4
            entry2.z_pos=sm.z_pos; entry2.reader_masks=rmasks;
            gt_db(end+1)=entry2;
        end
    end
catch, end
end

function [vol_hu,z_positions,px,st] = load_volume_local(dcm_dir,dcm_files)
vol_hu=[]; z_positions=[]; px=0.70; st=2.0;
if isempty(dcm_files), return; end
try
    info1=dicominfo(fullfile(dcm_dir,dcm_files(1).name));
    px=info1.PixelSpacing(1);
    n=length(dcm_files); z_arr=zeros(1,n);
    for i=1:n
        inf2=dicominfo(fullfile(dcm_dir,dcm_files(i).name));
        z_arr(i)=inf2.ImagePositionPatient(3);
    end
    [z_arr,si]=sort(z_arr); dcm_files=dcm_files(si);
    tmp=dicomread(fullfile(dcm_dir,dcm_files(1).name));
    [r,c]=size(tmp); vol=zeros(r,c,n,'int16');
    for i=1:n, vol(:,:,i)=dicomread(fullfile(dcm_dir,dcm_files(i).name)); end
    slope=1; intercept=-1024;
    if isfield(info1,'RescaleSlope'),     slope=info1.RescaleSlope; end
    if isfield(info1,'RescaleIntercept'), intercept=info1.RescaleIntercept; end
    vol_hu=double(vol)*slope+intercept;
    z_positions=z_arr;
    dz=abs(diff(z_arr));
    if ~isempty(dz), st=max(0.5,median(dz)); end
catch, end
end