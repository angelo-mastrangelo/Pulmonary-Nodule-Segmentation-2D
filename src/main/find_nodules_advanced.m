function pred_mask = find_nodules_advanced(hu_slice, hu_above, hu_below, pixel_spacing)

%% INIZIALIZZAZIONE
[rows, cols] = size(hu_slice);   % legge le dimensioni (righe × colonne) della slice 2D
pred_mask = false(rows, cols);   % inizializza la maschera di output a tutto 'falso' (nessun nodulo)

%% PARAMETRI GEOMETRICI E SOGLIE DIPENDENTI DALLA SCALA

% Raggi degli elementi strutturanti morfologici (espressi in mm)
dil_mm   = 5.0;   % raggio di dilatazione della maschera polmonare per cercare noduli pleurali
close_mm = 1.0;   % raggio di chiusura morfologica per regolarizzare i candidati

% Conversione mm → pixel: divide mm per la risoluzione, prende il massimo con 1
dil_px_xy   = max(1, round(dil_mm   / pixel_spacing));   % dilatazione in pixel
close_px_xy = max(1, round(close_mm / pixel_spacing));   % chiusura in pixel

% Vincoli di area derivati da diametri circolari equivalenti
% Formula: area = π * r² = π * (diametro/2 / pixel_spacing)²
min_area     = pi * (1.5  / pixel_spacing)^2;   % area minima ≈ nodulo da 3 mm di diametro
max_area     = pi * (40.0 / pixel_spacing)^2;   % area massima ≈ nodulo da 80 mm (esclude strutture enormi)
juxta_cand_A = pi * (15.0 / pixel_spacing)^2;   % area massima per i candidati juxta-pleurali (30 mm diam.)
area_10mm    = pi * (5.0  / pixel_spacing)^2;   % area equivalente a un cerchio da 10 mm
area_20mm    = pi * (10.0 / pixel_spacing)^2;   % area equivalente a un cerchio da 20 mm
juxta_sol_A  = pi * (25.0 / pixel_spacing)^2;   % soglia di area per rilassare la solidità nei juxta (50 mm)

% Scala per la stima del background locale (usata per i GGO)
sig_bg = max(3, round(8 / pixel_spacing));   % sigma del filtro gaussiano (≥3 pixel, ≈8 mm)

% Spessore dell'anello contestuale attorno ad ogni candidato (per misurare l'isolamento)
halo_inner_px = max(1, round(3.0 / pixel_spacing));   % raggio interno dell'anello (≈3 mm)
halo_outer_px = max(2, round(5.0 / pixel_spacing));   % raggio esterno dell'anello (≈5 mm)

% Soglie principali di decisione
HALO_THR_ISOLATED = 0.85;   % sopra questa soglia il candidato è considerato "isolato" nell'aria
LUNG_DIST_JUXTA   = 15;     % distanza massima dal bordo polmonare (pixel) per classificare come juxta

%% STEP A — ESTRAZIONE DELLA MASCHERA POLMONARE E DEFINIZIONE DELLA REGIONE DI RICERCA

% Chiama la funzione ausiliaria locale che segmenta il parenchima polmonare
[lung_mask, lobes] = extract_lung_mask_2d(hu_slice, rows, cols);
if ~any(lung_mask(:))   % se nessun pixel polmonare è trovato (slice senza polmone)
    return;             % esce immediatamente senza fare nulla
end

% Creazione degli elementi strutturanti morfologici a forma di disco
se_dil   = strel('disk', dil_px_xy);    % disco per dilatare la maschera polmonare
se_close = strel('disk', close_px_xy);  % disco per la chiusura morfologica dei candidati
se_inner = strel('disk', halo_inner_px); % disco interno per costruire l'anello halo
se_outer = strel('disk', halo_outer_px); % disco esterno per costruire l'anello halo

% Maschera di ricerca allargata: include la zona pleurale (noduli juxtapleurali)
search_mask = imdilate(lung_mask, se_dil);   % dilata la maschera polmonare di ~5 mm

% Stima del background locale con filtro gaussiano (usata per rilevare GGO)
local_bg = imgaussfilt(hu_slice, double(sig_bg));   % immagine sfumata = background locale

% Flag per il filtro 2.5D: attivo solo se entrambe le slice adiacenti esistono
use_tube_filter = ~isempty(hu_above) && ~isempty(hu_below);

%% STEP B — GENERAZIONE MULTI-CLASSE DEI CANDIDATI

% CLASSE 1 – Noduli solidi / parzialmente solidi
% Seleziona i pixel con HU tra -100 e +200 (tessuto molle) nella zona di ricerca
solid_c = (hu_slice >= -100) & (hu_slice <= 200) & search_mask;

% CLASSE 2 – Lesioni calcificate
% Seleziona i pixel con HU tra 200 e 1500 (calcio) solo dentro il polmone
calc_c = (hu_slice > 200) & (hu_slice <= 1500) & lung_mask;

% CLASSE 3 – GGO (Ground-Glass Opacity): regioni a vetro smerigliato
ggo_c    = false(rows, cols);   % inizializza la mappa GGO vuota
% Zona di attenzione GGO: pixel tra -800 e -100 HU dentro il polmone (più densi dell'aria pura)
ggo_zone = (hu_slice >= -800) & (hu_slice < -100) & lung_mask;
if any(ggo_zone(:))   % se esiste almeno un pixel nella zona GGO
    % Considera GGO solo dove il valore HU supera il background locale di almeno 25 HU
    ggo_c = ggo_zone & (hu_slice - local_bg >= 25);
end

% CLASSE 4 – Noduli juxta-pleurali (adesi alla pleura)
% Vengono estratti come "ammaccature" tra il contorno del polmone e il suo guscio convesso
juxta_c       = false(rows, cols);    % mappa candidati juxta
per_lobe_dent = false(rows, cols);    % mappa delle ammaccature per lobo (usata anche per il globale)

for li = 1:length(lobes)       % ciclo su ogni lobo polmonare (massimo 2: dx e sx)
    lobe = lobes{li};           % maschera binaria del lobo corrente
    if ~any(lobe(:))            % salta se il lobo è vuoto
        continue;
    end

    hull = bwconvhull(lobe);    % calcola il guscio convesso del lobo
    dent = hull & ~lobe;        % l'ammaccatura è la parte del guscio non coperta dal lobo
    if ~any(dent(:))            % salta se non ci sono ammaccature
        continue;
    end

    per_lobe_dent = per_lobe_dent | dent;   % accumula le ammaccature di tutti i lobi

    cc_d = bwconncomp(dent, 8);             % analisi delle componenti connesse dell'ammaccatura (8-connettività)
    sd   = regionprops(cc_d, 'Area');       % misura l'area di ogni componente

    for k = 1:cc_d.NumObjects              % ciclo su ogni ammaccatura
        % Accetta l'ammaccatura solo se la sua area è nel range plausibile per un nodulo juxta
        if sd(k).Area >= min_area && sd(k).Area <= juxta_cand_A
            juxta_c(cc_d.PixelIdxList{k}) = true;   % segna l'ammaccatura come candidato juxta
        end
    end
end

% Ammaccature globali (ilari): differenza tra il guscio convesso del polmone intero
% e le ammaccature per-lobo già trovate
if any(lung_mask(:))    % solo se esiste una maschera polmonare
    try
        combined_hull = bwconvhull(lung_mask);         % guscio convesso dell'intero polmone
        hilar_dent    = combined_hull & ~lung_mask & ~per_lobe_dent; % ammaccature non ancora classificate

        if any(hilar_dent(:))
            cc_h = bwconncomp(hilar_dent, 8);          % analisi componenti connesse
            sh   = regionprops(cc_h, 'Area');           % misura le aree

            for k = 1:cc_h.NumObjects
                if sh(k).Area >= min_area && sh(k).Area <= juxta_cand_A
                    juxta_c(cc_h.PixelIdxList{k}) = true;   % aggiunge come candidato juxta
                end
            end
        end
    catch
        % bwconvhull può fallire su maschere degeneri; in quel caso si ignora
        % questa fase senza interrompere il pipeline
    end
end

% Restringe i candidati juxta alla finestra HU del tessuto molle
juxta_c = juxta_c & (hu_slice >= -100) & (hu_slice <= 200);

% Mappa globale dei candidati: unione delle 4 classi
cands = solid_c | calc_c | ggo_c | juxta_c;
if ~any(cands(:))   % nessun candidato trovato: esce
    return;
end

% Regolarizzazione morfologica: chiude i piccoli buchi e riempie le cavità interne
cands = imclose(cands, se_close);    % chiusura morfologica (connette pixel vicini)
cands = imfill(cands, 'holes');      % riempie i buchi interni a ogni regione

%% STEP C — FILTRAGGIO, ANALISI CONTESTUALE E SCORING

% Analisi delle componenti connesse della mappa candidati finale
cc2 = bwconncomp(cands, 8);    % ogni componente = una regione candidata separata
if cc2.NumObjects == 0          % se non ci sono componenti, esce
    return;
end

% Calcola le proprietà radiologiche e geometriche di ogni componente
sp = regionprops(cc2, hu_slice, ...
    'Area', 'Eccentricity', 'Solidity', 'Centroid', 'MeanIntensity');
%   Area          = numero di pixel
%   Eccentricity  = quanto è allungata (0=cerchio, 1=linea)
%   Solidity      = rapporto area/guscio convesso (1=molto compatto)
%   Centroid      = coordinata (x,y) del centro geometrico
%   MeanIntensity = valore HU medio dei pixel nella regione

% Mappa delle distanze dal bordo del polmone (usata per classificare juxta vs centrale)
lung_dist = bwdist(~lung_mask);   % ogni pixel → distanza (pixel) dal bordo polmonare più vicino

% Array di supporto per lo scoring e la selezione finale
scores          = zeros(1, cc2.NumObjects);   % score Tier-1 per ogni componente
keep_it         = false(1, cc2.NumObjects);   % flag di accettazione Tier-1
scores_fallback = zeros(1, cc2.NumObjects);   % score Tier-2 (fallback per candidati irregolari)
keep_fallback   = false(1, cc2.NumObjects);   % flag di accettazione Tier-2

for k = 1:cc2.NumObjects   % ciclo su ogni componente candidata
    s    = sp(k);                   % proprietà della k-esima componente
    px_k = cc2.PixelIdxList{k};    % indici lineari dei pixel della componente

    % ── C.1 GATE DI ESCLUSIONE DURE ─────────────────────────────────────────
    % Elimina subito i componenti troppo piccoli, troppo grandi o con HU impossibile

    if s.Area < min_area || s.Area > max_area
        continue;   % area fuori range plausibile per un nodulo: scarta
    end

    if s.MeanIntensity < -850
        continue;   % HU media tipica dell'aria pura: non può essere un nodulo, scarta
    end

    % ── C.2 CLASSIFICAZIONE DEL CONTESTO ANATOMICO ──────────────────────────
    % Determina se il candidato è juxta-pleurale

    cy = max(1, min(rows, round(s.Centroid(2))));   % coordinata riga del centroide
    cx = max(1, min(cols, round(s.Centroid(1))));   % coordinata colonna del centroide

    % Distanza del centroide dal bordo polmonare (in pixel)
    dist_to_wall = lung_dist(cy, cx);

    % Juxta: centroide vicino al bordo o area piccola vicina al bordo
    is_juxta = dist_to_wall < LUNG_DIST_JUXTA;

    % GGO grande: area significativa, HU nella fascia GGO (-600 < HU < -100)
    is_large_ggo = s.Area >= area_20mm && ...
                   s.MeanIntensity < -100 && ...
                   s.MeanIntensity > -600;

    % ── C.3 FILTRO DI ECCENTRICITÀ ADATTIVO ─────────────────────────────────
    % Soglia di eccentricità più permissiva per strutture che tendono ad essere allungate

    if is_juxta
        ecc_thr = 0.98;      % juxta molto allungati sono accettabili
    elseif is_large_ggo
        ecc_thr = 0.95;      % GGO grandi possono essere ovali
    else
        ecc_thr = 0.93;      % candidati standard: forma quasi rotonda obbligatoria
    end

    if s.Eccentricity > ecc_thr
        continue;   % troppo allungato → probabilmente un vaso o struttura tubulare: scarta
    end

    % ── C.4 STIMA DELLO SCORE DI HALO (ANELLO CONTESTUALE) ─────────────────
    % Misura quanta aria circonda il candidato: score alto = nodulo isolato nel parenchima

    blob_bin = false(rows, cols);   % maschera binaria temporanea solo per il k-esimo candidato
    blob_bin(px_k) = true;          % accende i pixel del candidato

    % Costruisce l'anello: dilata con disco esterno, rimuove la dilatazione interna e i candidati stessi
    ring_mask = imdilate(blob_bin, se_outer) & ...   % zona esterna
                ~imdilate(blob_bin, se_inner) & ...   % rimuove la zona interna
                ~cands;                               % esclude altri candidati adiacenti

    if any(ring_mask(:))
        ring_HU = mean(hu_slice(ring_mask));   % HU medio dell'anello contestuale
    else
        ring_HU = -550;   % valore di default se l'anello è vuoto (≈ aria moderata)
    end

    % Normalizza il valore HU dell'anello in [0,1]:
    % ring_HU molto negativo (aria pura) → halo_score vicino a 1
    % ring_HU positivo (tessuto molle) → halo_score vicino a 0
    halo_score = max(0.01, min(1.0, (-ring_HU - 100.0) / 600.0));

    % ── C.5 SOPPRESSIONE DELL'ENFISEMA ──────────────────────────────────────
    % Regioni molto isolate, a bassa densità e poco compatte sono probabilmente enfisema

    if ~is_juxta && halo_score > 0.85 && ...
       s.MeanIntensity < -650 && s.Solidity < 0.80
        continue;   % classificato come bolla enfisematosa: scarta
    end

    % ── C.6 FILTRI PER BLOB PICCOLI TROPPO "PERFETTI" ───────────────────────
    % Strutture piccole, rotonde, molto solide ma scarsamente isolate sono probabilmente vasi

    % Caso 1: piccolo, solido, rotondo ma poco circondato da aria
    if s.Area < area_10mm && ...
       s.Solidity > 0.80 && ...
       s.Eccentricity < 0.35 && ...
       halo_score < 0.35 && ...
       ~is_juxta
        continue;   % falso positivo vascolare: scarta
    end

    % Caso 2: qualunque dimensione, molto solido, molto rotondo, poco isolato
    if s.Solidity > 0.92 && ...
       s.Eccentricity < 0.35 && ...
       halo_score < 0.30 && ...
       ~is_juxta
        continue;   % probabile vaso sezione trasversale: scarta
    end

    % ── C.7 CALCOLO DELLO SCORE BASE ────────────────────────────────────────
    % Score = contributo di area × morfologia × variabilità HU × peso contestuale

    % Limita il contributo dell'area a quello di un cerchio da 20 mm
    % (evita che noduli enormi dominino sempre il punteggio)
    capped_area = min(s.Area, pi * (10 / pixel_spacing)^2);

    std_hu    = std(double(hu_slice(px_k)));          % deviazione standard HU dentro il candidato
    var_bonus = min(2.0, 1.0 + std_hu / 100.0);      % bonus per variabilità HU (eterogeneo = più nodulo)

    if s.MeanIntensity > 300
        var_bonus = max(var_bonus, 2.0);   % calcificazione: garantisce bonus massimo per var_bonus
    end

    % Formula dello score base
    base_score = capped_area * (s.Solidity^1.5) * var_bonus;
    % Solidity^1.5: penalizza le strutture irregolari più di un esponente lineare

    % Moltiplicatore contestuale in base all'isolamento/posizione
    if halo_score > HALO_THR_ISOLATED
        base_score = base_score * 13;    % molto isolato: amplifica fortemente (quasi certamente nodulo)
    elseif is_juxta
        base_score = base_score * 15;   % juxta-pleurale: massima amplificazione (classe prioritaria)
    else
        base_score = base_score * max(0.10, halo_score);   % candidato standard: peso proporzionale all'isolamento
    end

    % ── C.8 PENALITÀ 2.5D PER STRUTTURE TUBULARI (VASI) ────────────────────
    % I vasi mantengono la loro sezione trasversale nelle slice adiacenti;
    % i noduli sferici invece appaiono più piccoli spostandosi dall'equatore

    large_juxta_solid = is_juxta && ...
                        s.MeanIntensity > -100 && ...
                        s.Area > area_20mm;   % juxta grande e solido = probabile pleura/parete toracica

    if use_tube_filter && (s.Area < area_20mm || large_juxta_solid)
        % Frazione di pixel del candidato che risulta "non-aria" nella slice sopra e sotto
        sA = sum(hu_above(px_k) > -600) / length(px_k);   % persistenza nella slice superiore
        sB = sum(hu_below(px_k) > -600) / length(px_k);   % persistenza nella slice inferiore

        thr_tube = 0.85;          % soglia standard: >85% di pixel persistenti = probabile vaso
        if large_juxta_solid
            thr_tube = 0.93;      % soglia più severa per i juxta grandi
        end

        if sA > thr_tube && sB > thr_tube
            base_score = base_score * 0.10;   % penalizza fortemente (-90%) se la struttura è tubulare
        end
    end

    % ── C.9 STORAGE NEL TIER FALLBACK ───────────────────────────────────────
    % Candidati irregolari ma plausibili vengono tenuti in riserva

    if s.Solidity >= 0.40 && halo_score > 0.50
        keep_fallback(k)   = true;          % accetta nel tier-2 (fallback)
        scores_fallback(k) = base_score;    % salva il suo score
    end

    % ── C.10 GATE ADATTIVO DI SOLIDITÀ ──────────────────────────────────────
    % Soglia di solidità variabile: più permissiva per juxta e GGO

    sol_thr = 0.65;   % soglia standard (forma abbastanza compatta richiesta)

    if is_juxta && s.Area < juxta_sol_A
        sol_thr = 0.10;   % juxta piccoli: ammessa qualsiasi forma (adesi alla pleura = deformati)
    elseif is_large_ggo
        sol_thr = 0.25;   % GGO grandi: forma irregolare accettabile
    elseif halo_score > 0.50
        sol_thr = 0.50;   % candidati moderatamente isolati: soglia intermedia
    end

    if s.Solidity < sol_thr
        continue;   % troppo irregolare per il tier-1: scarta
    end

    % Accettazione Tier-1
    keep_it(k) = true;    % segna il candidato come valido nel tier principale
    scores(k)  = base_score;   % salva il suo score
end


%% STEP D — SELEZIONE DEL VINCITORE

% Azzeramento degli score dei candidati non accettati (così max() non li seleziona)
valid_scores = scores;
valid_scores(~keep_it) = 0;

if any(valid_scores > 0)
    [~, best_k] = max(valid_scores);                    % indice del candidato con score massimo (Tier-1)
    pred_mask(cc2.PixelIdxList{best_k}) = true;         % accende i pixel del candidato vincitore
    return;                                              % ritorna la maschera finale
end

% Fallback: se il Tier-1 è vuoto usa il miglior candidato del Tier-2
valid_fallback = scores_fallback;
valid_fallback(~keep_fallback) = 0;

if any(valid_fallback > 0)
    [~, best_k] = max(valid_fallback);                  % miglior candidato Tier-2
    pred_mask(cc2.PixelIdxList{best_k}) = true;         % accende i suoi pixel
end

end   % fine della funzione principale


%% FUNZIONE AUSILIARIA — ESTRAZIONE MASCHERA POLMONARE 2D

function [lung_mask, lobes] = extract_lung_mask_2d(hu_slice, rows, cols)

lung_mask = false(rows, cols);   % inizializza la maschera a tutto 'falso'
lobes     = {};                  % cell array vuoto per i lobi trovati

% Soglia per l'aria: tutto sotto -320 HU è considerato aria (polmoni + sfondo esterno)
air = hu_slice < -320;

% Analisi delle componenti connesse con 4-connettività (più conservativa)
cc = bwconncomp(air, 4);
if cc.NumObjects == 0   % nessuna regione d'aria trovata (slice anomala)
    return;
end

stats = regionprops(cc, 'Area');              % misura l'area di ogni componente
[~, si] = sort([stats.Area], 'descend');      % ordina per area decrescente (più grandi prima)

n_lobes = 0;   % contatore dei lobi accettati

for k = 1:length(si)
    comp = false(rows, cols);                 % maschera temporanea
    comp(cc.PixelIdxList{si(k)}) = true;      % accende i pixel del k-esimo componente

    % Scarta le regioni che toccano il bordo dell'immagine
    % (aria esterna al paziente, lettino TC, sfondo)
    if any(comp(1,:)) || any(comp(end,:)) || any(comp(:,1)) || any(comp(:,end))
        continue;   % componente a bordo: non è un polmone, salta
    end

    % Riempie i buchi interni (vasi, bronchi) per ottenere una regione compatta
    cf = imfill(comp, 'holes');

    lung_mask = lung_mask | cf;     % aggiunge questo lobo alla maschera globale
    lobes{end+1} = cf;              % memorizza il lobo nel cell array

    n_lobes = n_lobes + 1;
    if n_lobes >= 2                 % prende al massimo 2 lobi (polmone dx e sx)
        break;
    end
end

end