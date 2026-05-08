%   Author: Angelo Mastrangelo

% Esegui questo script DOPO main_nodules.m
%
% Richiede in workspace:
%   all_dice
%   all_labels
%
% Output:
%   - 01_waterfall_thesis_final.png
%   - 02_zero_autopsy_barchart.png
%   - 03_success_vs_limit_clean.png

clearvars -except all_dice all_labels
clc; close all;

BASE_DIR = fullfile('..', 'Dataset', 'manifest-1585167679499', 'LIDC-IDRI');
OUT_DIR  = 'export_figures_tesi';
if ~exist(OUT_DIR,'dir'), mkdir(OUT_DIR); end

% Impostazioni globali
set(0, 'DefaultAxesColor', 'w');
set(0, 'DefaultTextColor', 'k');
set(0, 'DefaultAxesXColor', 'k');
set(0, 'DefaultAxesYColor', 'k');

if ~exist('all_dice','var') || ~exist('all_labels','var')
    error('Variabili all_dice / all_labels non trovate nel workspace.');
end

n_tot    = length(all_dice);
n_zeros  = sum(all_dice == 0);
mean_d   = mean(all_dice);
med_d    = median(all_dice);

% Palette
c_green = [0.15 0.70 0.35];
c_red   = [0.85 0.20 0.25];
c_orng  = [0.95 0.60 0.20];
c_blue  = [0.18 0.53 0.96];
c_gray  = [0.40 0.40 0.40];
c_purp  = [0.55 0.35 0.80];
c_teal  = [0.10 0.65 0.65];

fprintf('\n');
fprintf('══════════════════════════════════════════════════════\n');
fprintf('  FIGURE FINALI TESI\n');
fprintf('══════════════════════════════════════════════════════\n');
fprintf('  Casi totali   : %d\n', n_tot);
fprintf('  Zeri          : %d\n', n_zeros);
fprintf('  Mean Dice     : %.4f\n', mean_d);
fprintf('  Median Dice   : %.4f\n', med_d);

%% =========================================================================
%% 1. WATERFALL PLOT
%% =========================================================================
fig1 = figure('Position',[100 100 900 650], 'Color','w', 'Name','Waterfall');
ax1 = axes('Parent',fig1); hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on');
set(ax1, 'GridColor', [0.85 0.85 0.85], 'GridAlpha', 0.6, 'FontSize', 13);

sorted_dice = sort(all_dice, 'descend');
x_vals = 1:n_tot;

area(ax1, x_vals, sorted_dice, ...
    'FaceColor', c_green, 'FaceAlpha', 0.20, 'EdgeColor', 'none');
plot(ax1, x_vals, sorted_dice, 'LineWidth', 4, 'Color', c_green);

xlabel('Pazienti Ordinati dal Migliore al Peggiore', 'FontSize', 15, 'FontWeight', 'bold');
ylabel('Dice Score', 'FontSize', 15, 'FontWeight', 'bold');
title('Curva di Robustezza (Waterfall Plot)', 'FontSize', 18, 'FontWeight', 'bold');
xlim([1 n_tot]);
ylim([0 1.05]);

% Box statistiche in alto a destra, senza zeri
stats_str = {sprintf('\\bfMedia:\\rm %.4f', mean_d), ...
             sprintf('\\bfMediana:\\rm %.4f', med_d)};
annotation('textbox', [0.72 0.73 0.16 0.11], 'String', stats_str, ...
    'BackgroundColor', 'w', 'EdgeColor', 'k', 'LineWidth', 1, ...
    'FontSize', 12, 'Color', 'k', 'FitBoxToText', 'on');

exportgraphics(fig1, fullfile(OUT_DIR, '01_waterfall_thesis_final.png'), 'Resolution', 300);
fprintf('  Waterfall Plot salvato.\n');

%% =========================================================================
%% 2. AUTOPSIA DETTAGLIATA DEGLI ZERI
%% =========================================================================
fprintf('\nInizio autopsia dettagliata dei %d casi con Dice = 0...\n', n_zeros);

zero_idx = find(all_dice == 0);

cat_names = { ...
    'Nessuna predizione / detection mancata', ...
    'FP vessel-like isolato', ...
    'FP juxtapleurico / parete toracica', ...
    'FP ilare / strutturale centrale', ...
    'GT juxtapleurico irregolare', ...
    'GT GGO / margini diffusi', ...
    'Errore di ranking contro FP compatti', ...
    'Altro / ambiguo in 2D'};

cat_counts = zeros(1, numel(cat_names));
case_log = cell(length(zero_idx), 1);

for ii = 1:length(zero_idx)
    idx = zero_idx(ii);
    label_i = all_labels{idx};

    try
        [pid, gt_z_pos] = decode_case_label(label_i);

        patient_path = fullfile(BASE_DIR, pid);
        [xml_file, dcm_dir] = find_series_supporting(patient_path);
        if isempty(xml_file)
            cat_counts(8) = cat_counts(8) + 1;
            case_log{ii} = sprintf('%s -> Altro/Ambiguo (XML mancante)', label_i);
            continue;
        end

        gt_db = parse_lidc_xml_supporting(xml_file, 3);
        if isempty(gt_db)
            cat_counts(8) = cat_counts(8) + 1;
            case_log{ii} = sprintf('%s -> Altro/Ambiguo (GT non parsata)', label_i);
            continue;
        end

        gt_z_all = [gt_db.z_pos];
        if isnan(gt_z_pos)
            best_gt_idx = 1;
        else
            [~, best_gt_idx] = min(abs(gt_z_all - gt_z_pos));
        end
        gt_entry = gt_db(best_gt_idx);

        dcm_files = dir(fullfile(dcm_dir, '*.dcm'));
        [vol_hu, z_pos_vol, px, st] = load_volume_supporting(dcm_dir, dcm_files);
        if isempty(vol_hu)
            cat_counts(8) = cat_counts(8) + 1;
            case_log{ii} = sprintf('%s -> Altro/Ambiguo (volume non caricato)', label_i);
            continue;
        end

        [~, z_idx] = min(abs(z_pos_vol - gt_entry.z_pos));
        hu_slice = vol_hu(:,:,z_idx);

        n_vol = size(vol_hu, 3);
        z_offset_slices = max(1, round(5.0 / st));
        z_above = z_idx - z_offset_slices;
        z_below = z_idx + z_offset_slices;
        if z_above >= 1 && z_below <= n_vol
            hu_above = vol_hu(:,:,z_above);
            hu_below = vol_hu(:,:,z_below);
        else
            hu_above = [];
            hu_below = [];
        end

        rows = size(hu_slice,1);
        cols = size(hu_slice,2);
        gt_union = false(rows, cols);
        for r = 1:length(gt_entry.reader_masks)
            rm = gt_entry.reader_masks{r};
            if size(rm,1) ~= rows || size(rm,2) ~= cols
                rm = imresize(rm, [rows cols], 'nearest') > 0;
            end
            gt_union = gt_union | rm;
        end

        cc_gt = bwconncomp(gt_union, 8);
        if cc_gt.NumObjects == 0
            cat_counts(8) = cat_counts(8) + 1;
            case_log{ii} = sprintf('%s -> Altro/Ambiguo (GT vuota)', label_i);
            continue;
        end

        gt_stats = regionprops(cc_gt, hu_slice, 'Area', 'Solidity', 'Eccentricity', 'Centroid', 'MeanIntensity');
        [~, gt_main_idx] = max([gt_stats.Area]);
        gt_s = gt_stats(gt_main_idx);
        gt_mask_main = false(rows, cols);
        gt_mask_main(cc_gt.PixelIdxList{gt_main_idx}) = true;

        [lung_mask, ~] = extract_lung_mask_supporting(hu_slice, rows, cols);
        lung_dist = bwdist(~lung_mask);

        cy_gt = max(1, min(rows, round(gt_s.Centroid(2))));
        cx_gt = max(1, min(cols, round(gt_s.Centroid(1))));
        gt_dist = lung_dist(cy_gt, cx_gt);

        gt_halo = compute_halo_score_supporting(gt_mask_main, hu_slice);
        gt_is_juxta = gt_dist < 15;
        gt_is_ggo   = gt_s.MeanIntensity < -100;

        pred_mask = find_nodules_advanced(hu_slice, hu_above, hu_below, px);

        if ~any(pred_mask(:))
            if gt_is_ggo
                cat_counts(6) = cat_counts(6) + 1;
                case_log{ii} = sprintf('%s -> GT GGO / margini diffusi', label_i);
            elseif gt_is_juxta && gt_s.Solidity < 0.35
                cat_counts(5) = cat_counts(5) + 1;
                case_log{ii} = sprintf('%s -> GT juxtapleurico irregolare', label_i);
            else
                cat_counts(1) = cat_counts(1) + 1;
                case_log{ii} = sprintf('%s -> Nessuna predizione / detection mancata', label_i);
            end
            continue;
        end

        cc_pr = bwconncomp(pred_mask, 8);
        pr_stats = regionprops(cc_pr, hu_slice, 'Area', 'Solidity', 'Eccentricity', 'Centroid', 'MeanIntensity');
        [~, pr_main_idx] = max([pr_stats.Area]);
        pr_s = pr_stats(pr_main_idx);

        pr_mask_main = false(rows, cols);
        pr_mask_main(cc_pr.PixelIdxList{pr_main_idx}) = true;

        cy_pr = max(1, min(rows, round(pr_s.Centroid(2))));
        cx_pr = max(1, min(cols, round(pr_s.Centroid(1))));
        pr_dist = lung_dist(cy_pr, cx_pr);
        pr_halo = compute_halo_score_supporting(pr_mask_main, hu_slice);

        pr_is_juxta = pr_dist < 15;
        pr_is_isolated = pr_halo > 0.85;

        inter = sum(pr_mask_main(:) & gt_mask_main(:));
        ov_gt = inter / (sum(gt_mask_main(:)) + eps);

        if gt_is_juxta && gt_s.Solidity < 0.35 && ov_gt < 0.20
            cat_counts(5) = cat_counts(5) + 1;
            case_log{ii} = sprintf('%s -> GT juxtapleurico irregolare', label_i);
            continue;
        end

        if gt_is_ggo && (gt_s.Solidity < 0.45 || gt_halo < 0.70) && ov_gt < 0.20
            cat_counts(6) = cat_counts(6) + 1;
            case_log{ii} = sprintf('%s -> GT GGO / margini diffusi', label_i);
            continue;
        end

        if pr_is_isolated && ~pr_is_juxta && pr_s.Eccentricity > 0.88 && pr_s.MeanIntensity < -150
            cat_counts(2) = cat_counts(2) + 1;
            case_log{ii} = sprintf('%s -> FP vessel-like isolato', label_i);
            continue;
        end

        if pr_is_juxta && pr_s.MeanIntensity > -150 && pr_s.Area > (pi * (10 / px)^2)
            cat_counts(3) = cat_counts(3) + 1;
            case_log{ii} = sprintf('%s -> FP juxtapleurico / parete toracica', label_i);
            continue;
        end

        if pr_is_juxta && pr_s.MeanIntensity > -150
            cat_counts(4) = cat_counts(4) + 1;
            case_log{ii} = sprintf('%s -> FP ilare / strutturale centrale', label_i);
            continue;
        end

        if (pr_s.Solidity - gt_s.Solidity > 0.20) || ...
           (pr_s.Area < 0.80 * gt_s.Area && pr_s.Solidity > gt_s.Solidity)
            cat_counts(7) = cat_counts(7) + 1;
            case_log{ii} = sprintf('%s -> Errore di ranking contro FP compatti', label_i);
            continue;
        end

        cat_counts(8) = cat_counts(8) + 1;
        case_log{ii} = sprintf('%s -> Altro / ambiguo in 2D', label_i);

    catch ME
        cat_counts(8) = cat_counts(8) + 1;
        case_log{ii} = sprintf('%s -> Altro/Ambiguo (errore: %s)', label_i, ME.message);
    end
end

valid = cat_counts > 0;
cat_names_plot  = cat_names(valid);
cat_counts_plot = cat_counts(valid);

fprintf('\n');
fprintf('══════════════════════════════════════════════════════\n');
fprintf('  DETTAGLIO AUTOPSIA ZERI\n');
fprintf('══════════════════════════════════════════════════════\n');
for i = 1:length(cat_names)
    if cat_counts(i) > 0
        fprintf('  %-38s : %2d (%.1f%%)\n', ...
            cat_names{i}, cat_counts(i), 100*cat_counts(i)/sum(cat_counts));
    end
end

fprintf('\n  Log caso per caso:\n');
for i = 1:length(case_log)
    fprintf('    %s\n', case_log{i});
end

%% =========================================================================
%% 3. BAR CHART ORIZZONTALE AUTOPSIA
%% =========================================================================
fig2 = figure('Position',[150 120 1100 700], 'Color','w', 'Name','Autopsy Zero Bar Chart');
ax2 = axes('Parent',fig2); hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on');
set(ax2, 'GridColor', [0.88 0.88 0.88], 'GridAlpha', 0.7, 'FontSize', 13);

y = 1:length(cat_counts_plot);

bar_colors = [
    c_gray;
    c_red;
    c_orng;
    c_purp;
    c_blue;
    c_teal;
    [0.30 0.30 0.30];
    [0.60 0.60 0.60]
];
bar_colors = bar_colors(1:length(cat_counts_plot), :);

for i = 1:length(cat_counts_plot)
    barh(ax2, y(i), cat_counts_plot(i), 0.72, ...
        'FaceColor', bar_colors(i,:), 'EdgeColor', 'k', 'LineWidth', 1.2);
    pct = 100 * cat_counts_plot(i) / sum(cat_counts_plot);
    text(cat_counts_plot(i) + 0.10, y(i), sprintf('(%.1f%%)', pct), ...
        'VerticalAlignment', 'middle', 'FontSize', 13, 'FontWeight', 'bold', 'Color', 'k');
end

set(ax2, 'YTick', y, 'YTickLabel', cat_names_plot, 'YDir', 'reverse');
xlabel('Numero di casi', 'FontSize', 15, 'FontWeight', 'bold');
title('Autopsia dettagliata dei fallimenti', 'FontSize', 18, 'FontWeight', 'bold');
xlim([0 max(cat_counts_plot) + 1.8]);

exportgraphics(fig2, fullfile(OUT_DIR, '02_zero_autopsy_barchart.png'), 'Resolution', 300);
fprintf('\n  Bar chart autopsia zeri salvato.\n');

%% =========================================================================
%% 4. CONFRONTO QUALITATIVO: SUCCESSO VS LIMITE 2D
%% =========================================================================
fprintf('\nCostruzione figura qualitativa di confronto...\n');

[~, idx_best] = max(all_dice);
label_best = all_labels{idx_best};
[hu_best, gt_best, pred_best, pid_best] = recover_case_masks(label_best, BASE_DIR);

idx_zero_candidates = find(all_dice == 0);
idx_bad = [];

for kk = 1:length(idx_zero_candidates)
    test_idx = idx_zero_candidates(kk);
    test_label = all_labels{test_idx};
    try
        [hu_tmp, gt_tmp, pred_tmp, ~] = recover_case_masks(test_label, BASE_DIR);
        if any(pred_tmp(:))
            idx_bad = test_idx;
            hu_bad  = hu_tmp;
            gt_bad  = gt_tmp;
            pred_bad = pred_tmp;
            break;
        end
    catch
    end
end

if isempty(idx_bad)
    idx_bad = idx_zero_candidates(1);
    label_bad = all_labels{idx_bad};
    [hu_bad, gt_bad, pred_bad, pid_bad] = recover_case_masks(label_bad, BASE_DIR);
else
    pid_bad = extract_pid_from_label(all_labels{idx_bad});
end

dice_best = all_dice(idx_best);
dice_bad  = all_dice(idx_bad);

fig3 = figure('Position',[100 100 1300 560], 'Color','w', 'Name','Success vs Limit 2D');

subplot(1,2,1);
imshow(mat2gray(hu_best, [-1000 400])); hold on;
overlay_contour_supporting(gt_best,  c_green, 2.8);
overlay_contour_supporting(pred_best, c_red,  2.8);
title(sprintf('%s  |  Dice = %.4f', pid_best, dice_best), ...
    'FontSize', 15, 'FontWeight', 'bold', 'Color', c_green);
text(size(hu_best,2)/2, size(hu_best,1)+18, ...
    'Caso di successo', ...
    'HorizontalAlignment', 'center', 'FontSize', 13, 'Color', 'k');

subplot(1,2,2);
imshow(mat2gray(hu_bad, [-1000 400])); hold on;
overlay_contour_supporting(gt_bad,  c_green, 2.8);
overlay_contour_supporting(pred_bad, c_red,  2.8);
title(sprintf('%s  |  Dice = %.4f', pid_bad, dice_bad), ...
    'FontSize', 15, 'FontWeight', 'bold', 'Color', c_red);
text(size(hu_bad,2)/2, size(hu_bad,1)+18, ...
    'Caso limite: ambiguità morfologica in 2D', ...
    'HorizontalAlignment', 'center', 'FontSize', 13, 'Color', 'k');

% Titolo principale nero
sgtitle('Confronto qualitativo: caso di successo vs limite fisico 2D', ...
    'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');

% Legenda più piccola e più bassa
annotation('textbox', [0.37 0.005 0.26 0.045], ...
    'String', 'GT (verde)  |  Algoritmo (rosso)', ...
    'HorizontalAlignment', 'center', ...
    'BackgroundColor', 'w', 'EdgeColor', [0.7 0.7 0.7], ...
    'FontSize', 10, 'Color', 'k');

exportgraphics(fig3, fullfile(OUT_DIR, '03_success_vs_limit_clean.png'), 'Resolution', 300);
fprintf('  Figura confronto successo-vs-limite salvata.\n');
fprintf('  Figure salvate in: %s/\n', OUT_DIR);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOCAL HELPERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pid, z_pos] = decode_case_label(label_str)
parts = strsplit(label_str, '_');
pid = parts{1};
z_pos = NaN;
end

function pid = extract_pid_from_label(label_str)
parts = strsplit(label_str, '_');
pid = parts{1};
end

function [hu_slice, gt_mask_main, pred_mask, pid] = recover_case_masks(label_str, BASE_DIR)
pid = extract_pid_from_label(label_str);

patient_path = fullfile(BASE_DIR, pid);
[xml_file, dcm_dir] = find_series_supporting(patient_path);
if isempty(xml_file)
    error('XML non trovato per %s', pid);
end

gt_db = parse_lidc_xml_supporting(xml_file, 3);
if isempty(gt_db)
    error('GT non parsata per %s', pid);
end

areas = zeros(1, length(gt_db));
for i = 1:length(gt_db)
    gt_u = false(512,512);
    for r = 1:length(gt_db(i).reader_masks)
        gt_u = gt_u | gt_db(i).reader_masks{r};
    end
    areas(i) = sum(gt_u(:));
end
[~, best_gt_idx] = max(areas);
gt_entry = gt_db(best_gt_idx);

dcm_files = dir(fullfile(dcm_dir, '*.dcm'));
[vol_hu, z_pos_vol, px, st] = load_volume_supporting(dcm_dir, dcm_files);
if isempty(vol_hu)
    error('Volume non caricato per %s', pid);
end

[~, z_idx] = min(abs(z_pos_vol - gt_entry.z_pos));
hu_slice = vol_hu(:,:,z_idx);

n_vol = size(vol_hu,3);
z_offset_slices = max(1, round(5.0 / st));
z_above = z_idx - z_offset_slices;
z_below = z_idx + z_offset_slices;
if z_above >= 1 && z_below <= n_vol
    hu_above = vol_hu(:,:,z_above);
    hu_below = vol_hu(:,:,z_below);
else
    hu_above = [];
    hu_below = [];
end

rows = size(hu_slice,1);
cols = size(hu_slice,2);

gt_union = false(rows, cols);
for r = 1:length(gt_entry.reader_masks)
    rm = gt_entry.reader_masks{r};
    if size(rm,1) ~= rows || size(rm,2) ~= cols
        rm = imresize(rm, [rows cols], 'nearest') > 0;
    end
    gt_union = gt_union | rm;
end

cc_gt = bwconncomp(gt_union, 8);
if cc_gt.NumObjects == 0
    gt_mask_main = false(rows, cols);
else
    gt_stats = regionprops(cc_gt, 'Area');
    [~, mx] = max([gt_stats.Area]);
    gt_mask_main = false(rows, cols);
    gt_mask_main(cc_gt.PixelIdxList{mx}) = true;
end

pred_mask = find_nodules_advanced(hu_slice, hu_above, hu_below, px);
end

function [xml_file, dcm_dir] = find_series_supporting(patient_path)
xml_file = '';
dcm_dir  = '';

date_entries = dir(patient_path);
date_entries = date_entries([date_entries.isdir] & ...
    ~strcmp({date_entries.name},'.') & ~strcmp({date_entries.name},'..'));

for d = 1:length(date_entries)
    date_path = fullfile(patient_path, date_entries(d).name);
    series_entries = dir(date_path);
    series_entries = series_entries([series_entries.isdir] & ...
        ~strcmp({series_entries.name},'.') & ~strcmp({series_entries.name},'..'));
    for s = 1:length(series_entries)
        series_path = fullfile(date_path, series_entries(s).name);
        xml_list = dir(fullfile(series_path,'*.xml'));
        if ~isempty(xml_list)
            xml_file = fullfile(series_path, xml_list(1).name);
            dcm_dir  = series_path;
            return;
        end
    end
end
end

function gt_db = parse_lidc_xml_supporting(xml_file, min_edge_pts)
gt_db = struct('z_pos',{},'reader_masks',{});

try
    doc = xmlread(xml_file);
catch
    return;
end

IMG_SIZE = 512;
sessions = doc.getElementsByTagName('readingSession');
n_sessions = sessions.getLength();

sess_pools = cell(1, n_sessions);
for si = 1:n_sessions
    sess_pools{si} = containers.Map('KeyType','char','ValueType','any');
end

z_pool = containers.Map('KeyType','char','ValueType','double');

for si = 0:n_sessions-1
    session = sessions.item(si);
    nodules = session.getElementsByTagName('unblindedReadNodule');
    pool = sess_pools{si+1};

    for ni = 0:nodules.getLength()-1
        nodule = nodules.item(ni);
        rois = nodule.getElementsByTagName('roi');

        for ri = 0:rois.getLength()-1
            roi = rois.item(ri);

            sop_uid = xml_text_supporting(roi, 'imageSOP_UID');
            if isempty(sop_uid), continue; end

            z_str = xml_text_supporting(roi, 'imageZposition');
            if isempty(z_str), continue; end
            z_pos = str2double(z_str);

            inc_str = xml_text_supporting(roi, 'inclusion');
            if ~strcmpi(strtrim(inc_str), 'TRUE'), continue; end

            edges = roi.getElementsByTagName('edgeMap');
            n_edge = edges.getLength();
            if n_edge < min_edge_pts, continue; end

            xv = zeros(1,n_edge);
            yv = zeros(1,n_edge);
            ok = true;

            for ei = 0:n_edge-1
                em = edges.item(ei);
                xn = em.getElementsByTagName('xCoord');
                yn = em.getElementsByTagName('yCoord');
                if xn.getLength()==0 || yn.getLength()==0
                    ok = false;
                    break;
                end
                xv(ei+1) = str2double(char(xn.item(0).getTextContent()));
                yv(ei+1) = str2double(char(yn.item(0).getTextContent()));
            end
            if ~ok, continue; end

            entry.x = xv;
            entry.y = yv;

            if ~isKey(pool, sop_uid)
                pool(sop_uid) = {};
                if ~isKey(z_pool, sop_uid), z_pool(sop_uid) = z_pos; end
            end

            tmp = pool(sop_uid);
            tmp{end+1} = entry;
            pool(sop_uid) = tmp;
        end
    end

    sess_pools{si+1} = pool;
end

all_uids = {};
for si = 1:n_sessions
    all_uids = union(all_uids, keys(sess_pools{si}));
end

for i = 1:length(all_uids)
    uid = all_uids{i};
    reader_masks = {};

    for si = 1:n_sessions
        pool = sess_pools{si};
        if ~isKey(pool, uid), continue; end
        rois = pool(uid);

        mask = false(IMG_SIZE, IMG_SIZE);
        for j = 1:length(rois)
            r = rois{j};
            try
                m = poly2mask(min(max(r.x,1),IMG_SIZE), min(max(r.y,1),IMG_SIZE), IMG_SIZE, IMG_SIZE);
            catch
                continue;
            end
            mask = mask | m;
        end

        if any(mask(:))
            reader_masks{end+1} = mask; %#ok<AGROW>
        end
    end

    if isempty(reader_masks), continue; end
    if ~isKey(z_pool, uid), continue; end

    k = length(gt_db) + 1;
    gt_db(k).z_pos = z_pool(uid);
    gt_db(k).reader_masks = reader_masks;
end
end

function txt = xml_text_supporting(node, tag)
txt = '';
list = node.getElementsByTagName(tag);
if list.getLength()==0, return; end
txt = strtrim(char(list.item(0).getTextContent()));
end

function [vol_hu, z_positions, px, st] = load_volume_supporting(dcm_dir, dcm_files)
vol_hu = [];
z_positions = [];
px = 0.70;
st = 2.0;

if isempty(dcm_files), return; end

slices = struct('z',{},'img',{},'slope',{},'intercept',{},'px_sp',{});
for i = 1:length(dcm_files)
    fpath = fullfile(dcm_dir, dcm_files(i).name);
    try
        info = dicominfo(fpath);
        img  = double(dicomread(info));
    catch
        continue;
    end

    z = NaN;
    if isfield(info,'ImagePositionPatient') && numel(info.ImagePositionPatient)>=3
        z = double(info.ImagePositionPatient(3));
    elseif isfield(info,'SliceLocation')
        z = double(info.SliceLocation);
    end
    if isnan(z), continue; end

    slope = 1; intercept = 0;
    if isfield(info,'RescaleSlope'), slope = double(info.RescaleSlope); end
    if isfield(info,'RescaleIntercept'), intercept = double(info.RescaleIntercept); end

    px_sp = 0.70;
    if isfield(info,'PixelSpacing') && ~isempty(info.PixelSpacing)
        px_sp = double(info.PixelSpacing(1));
    end

    k = length(slices) + 1;
    slices(k).z = z;
    slices(k).img = img;
    slices(k).slope = slope;
    slices(k).intercept = intercept;
    slices(k).px_sp = px_sp;
end

if isempty(slices), return; end

[z_sorted, si] = sort([slices.z]);
slices = slices(si);

px = slices(1).px_sp;
if length(z_sorted) > 1
    st = max(0.5, median(abs(diff(z_sorted))));
end

[rows, cols] = size(slices(1).img);
n_slices = length(slices);

vol_hu = zeros(rows, cols, n_slices);
z_positions = z_sorted;

for i = 1:n_slices
    img_i = slices(i).img;
    if size(img_i,1) ~= rows || size(img_i,2) ~= cols
        img_i = imresize(img_i, [rows cols], 'nearest');
    end
    vol_hu(:,:,i) = img_i .* slices(i).slope + slices(i).intercept;
end
end

function [lung_mask, lobes] = extract_lung_mask_supporting(hu_slice, rows, cols)
lung_mask = false(rows, cols);
lobes = {};

air = hu_slice < -320;
cc  = bwconncomp(air, 4);
if cc.NumObjects == 0, return; end

stats = regionprops(cc, 'Area');
[~, si] = sort([stats.Area], 'descend');

n_lobes = 0;
for k = 1:length(si)
    comp = false(rows, cols);
    comp(cc.PixelIdxList{si(k)}) = true;

    if any(comp(1,:)) || any(comp(end,:)) || any(comp(:,1)) || any(comp(:,end))
        continue;
    end

    cf = imfill(comp, 'holes');
    lung_mask = lung_mask | cf;
    lobes{end+1} = cf; %#ok<AGROW>
    n_lobes = n_lobes + 1;
    if n_lobes >= 2, break; end
end
end

function halo_score = compute_halo_score_supporting(blob_mask, hu_slice)
se_inner = strel('disk', 3);
se_outer = strel('disk', 5);

ring_mask = imdilate(blob_mask, se_outer) & ~imdilate(blob_mask, se_inner);
ring_mask = ring_mask & ~blob_mask;

if any(ring_mask(:))
    ring_HU = mean(hu_slice(ring_mask));
else
    ring_HU = -550;
end
halo_score = max(0.01, min(1.0, (-ring_HU - 100.0) / 600.0));
end

function overlay_contour_supporting(mask, color, lw)
if ~any(mask(:)), return; end
B = bwboundaries(mask, 'noholes');
for k = 1:length(B)
    plot(B{k}(:,2), B{k}(:,1), 'Color', color, 'LineWidth', lw);
end
end