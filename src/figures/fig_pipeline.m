%   Author: Angelo Mastrangelo

%  Generates 5 separate figures for each representative patient:
%   (1) Raw CT — Lung Window
%   (2) Lung Mask Extraction
%   (3) Halo Score Map
%   (4) Prediction vs. Ground Truth
%   (5) Zoom — Prediction vs. GT

clear; close all; clc;

BASE_DIR = fullfile('..', 'Dataset', 'manifest-1585167679499', 'LIDC-IDRI');
OUT_DIR  = 'export_figures';
if ~exist(OUT_DIR,'dir'), mkdir(OUT_DIR); end

PATIENTS = {'LIDC-IDRI-0016', 'LIDC-IDRI-0094'};
SUFFIXES = {'good', 'juxta'};

for pi_idx = 1:length(PATIENTS)
    patient_id   = PATIENTS{pi_idx};
    patient_path = fullfile(BASE_DIR, patient_id);
    suffix       = SUFFIXES{pi_idx};

    fprintf('Processing %s...\n', patient_id);

    [xml_file, dcm_dir] = find_series(patient_path);
    if isempty(xml_file), continue; end

    gt_db = parse_lidc_xml(xml_file, 3);
    mask4 = arrayfun(@(x) length(x.reader_masks) == 4, gt_db);
    if ~any(mask4), continue; end
    gt_db = gt_db(mask4);

    dcm_files = dir(fullfile(dcm_dir, '*.dcm'));
    [vol_hu, z_positions, pixel_spacing, slice_thickness] = load_patient_volume(dcm_dir, dcm_files);
    if isempty(vol_hu), continue; end

    [rows, cols, n_vol] = size(vol_hu);
    selected = select_central_slices(gt_db, 1);
    if isempty(selected), continue; end

    gt_z_pos     = selected(1).z_pos;
    reader_masks = selected(1).reader_masks;
    [~, z_idx]   = min(abs(z_positions - gt_z_pos));
    hu_slice     = vol_hu(:,:,z_idx);

    z_off = max(1, round(5.0 / slice_thickness));
    if z_idx - z_off >= 1 && z_idx + z_off <= n_vol
        hu_above = vol_hu(:,:, z_idx - z_off);
        hu_below = vol_hu(:,:, z_idx + z_off);
    else
        hu_above = [];
        hu_below = [];
    end

    % Use largest reader mask as GT display reference
    best_gt = reader_masks{1};
    for r = 2:length(reader_masks)
        if sum(reader_masks{r}(:)) > sum(best_gt(:))
            best_gt = reader_masks{r};
        end
    end

    %% Recompute internal pipeline maps for visualization
    dil_px   = max(1, round(5.0 / pixel_spacing));
    cls_px   = max(1, round(1.0 / pixel_spacing));
    sig_bg   = max(3, round(8.0 / pixel_spacing));
    halo_in  = max(1, round(3.0 / pixel_spacing));
    halo_out = max(2, round(5.0 / pixel_spacing));

    min_area    = pi * (1.5 / pixel_spacing)^2;
    max_area    = pi * (40  / pixel_spacing)^2;
    juxta_cA    = pi * (15  / pixel_spacing)^2;
    area20mm    = pi * (10  / pixel_spacing)^2;
    juxta_sol_A = pi * (25  / pixel_spacing)^2;

    [lung_mask, lobes] = extract_lung_mask_2d_local(hu_slice, rows, cols);
    search_mask = imdilate(lung_mask, strel('disk', dil_px));
    local_bg    = imgaussfilt(hu_slice, double(sig_bg));

    solid_c = (hu_slice >= -100) & (hu_slice <= 200) & search_mask;
    calc_c  = (hu_slice > 200)   & (hu_slice <= 1500) & lung_mask;

    ggo_c = false(rows, cols);
    gz = (hu_slice >= -800) & (hu_slice < -100) & lung_mask;
    if any(gz(:))
        ggo_c = gz & (hu_slice - local_bg >= 25);
    end

    juxta_c = false(rows, cols);
    per_lobe_dent = false(rows, cols);

    for li = 1:length(lobes)
        lobe = lobes{li};
        if ~any(lobe(:)), continue; end

        hull = bwconvhull(lobe);
        dent = hull & ~lobe;
        if ~any(dent(:)), continue; end

        per_lobe_dent = per_lobe_dent | dent;
        cc_d = bwconncomp(dent, 8);
        sd   = regionprops(cc_d, 'Area');

        for k = 1:cc_d.NumObjects
            if sd(k).Area >= min_area && sd(k).Area <= juxta_cA
                juxta_c(cc_d.PixelIdxList{k}) = true;
            end
        end
    end

    if any(lung_mask(:))
        try
            ch = bwconvhull(lung_mask);
            hd = ch & ~lung_mask & ~per_lobe_dent;
            if any(hd(:))
                cc_h = bwconncomp(hd, 8);
                sh   = regionprops(cc_h, 'Area');
                for k = 1:cc_h.NumObjects
                    if sh(k).Area >= min_area && sh(k).Area <= juxta_cA
                        juxta_c(cc_h.PixelIdxList{k}) = true;
                    end
                end
            end
        catch
        end
    end

    juxta_c = juxta_c & (hu_slice >= -100) & (hu_slice <= 200);

    cands = solid_c | calc_c | ggo_c | juxta_c;
    cands = imclose(cands, strel('disk', cls_px));
    cands = imfill(cands, 'holes');

    cc2 = bwconncomp(cands, 8);
    sp  = regionprops(cc2, hu_slice, ...
        'Area', 'Eccentricity', 'Solidity', 'Centroid', 'MeanIntensity');

    lung_dist = bwdist(~lung_mask);
    halo_map  = zeros(rows, cols);
    scores    = zeros(1, cc2.NumObjects);
    keep_it   = false(1, cc2.NumObjects);

    for k = 1:cc2.NumObjects
        s    = sp(k);
        px_k = cc2.PixelIdxList{k};

        if s.Area < min_area || s.Area > max_area || s.MeanIntensity < -850
            continue;
        end

        cy = max(1, min(rows, round(s.Centroid(2))));
        cx = max(1, min(cols, round(s.Centroid(1))));
        is_juxta = (lung_dist(cy,cx) < 15) || any(juxta_c(px_k));
        is_large_ggo = ~is_juxta && s.Area >= area20mm && ...
                       s.MeanIntensity < -100 && s.MeanIntensity > -600;

        if is_juxta
            ecc_thr = 0.98;
        elseif is_large_ggo
            ecc_thr = 0.95;
        else
            ecc_thr = 0.93;
        end
        if s.Eccentricity > ecc_thr, continue; end

        blob_bin = false(rows, cols);
        blob_bin(px_k) = true;
        ring_mask = imdilate(blob_bin, strel('disk', halo_out)) & ...
                    ~imdilate(blob_bin, strel('disk', halo_in)) & ~cands;

        if any(ring_mask(:))
            rHU = mean(hu_slice(ring_mask));
        else
            rHU = -550;
        end

        hs = max(0.01, min(1.0, (-rHU - 100.0) / 600.0));
        halo_map(px_k) = hs;

        if ~is_juxta && hs > 0.85 && s.MeanIntensity < -650 && s.Solidity < 0.80
            continue;
        end

        sol_thr = 0.65;
        if is_juxta && s.Area < juxta_sol_A
            sol_thr = 0.10;
        elseif is_large_ggo
            sol_thr = 0.25;
        elseif hs > 0.50
            sol_thr = 0.50;
        end
        if s.Solidity < sol_thr, continue; end

        capped = min(s.Area, pi*(10/pixel_spacing)^2);
        std_hu = std(double(hu_slice(px_k)));
        vb = min(2.0, 1.0 + std_hu / 100.0);
        if s.MeanIntensity > 300
            vb = max(vb, 2.0);
        end

        base_s = capped * (s.Solidity^1.5) * vb;

        if hs > 0.85
            base_s = base_s * 13;
        elseif is_juxta
            base_s = base_s * 15;
        else
            base_s = base_s * max(0.10, hs);
        end

        large_juxta_solid = is_juxta && s.MeanIntensity > -100 && s.Area > area20mm;
        if ~isempty(hu_above) && ~isempty(hu_below) && (s.Area < area20mm || large_juxta_solid)
            sA = sum(hu_above(px_k) > -600) / length(px_k);
            sB = sum(hu_below(px_k) > -600) / length(px_k);
            thr_t = 0.85;
            if large_juxta_solid, thr_t = 0.95; end
            if sA > thr_t && sB > thr_t
                base_s = base_s * 0.10;
            end
        end

        scores(k)  = base_s;
        keep_it(k) = true;
    end

    valid_scores = scores;
    valid_scores(~keep_it) = 0;

    pred_mask = false(rows, cols);
    if any(valid_scores > 0)
        [~, best_k] = max(valid_scores);
        pred_mask(cc2.PixelIdxList{best_k}) = true;
    end

    d_val = compute_dice(best_gt, pred_mask);

    %% Display preparation
    wl = -600;
    ww = 1500;
    wmin = wl - ww/2;
    wmax = wl + ww/2;

    ct_disp = max(0, min(1, (double(hu_slice) - wmin) / (wmax - wmin)));
    ct_rgb  = repmat(ct_disp, 1, 1, 3);

    [rr, cc_roi] = find(lung_mask);
    if ~isempty(rr)
        r1 = max(1, min(rr)-60);
        r2 = min(rows, max(rr)+60);
        c1 = max(1, min(cc_roi)-60);
        c2 = min(cols, max(cc_roi)+60);
    else
        r1 = 1; r2 = rows; c1 = 1; c2 = cols;
    end

    [rg, cg] = find(best_gt | pred_mask);
    if ~isempty(rg)
        zr1 = max(1, min(rg)-80);
        zr2 = min(rows, max(rg)+80);
        zc1 = max(1, min(cg)-80);
        zc2 = min(cols, max(cg)+80);
    else
        zr1 = r1; zr2 = r2; zc1 = c1; zc2 = c2;
    end

    %% ============================================================
    %% STEP 1 — Raw CT
    %% ============================================================
    fig1 = figure('Position',[100 100 900 760], 'Color','w', 'Name',[patient_id ' step1']);
    ax1 = axes('Parent', fig1);
    imshow(ct_rgb(r1:r2, c1:c2, :), 'Parent', ax1);
    title(ax1, sprintf('Raw CT --- Lung Window (%s)', patient_id), ...
        'FontSize', 16, 'FontWeight', 'bold', 'Color', 'k');
    set(ax1, 'XTick', [], 'YTick', []);
    exportgraphics(fig1, fullfile(OUT_DIR, sprintf('fig_pipeline_%s_step1_raw.png', suffix)), 'Resolution', 300);
    close(fig1);

    %% ============================================================
    %% STEP 2 — Lung Mask Extraction
    %% ============================================================
    fig2 = figure('Position',[100 100 900 760], 'Color','w', 'Name',[patient_id ' step2']);
    ax2 = axes('Parent', fig2);

    lm_disp = ct_rgb(r1:r2, c1:c2, :);
    lm_crop = lung_mask(r1:r2, c1:c2);

    lm_disp(:,:,1) = lm_disp(:,:,1) .* (1 - 0.30*lm_crop);
    lm_disp(:,:,2) = lm_disp(:,:,2) .* (1 - 0.30*lm_crop);
    lm_disp(:,:,3) = min(1, lm_disp(:,:,3) + 0.35*lm_crop);

    imshow(lm_disp, 'Parent', ax2); hold(ax2, 'on');
    B_lung = bwboundaries(lm_crop, 8, 'noholes');
    for b = 1:length(B_lung)
        plot(ax2, B_lung{b}(:,2), B_lung{b}(:,1), '-', 'Color', [0 1 1], 'LineWidth', 1.8);
    end
    title(ax2, sprintf('Lung Mask Extraction (%s)', patient_id), ...
        'FontSize', 16, 'FontWeight', 'bold', 'Color', 'k');
    set(ax2, 'XTick', [], 'YTick', []);
    exportgraphics(fig2, fullfile(OUT_DIR, sprintf('fig_pipeline_%s_step2_lungmask.png', suffix)), 'Resolution', 300);
    close(fig2);

    %% ============================================================
    %% STEP 3 — Halo Score Map
    %% ============================================================
    fig3 = figure('Position',[100 100 980 760], 'Color','w', 'Name',[patient_id ' step3']);
    ax3 = axes('Parent', fig3);

    hm_crop   = halo_map(r1:r2, c1:c2);
    cand_crop = cands(r1:r2, c1:c2);
    hm_rgb    = ct_rgb(r1:r2, c1:c2, :);
    cm = parula(256);

    for row = 1:size(hm_crop,1)
        for col = 1:size(hm_crop,2)
            if cand_crop(row,col)
                idxc = max(1, min(256, round(hm_crop(row,col)*255)+1));
                hm_rgb(row,col,:) = 0.35*hm_rgb(row,col,:) + 0.65*reshape(cm(idxc,:),1,1,3);
            end
        end
    end

    imshow(hm_rgb, 'Parent', ax3);
    title(ax3, sprintf('Halo Score Map (%s)', patient_id), ...
        'FontSize', 16, 'FontWeight', 'bold', 'Color', 'k');
    set(ax3, 'XTick', [], 'YTick', []);
    colormap(ax3, parula);
    cb = colorbar(ax3);
    cb.Label.String = 'Halo Score';
    cb.Label.Color  = 'k';
    cb.Color        = 'k';
    clim(ax3, [0 1]);
    exportgraphics(fig3, fullfile(OUT_DIR, sprintf('fig_pipeline_%s_step3_halo.png', suffix)), 'Resolution', 300);
    close(fig3);

    %% ============================================================
    %% STEP 4 — Prediction vs Ground Truth
    %% ============================================================
    fig4 = figure('Position',[100 100 900 760], 'Color','w', 'Name',[patient_id ' step4']);
    ax4 = axes('Parent', fig4);

    ov4 = ct_rgb(r1:r2, c1:c2, :);
    pm4 = pred_mask(r1:r2, c1:c2);
    gt4 = best_gt(r1:r2, c1:c2);

    imshow(ov4, 'Parent', ax4); hold(ax4, 'on');
    Bgt4 = bwboundaries(gt4, 8, 'noholes');
    Bpd4 = bwboundaries(pm4, 8, 'noholes');

    for b = 1:length(Bgt4)
        plot(ax4, Bgt4{b}(:,2), Bgt4{b}(:,1), '-', 'Color', [0 1 0], 'LineWidth', 2.4);
    end
    for b = 1:length(Bpd4)
        plot(ax4, Bpd4{b}(:,2), Bpd4{b}(:,1), '--', 'Color', [1 0 0], 'LineWidth', 2.2);
    end

    text(ax4, 10, 18, sprintf('Dice = %.4f', d_val), ...
        'Color', 'w', 'FontSize', 11, 'FontWeight', 'bold', ...
        'BackgroundColor', [0 0 0], 'Margin', 4);

    title(ax4, sprintf('Prediction vs. Ground Truth (%s)', patient_id), ...
        'FontSize', 16, 'FontWeight', 'bold', 'Color', 'k');
    set(ax4, 'XTick', [], 'YTick', []);
    exportgraphics(fig4, fullfile(OUT_DIR, sprintf('fig_pipeline_%s_step4_pred_vs_gt.png', suffix)), 'Resolution', 300);
    close(fig4);

    %% ============================================================
    %% STEP 5 — Zoom
    %% ============================================================
    fig5 = figure('Position',[100 100 900 760], 'Color','w', 'Name',[patient_id ' step5']);
    ax5 = axes('Parent', fig5);

    ov5 = ct_rgb(zr1:zr2, zc1:zc2, :);
    pm5 = pred_mask(zr1:zr2, zc1:zc2);
    gt5 = best_gt(zr1:zr2, zc1:zc2);

    imshow(ov5, 'Parent', ax5); hold(ax5, 'on');
    Bgt5 = bwboundaries(gt5, 8, 'noholes');
    Bpd5 = bwboundaries(pm5, 8, 'noholes');

    for b = 1:length(Bgt5)
        plot(ax5, Bgt5{b}(:,2), Bgt5{b}(:,1), '-', 'Color', [0 1 0], 'LineWidth', 2.8);
    end
    for b = 1:length(Bpd5)
        plot(ax5, Bpd5{b}(:,2), Bpd5{b}(:,1), '--', 'Color', [1 0 0], 'LineWidth', 2.5);
    end

    title(ax5, sprintf('Zoom --- Prediction vs. GT (%s)', patient_id), ...
        'FontSize', 16, 'FontWeight', 'bold', 'Color', 'k');
    set(ax5, 'XTick', [], 'YTick', []);
    exportgraphics(fig5, fullfile(OUT_DIR, sprintf('fig_pipeline_%s_step5_zoom.png', suffix)), 'Resolution', 300);
    close(fig5);

    fprintf('  Saved 5 separate figures for %s\n', patient_id);
end

fprintf('\nDone.\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOCAL FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function d = compute_dice(gt, pr)
tp = sum(gt(:) & pr(:));
d  = 2 * tp / (sum(gt(:)) + sum(pr(:)) + eps);
end

function [xml_file, dcm_dir] = find_series(p_path)
xml_file = '';
dcm_dir  = '';

d1 = dir(p_path);
d1 = d1([d1.isdir] & ~strcmp({d1.name},'.') & ~strcmp({d1.name},'..'));

for d = 1:length(d1)
    d2 = fullfile(p_path, d1(d).name);
    s1 = dir(d2);
    s1 = s1([s1.isdir] & ~strcmp({s1.name},'.') & ~strcmp({s1.name},'..'));

    for s = 1:length(s1)
        sp2 = fullfile(d2, s1(s).name);
        xl = dir(fullfile(sp2,'*.xml'));
        if ~isempty(xl)
            xml_file = fullfile(sp2, xl(1).name);
            dcm_dir  = sp2;
            return;
        end
    end
end
end

function gt_db = parse_lidc_xml(xml_file, min_pts)
gt_db = struct('z_pos',{},'reader_masks',{});

try
    doc = xmlread(xml_file);
    sessions = doc.getElementsByTagName('readingSession');
    nodule_map = containers.Map('KeyType','char','ValueType','any');

    for si = 0:sessions.getLength-1
        sess = sessions.item(si);
        nods = sess.getElementsByTagName('unblindedReadNodule');

        for ni = 0:nods.getLength-1
            nd = nods.item(ni);
            nid_el = nd.getElementsByTagName('noduleID');
            if nid_el.getLength == 0, continue; end
            nid = strtrim(char(nid_el.item(0).getTextContent()));

            rois = nd.getElementsByTagName('roi');

            for ri = 0:rois.getLength-1
                roi = rois.item(ri);
                z_el   = roi.getElementsByTagName('imageZposition');
                inc_el = roi.getElementsByTagName('inclusion');
                if z_el.getLength == 0 || inc_el.getLength == 0, continue; end

                z_pos = str2double(char(z_el.item(0).getTextContent()));
                incl  = strtrim(char(inc_el.item(0).getTextContent()));
                if ~strcmpi(incl,'TRUE'), continue; end

                edges = roi.getElementsByTagName('edgeMap');
                if edges.getLength < min_pts, continue; end

                xs = zeros(edges.getLength,1);
                ys = zeros(edges.getLength,1);

                for ei = 0:edges.getLength-1
                    em  = edges.item(ei);
                    xel = em.getElementsByTagName('xCoord');
                    yel = em.getElementsByTagName('yCoord');
                    if xel.getLength == 0 || yel.getLength == 0, continue; end
                    xs(ei+1) = str2double(char(xel.item(0).getTextContent()));
                    ys(ei+1) = str2double(char(yel.item(0).getTextContent()));
                end

                key = sprintf('%s_%.4f', nid, z_pos);
                if nodule_map.isKey(key)
                    entry = nodule_map(key);
                else
                    entry = struct('z_pos', z_pos, 'xs_list', {{}}, 'ys_list', {{}});
                end
                entry.xs_list{end+1} = xs;
                entry.ys_list{end+1} = ys;
                nodule_map(key) = entry;
            end
        end
    end

    keys_all = nodule_map.keys;
    slice_map = containers.Map('KeyType','char','ValueType','any');

    for ki = 1:length(keys_all)
        k2  = keys_all{ki};
        ent = nodule_map(k2);
        zkey = sprintf('%.4f', ent.z_pos);

        if slice_map.isKey(zkey)
            sm = slice_map(zkey);
        else
            sm = struct('z_pos', ent.z_pos, 'readers', {{}});
        end

        sm.readers{end+1} = struct('xs_list',{ent.xs_list},'ys_list',{ent.ys_list});
        slice_map(zkey) = sm;
    end

    zkeys = slice_map.keys;
    for zi = 1:length(zkeys)
        sm = slice_map(zkeys{zi});
        if length(sm.readers) < 4, continue; end

        rmasks = {};
        for ri = 1:length(sm.readers)
            rd = sm.readers{ri};
            mask_r = false(512,512);

            for pi2 = 1:length(rd.xs_list)
                xv = rd.xs_list{pi2};
                yv = rd.ys_list{pi2};
                if length(xv) < 3, continue; end
                m2 = poly2mask(xv+1, yv+1, 512, 512);
                mask_r = mask_r | m2;
            end

            if any(mask_r(:))
                rmasks{end+1} = mask_r; %#ok<AGROW>
            end
        end

        if length(rmasks) >= 4
            entry2.z_pos = sm.z_pos;
            entry2.reader_masks = rmasks;
            gt_db(end+1) = entry2; %#ok<AGROW>
        end
    end
catch
end
end

function [vol_hu, z_positions, pixel_spacing, slice_thickness] = load_patient_volume(dcm_dir, dcm_files)
vol_hu = [];
z_positions = [];
pixel_spacing = 0.70;
slice_thickness = 2.0;

if isempty(dcm_files), return; end

try
    info1 = dicominfo(fullfile(dcm_dir, dcm_files(1).name));
    pixel_spacing = info1.PixelSpacing(1);
    n = length(dcm_files);
    z_arr = zeros(1,n);

    for i = 1:n
        inf2 = dicominfo(fullfile(dcm_dir, dcm_files(i).name));
        z_arr(i) = inf2.ImagePositionPatient(3);
    end

    [z_arr, si] = sort(z_arr);
    dcm_files = dcm_files(si);

    tmp = dicomread(fullfile(dcm_dir, dcm_files(1).name));
    [r,c] = size(tmp);
    vol = zeros(r,c,n,'int16');

    for i = 1:n
        vol(:,:,i) = dicomread(fullfile(dcm_dir, dcm_files(i).name));
    end

    slope = 1;
    intercept = -1024;
    if isfield(info1,'RescaleSlope'), slope = info1.RescaleSlope; end
    if isfield(info1,'RescaleIntercept'), intercept = info1.RescaleIntercept; end

    vol_hu = double(vol)*slope + intercept;
    z_positions = z_arr;

    dz = abs(diff(z_arr));
    if ~isempty(dz)
        slice_thickness = max(0.5, median(dz));
    end
catch
end
end

function sel = select_central_slices(gt_db, max_slices)
if length(gt_db) <= max_slices
    sel = gt_db;
    return;
end

for i = 1:length(gt_db)
    gt_db(i).gt_union = false(512,512);
    for r = 1:length(gt_db(i).reader_masks)
        gt_db(i).gt_union = gt_db(i).gt_union | gt_db(i).reader_masks{r};
    end
end

areas = arrayfun(@(x) sum(x.gt_union(:)), gt_db);
[~, idx] = sort(areas, 'descend');
sel = gt_db(idx(1:max_slices));
end

function [lung_mask, lobes] = extract_lung_mask_2d_local(hu_slice, rows, cols)
lung_mask = false(rows, cols);
lobes     = {};

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