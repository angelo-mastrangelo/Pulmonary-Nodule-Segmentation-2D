%   Author: Angelo Mastrangelo

%  Saves:
%    fig5_3d_lung_nolegend.png   — 3D isosurface rendering without legend
%    fig5_3d_lung_legend.png     — 3D isosurface rendering with legend
%    fig6_mpr_views.png          — Multi-view orthogonal slices (clean version)

clear; close all; clc;

BASE_DIR = fullfile('..', 'Dataset', 'manifest-1585167679499', 'LIDC-IDRI');
OUT_DIR  = 'export_figures';
if ~exist(OUT_DIR,'dir'), mkdir(OUT_DIR); end

PATIENT_ID = 'LIDC-IDRI-0016';
fprintf('Loading %s...\n', PATIENT_ID);
patient_path = fullfile(BASE_DIR, PATIENT_ID);

[xml_file, dcm_dir] = find_series(patient_path);
if isempty(xml_file), error('No XML for %s', PATIENT_ID); end

gt_db = parse_lidc_xml(xml_file, 3);
mask4 = arrayfun(@(x) length(x.reader_masks) == 4, gt_db);
gt_db = gt_db(mask4);

dcm_files = dir(fullfile(dcm_dir,'*.dcm'));
[vol_hu, z_positions, pixel_spacing, slice_thickness] = load_patient_volume(dcm_dir, dcm_files);
[rows, cols, n_slices] = size(vol_hu);
fprintf('Volume: %dx%dx%d  px=%.3fmm  dz=%.2fmm\n', rows, cols, n_slices, pixel_spacing, slice_thickness);

selected = select_central_slices(gt_db, 1);
gt_z_pos     = selected(1).z_pos;
reader_masks = selected(1).reader_masks;
[~, z_idx]   = min(abs(z_positions - gt_z_pos));
hu_slice     = vol_hu(:,:, z_idx);

% Best GT mask
best_gt = reader_masks{1};
for r = 2:length(reader_masks)
    if sum(reader_masks{r}(:)) > sum(best_gt(:))
        best_gt = reader_masks{r};
    end
end
fprintf('Equatorial slice: z_idx=%d, z=%.1fmm\n', z_idx, gt_z_pos);

%% ── Build 3D lung mask (fast, ~50 slices around equator) ─────────────────
z_range_mm = 80;
z_range_sl = max(1, round(z_range_mm / slice_thickness));
iz1 = max(1, z_idx - z_range_sl);
iz2 = min(n_slices, z_idx + z_range_sl);
fprintf('Building lung mask for slices %d:%d...\n', iz1, iz2);

lung3d = false(rows, cols, iz2-iz1+1);
for iz = iz1:iz2
    sl = vol_hu(:,:,iz);
    [lm, ~] = extract_lung_mask_2d_local(sl, rows, cols);
    lung3d(:,:,iz-iz1+1) = lm;
end

%% ── Build nodule GT 3D mask (±20mm around equator) ───────────────────────
nod_range_sl = max(1, round(20 / slice_thickness));
niz1 = max(iz1, z_idx - nod_range_sl);
niz2 = min(iz2, z_idx + nod_range_sl);
nodule3d = false(rows, cols, iz2-iz1+1);

[rgt, cgt] = find(best_gt);
if ~isempty(rgt)
    r_cen = round(mean(rgt));
    c_cen = round(mean(cgt));
    rng_r = max(1,r_cen-20):min(rows,r_cen+20);
    rng_c = max(1,c_cen-20):min(cols,c_cen+20);
    nod_hu_vals = hu_slice(rng_r, rng_c);
    nod_hu_mu   = mean(nod_hu_vals(:));
    nod_hu_tol  = 150;

    for iz = niz1:niz2
        sl = vol_hu(:,:,iz);
        rr = max(1,r_cen-25):min(rows,r_cen+25);
        rc = max(1,c_cen-25):min(cols,c_cen+25);

        patch_sl = false(rows,cols);
        hu_patch = sl(rr,rc);
        patch_sl(rr,rc) = (hu_patch >= nod_hu_mu - nod_hu_tol) & ...
                          (hu_patch <= nod_hu_mu + nod_hu_tol);
        patch_sl = imfill(patch_sl,'holes');
        nodule3d(:,:,iz-iz1+1) = patch_sl;
    end
end

%% ── Downsample for isosurface rendering ──────────────────────────────────
ds = 2;
lung3d_ds   = lung3d(1:ds:end, 1:ds:end, :);
nodule3d_ds = nodule3d(1:ds:end, 1:ds:end, :);

vs_xy = pixel_spacing * ds;
vs_z  = slice_thickness;

%% ═══════════════════════════════════════════════════════════════════════════
%% FIGURE 5A — 3D Lung Isosurface + Nodule (NO LEGEND)
%% ═══════════════════════════════════════════════════════════════════════════
fig5a = figure('Position',[50 50 1200 700],'Color','black','Name','3D Lung No Legend');

lung_smooth = smooth3(double(lung3d_ds), 'gaussian', 5);
[faces_L, verts_L] = isosurface(lung_smooth, 0.5);
verts_L(:,1) = verts_L(:,1) * vs_xy;
verts_L(:,2) = verts_L(:,2) * vs_xy;
verts_L(:,3) = verts_L(:,3) * vs_z;

p_lung = patch('Faces', faces_L, 'Vertices', verts_L, ...
    'FaceColor', [0.65 0.78 0.90], 'EdgeColor', 'none', ...
    'FaceAlpha', 0.18, 'FaceLighting', 'gouraud');
hold on;

p_nod = [];
if any(nodule3d_ds(:))
    nod_smooth = smooth3(double(nodule3d_ds), 'gaussian', 3);
    [faces_N, verts_N] = isosurface(nod_smooth, 0.5);
    verts_N(:,1) = verts_N(:,1) * vs_xy;
    verts_N(:,2) = verts_N(:,2) * vs_xy;
    verts_N(:,3) = verts_N(:,3) * vs_z;

    p_nod = patch('Faces', faces_N, 'Vertices', verts_N, ...
        'FaceColor', [0.95 0.35 0.15], 'EdgeColor', 'none', ...
        'FaceAlpha', 0.90, 'FaceLighting', 'gouraud');

    nod_cx = mean(verts_N(:,1));
    nod_cy = mean(verts_N(:,2));
    nod_cz = mean(verts_N(:,3));

    text(nod_cx+30, nod_cy+30, nod_cz+30, 'Nodule', ...
        'Color', [1 0.8 0.2], 'FontSize', 13, 'FontWeight', 'bold');
    plot3([nod_cx nod_cx+28], [nod_cy nod_cy+28], [nod_cz nod_cz+28], ...
        '-', 'Color', [1 0.8 0.2], 'LineWidth', 1.5);
end

camlight('headlight');
camlight('left');
lighting gouraud;
axis equal;
axis off;
view([-55 28]);
set(gca,'Color','black');

exportgraphics(fig5a, fullfile(OUT_DIR,'fig5_3d_lung_nolegend.png'), 'Resolution', 220);
fprintf('Saved fig5_3d_lung_nolegend.png\n');

%% ═══════════════════════════════════════════════════════════════════════════
%% FIGURE 5B — 3D Lung Isosurface + Nodule (WITH LEGEND)
%% ═══════════════════════════════════════════════════════════════════════════
fig5b = figure('Position',[50 50 1200 700],'Color','black','Name','3D Lung With Legend');

p_lung2 = patch('Faces', faces_L, 'Vertices', verts_L, ...
    'FaceColor', [0.65 0.78 0.90], 'EdgeColor', 'none', ...
    'FaceAlpha', 0.18, 'FaceLighting', 'gouraud');
hold on;

p_nod2 = [];
if any(nodule3d_ds(:))
    p_nod2 = patch('Faces', faces_N, 'Vertices', verts_N, ...
        'FaceColor', [0.95 0.35 0.15], 'EdgeColor', 'none', ...
        'FaceAlpha', 0.90, 'FaceLighting', 'gouraud');

    nod_cx = mean(verts_N(:,1));
    nod_cy = mean(verts_N(:,2));
    nod_cz = mean(verts_N(:,3));

    text(nod_cx+30, nod_cy+30, nod_cz+30, 'Nodule', ...
        'Color', [1 0.8 0.2], 'FontSize', 13, 'FontWeight', 'bold');
    plot3([nod_cx nod_cx+28], [nod_cy nod_cy+28], [nod_cz nod_cz+28], ...
        '-', 'Color', [1 0.8 0.2], 'LineWidth', 1.5);
end

camlight('headlight');
camlight('left');
lighting gouraud;
axis equal;
axis off;
view([-55 28]);
set(gca,'Color','black');

if ~isempty(p_nod2)
    legend([p_lung2, p_nod2], {'Lung parenchyma', 'Detected nodule'}, ...
        'Location','northeast', 'TextColor','w', 'Color',[0.12 0.12 0.12], ...
        'FontSize', 11);
else
    legend([p_lung2], {'Lung parenchyma'}, ...
        'Location','northeast', 'TextColor','w', 'Color',[0.12 0.12 0.12], ...
        'FontSize', 11);
end

exportgraphics(fig5b, fullfile(OUT_DIR,'fig5_3d_lung_legend.png'), 'Resolution', 220);
fprintf('Saved fig5_3d_lung_legend.png\n');

%% ═══════════════════════════════════════════════════════════════════════════
%% FIGURE 6 — Multi-view Orthogonal Slices (MPR) CLEAN
%% ═══════════════════════════════════════════════════════════════════════════
fig6 = figure('Position',[50 50 1400 520],'Color','white','Name','MPR Views');

wl = -600;
ww = 1500;
wmin = wl - ww/2;
wmax = wl + ww/2;
win = @(x) max(0,min(1,(double(x)-wmin)/(wmax-wmin)));

if ~isempty(rgt)
    cen_r = r_cen;
    cen_c = c_cen;
    cen_z = z_idx;
else
    cen_r = round(rows/2);
    cen_c = round(cols/2);
    cen_z = z_idx;
end

pad_xy = 60;
pad_z  = 8;

rng_r2 = max(1,cen_r-pad_xy):min(rows,cen_r+pad_xy);
rng_c2 = max(1,cen_c-pad_xy):min(cols,cen_c+pad_xy);
rng_z2 = max(1,cen_z-pad_z):min(n_slices,cen_z+pad_z);

axial_hu = vol_hu(rng_r2, rng_c2, cen_z);
axial_gt = best_gt(rng_r2, rng_c2);

sag_hu = squeeze(vol_hu(rng_r2, cen_c, rng_z2));
cor_hu = squeeze(vol_hu(cen_r, rng_c2, rng_z2));

pred_mask = find_nodules_advanced(hu_slice, [], [], pixel_spacing);
axial_pred = pred_mask(rng_r2, rng_c2);

% ── Axial XY: overlay semplice, pochi colori
ax6a = subplot(1,3,1);
axial_rgb = repmat(win(axial_hu), 1, 1, 3);
imshow(axial_rgb); hold on;

Bgt_a = bwboundaries(axial_gt,8,'noholes');
Bpd_a = bwboundaries(axial_pred,8,'noholes');

for b = 1:length(Bgt_a)
    plot(Bgt_a{b}(:,2), Bgt_a{b}(:,1), 'y-', 'LineWidth', 2.2);
end
for b = 1:length(Bpd_a)
    plot(Bpd_a{b}(:,2), Bpd_a{b}(:,1), 'r--', 'LineWidth', 2.0);
end

plot([size(axial_hu,2)/2 size(axial_hu,2)/2], [1 size(axial_hu,1)], ...
    '--', 'Color', [0.8 0.8 0.8], 'LineWidth', 0.8);
plot([1 size(axial_hu,2)], [size(axial_hu,1)/2 size(axial_hu,1)/2], ...
    '--', 'Color', [0.8 0.8 0.8], 'LineWidth', 0.8);

title('Axial (XY)', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
xlabel(sprintf('GT (yellow) | Pred (red-dash) | Dice = %.3f', compute_dice(best_gt,pred_mask)), ...
    'FontSize', 9, 'Color', 'k');
set(ax6a,'XTick',[],'YTick',[]);

% ── Sagittal XZ
ax6b = subplot(1,3,2);
sag_disp = win(sag_hu);
imshow(sag_disp, 'XData', [0 size(sag_hu,2)*slice_thickness], ...
                 'YData', [0 size(sag_hu,1)*pixel_spacing]);
hold on;
[~, eq_z_local] = min(abs(rng_z2 - cen_z));
gt_z_pix = (eq_z_local - 0.5) * slice_thickness;
yline_max = size(sag_hu,1) * pixel_spacing;
plot([gt_z_pix gt_z_pix], [0 yline_max], 'y-', 'LineWidth', 1.5);
set(ax6b, 'YDir', 'normal');
title('Sagittal (XZ)', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
xlabel('z (mm)', 'Color', 'k');
ylabel('row (mm)', 'Color', 'k');
set(ax6b, 'XColor', 'k', 'YColor', 'k');

% ── Coronal YZ
ax6c = subplot(1,3,3);
cor_disp = win(cor_hu);
imshow(cor_disp, 'XData', [0 size(cor_hu,2)*slice_thickness], ...
                 'YData', [0 size(cor_hu,1)*pixel_spacing]);
hold on;
[~, eq_z_loc2] = min(abs(rng_z2 - cen_z));
gt_z_pix2 = (eq_z_loc2 - 0.5) * slice_thickness;
plot([gt_z_pix2 gt_z_pix2], [0 size(cor_hu,1)*pixel_spacing], 'y-', 'LineWidth', 1.5);
set(ax6c, 'YDir', 'normal');
title('Coronal (YZ)', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
xlabel('z (mm)', 'Color', 'k');
ylabel('col (mm)', 'Color', 'k');
set(ax6c, 'XColor', 'k', 'YColor', 'k');

sgtitle('Multi-Planar Reconstruction', 'FontSize', 13, 'FontWeight', 'bold', 'Color', 'k');

exportgraphics(fig6, fullfile(OUT_DIR,'fig6_mpr_views.png'), 'Resolution', 220);
fprintf('Saved fig6_mpr_views.png\n');

close all;
fprintf('\nDone.\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOCAL HELPERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = compute_dice(gt, pr)
tp = sum(gt(:)&pr(:));
d = 2*tp/(sum(gt(:))+sum(pr(:))+eps);
end

function [xml_file, dcm_dir] = find_series(p_path)
xml_file=''; dcm_dir='';
d1=dir(p_path);
d1=d1([d1.isdir]&~strcmp({d1.name},'.')&~strcmp({d1.name},'..'));
for d=1:length(d1)
    d2=fullfile(p_path,d1(d).name);
    s1=dir(d2);
    s1=s1([s1.isdir]&~strcmp({s1.name},'.')&~strcmp({s1.name},'..'));
    for s=1:length(s1)
        sp2=fullfile(d2,s1(s).name);
        xl=dir(fullfile(sp2,'*.xml'));
        if ~isempty(xl)
            xml_file=fullfile(sp2,xl(1).name);
            dcm_dir=sp2;
            return;
        end
    end
end
end

function gt_db = parse_lidc_xml(xml_file, min_pts)
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
                if nodule_map.isKey(key)
                    entry=nodule_map(key);
                else
                    entry=struct('z_pos',z_pos,'xs_list',{{}},'ys_list',{{}});
                end
                entry.xs_list{end+1}=xs;
                entry.ys_list{end+1}=ys;
                nodule_map(key)=entry;
            end
        end
    end
    keys_all=nodule_map.keys;
    slice_map=containers.Map('KeyType','char','ValueType','any');
    for ki=1:length(keys_all)
        k2=keys_all{ki};
        ent=nodule_map(k2);
        zkey=sprintf('%.4f',ent.z_pos);
        if slice_map.isKey(zkey)
            sm=slice_map(zkey);
        else
            sm=struct('z_pos',ent.z_pos,'readers',{{}});
        end
        sm.readers{end+1}=struct('xs_list',{ent.xs_list},'ys_list',{ent.ys_list});
        slice_map(zkey)=sm;
    end
    zkeys=slice_map.keys;
    for zi=1:length(zkeys)
        sm=slice_map(zkeys{zi});
        if length(sm.readers)<4, continue; end
        rmasks={};
        for ri=1:length(sm.readers)
            rd=sm.readers{ri};
            mask_r=false(512,512);
            for pi2=1:length(rd.xs_list)
                xv=rd.xs_list{pi2};
                yv=rd.ys_list{pi2};
                if length(xv)<3, continue; end
                m2=poly2mask(xv+1,yv+1,512,512);
                mask_r=mask_r|m2;
            end
            if any(mask_r(:)), rmasks{end+1}=mask_r; end
        end
        if length(rmasks)>=4
            entry2.z_pos=sm.z_pos;
            entry2.reader_masks=rmasks;
            gt_db(end+1)=entry2;
        end
    end
catch
end
end

function [vol_hu,z_positions,pixel_spacing,slice_thickness] = load_patient_volume(dcm_dir,dcm_files)
vol_hu=[]; z_positions=[]; pixel_spacing=0.70; slice_thickness=2.0;
if isempty(dcm_files), return; end
try
    info1=dicominfo(fullfile(dcm_dir,dcm_files(1).name));
    pixel_spacing=info1.PixelSpacing(1);
    n=length(dcm_files);
    z_arr=zeros(1,n);
    for i=1:n
        inf2=dicominfo(fullfile(dcm_dir,dcm_files(i).name));
        z_arr(i)=inf2.ImagePositionPatient(3);
    end
    [z_arr,si]=sort(z_arr);
    dcm_files=dcm_files(si);
    tmp=dicomread(fullfile(dcm_dir,dcm_files(1).name));
    [r,c]=size(tmp);
    vol=zeros(r,c,n,'int16');
    for i=1:n
        vol(:,:,i)=dicomread(fullfile(dcm_dir,dcm_files(i).name));
    end
    slope=1; intercept=-1024;
    if isfield(info1,'RescaleSlope'), slope=info1.RescaleSlope; end
    if isfield(info1,'RescaleIntercept'), intercept=info1.RescaleIntercept; end
    vol_hu=double(vol)*slope+intercept;
    z_positions=z_arr;
    dz=abs(diff(z_arr));
    if ~isempty(dz), slice_thickness=max(0.5,median(dz)); end
catch
end
end

function sel = select_central_slices(gt_db, max_slices)
if length(gt_db) <= max_slices
    sel = gt_db;
    return;
end
for i=1:length(gt_db)
    gt_db(i).gt_union = false(512,512);
    for r=1:length(gt_db(i).reader_masks)
        gt_db(i).gt_union = gt_db(i).gt_union | gt_db(i).reader_masks{r};
    end
end
areas = arrayfun(@(x) sum(x.gt_union(:)), gt_db);
[~, idx] = sort(areas, 'descend');
sel = gt_db(idx(1:max_slices));
end

function [lung_mask, lobes] = extract_lung_mask_2d_local(hu_slice, rows, cols)
lung_mask = false(rows, cols);
lobes = {};
air = hu_slice < -320;
cc = bwconncomp(air, 4);
if cc.NumObjects == 0, return; end
stats = regionprops(cc, 'Area');
[~, si] = sort([stats.Area], 'descend');
n_lobes = 0;
for k = 1:length(si)
    comp = false(rows, cols);
    comp(cc.PixelIdxList{si(k)}) = true;
    if any(comp(1,:))||any(comp(end,:))||any(comp(:,1))||any(comp(:,end))
        continue;
    end
    cf = imfill(comp, 'holes');
    lung_mask = lung_mask | cf;
    lobes{end+1} = cf;
    n_lobes = n_lobes + 1;
    if n_lobes >= 2, break; end
end
end