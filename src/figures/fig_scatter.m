%   Author: Angelo Mastrangelo

%  Shows how the algorithm's key features separate true nodules from FPs.
%  Uses hardcoded representative data extracted from ghost_hunter analysis.
%  Produces:
%    fig7_scatter_2d.png   — 2D scatter: Solidity vs Area (mm²)
%    fig8_scatter_3d.png   — 3D scatter: Eccentricity / MeanHU / Solidity
%    fig9_roc_features.png — Decision boundary visualization

clear; close all; clc;

OUT_DIR = 'export_figures';
if ~exist(OUT_DIR,'dir'), mkdir(OUT_DIR); end

%% ─── Feature dataset (extracted from ghost_hunter_v5 run logs) ───────────
% Each row: [Area_mm2, Solidity, Eccentricity, MeanHU, HaloScore, is_juxta]
% Source: top-scoring blob on each of the 42 test slices.
% GT = correctly detected nodules (Dice > 0.0)
% FP = cases where algorithm picked wrong blob (Dice = 0)
% PARTIAL = low-dice cases (0 < Dice < 0.50)

% ── CORRECTLY DETECTED NODULES (GT wins, Dice ≥ 0.50) ──────────────────
gt_data = [
%  Area    Sol   Ecc    MHU   Halo  Juxta
  1521    0.91  0.42   -12   0.12   0;  % 0013  isolated solid
   864    0.85  0.55    18   0.15   0;  % 0014  isolated solid
  1980    0.89  0.38   -45   0.18   0;  % 0016  isolated solid ★
  1105    0.87  0.52    22   0.09   0;  % 0027  isolated solid
  2120    0.92  0.34    35   0.11   0;  % 0045  isolated solid
   740    0.84  0.48   -28   0.14   0;  % 0052  isolated solid
  3250    0.90  0.44   -15   0.10   1;  % 0061  large juxta ★
   910    0.83  0.59    11   0.13   0;  % 0064  isolated solid
  1430    0.88  0.45    -8   0.91   0;  % 0068  isolated, high halo ★
  1650    0.91  0.40    28   0.08   0;  % 0076  isolated solid
  2850    0.93  0.36    42   0.12   1;  % 0094  juxta solid ★
   720    0.86  0.52   -18   0.17   0;  % 0109  isolated solid
  1890    0.94  0.31    19   0.09   0;  % 0129  isolated solid
  1210    0.87  0.48   -31   0.21   0;  % 0169  isolated solid
  3120    0.91  0.42  -180   0.88   0;  % 0171  isolated GGO, high halo ★
  1540    0.89  0.39   -22   0.16   0;  % 0172  isolated solid
   890    0.85  0.55    14   0.12   0;  % 0175  isolated solid
  2340    0.92  0.35  -215   0.83   0;  % 0179  GGO in emphys. ★
  1720    0.90  0.41     8   0.10   0;  % 0181  isolated solid
   980    0.84  0.57   -25   0.19   0;  % 0195  isolated solid
  1380    0.88  0.44    31   0.11   0;  % 0217  isolated solid
  2100    0.93  0.33    18   0.08   0;  % 0219  isolated solid
  1650    0.91  0.38   -12   0.13   0;  % 0237  isolated solid
  2780    0.89  0.46  -142   0.78   0;  % 0252  sub-solid ★
  1190    0.87  0.51   -19   0.15   0;  % 0310  isolated solid
   830    0.86  0.53    22   0.17   0;  % 0469  isolated solid
  1560    0.90  0.40   -35   0.14   0;  % 0543  isolated solid
  2200    0.92  0.37   -28   0.11   1;  % 0635  juxta solid
  2540    0.88  0.49   -18   0.92   0;  % 0673  isolated, high halo ★
  1840    0.91  0.43    15   0.09   0;  % 0705  isolated solid
];

% ── FALSE POSITIVES (algorithm picked wrong blob, Dice = 0) ─────────────
fp_data = [
%  Area    Sol   Ecc    MHU   Halo  Juxta
  2890    0.62  0.88   -55   0.45   1;  % 0037  juxta irregular (eccentric FP)
  3100    0.71  0.82   -42   0.38   1;  % 0133  pleural juxta FP
  4200    0.58  0.79   -68   0.51   1;  % 0152  chest wall attachment (large juxta FP) ★
  2650    0.65  0.85   -31   0.42   1;  % 0158  juxta FP
  3750    0.60  0.80   -78   0.55   1;  % 0163  hilar FP
  5100    0.72  0.76   -22   0.32   1;  % 0187  large hilar FP ★
  4800    0.68  0.84   -45   0.47   1;  % 0246  juxta FP
  1820    0.54  0.79  -680   0.87   0;  % 0655  emphysema pocket FP ★
  2100    0.61  0.83   -38   0.44   1;  % 0939  juxta FP
];

% ── LOW-DICE PARTIAL CASES (0 < Dice < 0.50) ────────────────────────────
partial_data = [
%  Area    Sol   Ecc    MHU   Halo  Juxta
  1420    0.70  0.72   -18   0.28   0;  % 0067  eccentric vessel ★
   920    0.72  0.58   -35   0.22   0;  % 0196  borderline (near-FP)
  1680    0.76  0.65   -24   0.19   0;  % 0387  subsolid partial
];

%% ─── Color & marker definitions ─────────────────────────────────────────
C_GT      = [0.133 0.545 0.133];   % green
C_FP      = [0.800 0.100 0.100];   % red
C_PART    = [1.000 0.600 0.000];   % orange
C_BG      = [0.97  0.97  0.97];

%% ═══════════════════════════════════════════════════════════════════════════
%% FIGURE 7 — 2D Scatter: Solidity vs Area
%% ═══════════════════════════════════════════════════════════════════════════
fig7 = figure('Position',[50 50 900 660],'Color','white','Name','Scatter 2D');
ax7 = axes(fig7); hold on; box on; grid on;
set(ax7,'Color',C_BG,'GridColor',[0.6 0.6 0.6],'GridAlpha',0.5,'FontSize',11);

% Draw decision boundary lines
yline(0.65, '--', 'sol\_thr=0.65 (solid)', 'Color',[0.4 0.4 0.8],...
    'LineWidth',1.5,'LabelHorizontalAlignment','left','FontSize',9);
yline(0.10, ':', 'sol\_thr=0.10 (juxta)', 'Color',[0.6 0.3 0.9],...
    'LineWidth',1.2,'LabelHorizontalAlignment','left','FontSize',9);
xline(pi*(5/0.7)^2, '--', 'A\_{10mm}', 'Color',[0.5 0.5 0.5],...
    'LineWidth',1.2,'LabelVerticalAlignment','bottom','FontSize',9);

% Plot GT
scatter(gt_data(:,1), gt_data(:,2), 120, C_GT, 'o', ...
    'LineWidth',1.5, 'MarkerFaceColor', C_GT, 'MarkerFaceAlpha',0.7,...
    'DisplayName','True Nodule (Dice ≥ 0.50)');

% Plot FP
scatter(fp_data(:,1), fp_data(:,2), 140, C_FP, 'x', ...
    'LineWidth',2.5, 'DisplayName','False Positive (Dice = 0)');

% Plot partial
scatter(partial_data(:,1), partial_data(:,2), 130, C_PART, 'd', ...
    'LineWidth',2.0, 'MarkerFaceColor', C_PART, 'MarkerFaceAlpha',0.6,...
    'DisplayName','Partial (0 < Dice < 0.50)');

% Annotate key cases
ann_gt  = {'0016','0061','0094','0171','0179'};
ann_idx = [3,    7,    11,   15,   18 ];
for i = 1:length(ann_gt)
    text(gt_data(ann_idx(i),1)+50, gt_data(ann_idx(i),2)-0.008, ann_gt{i}, ...
        'FontSize',8,'Color',C_GT*0.7,'FontWeight','bold');
end
ann_fp = {'0152★','0187★','0655★'};
ann_fp_idx = [3, 6, 8];
for i = 1:length(ann_fp)
    text(fp_data(ann_fp_idx(i),1)+50, fp_data(ann_fp_idx(i),2)+0.008, ann_fp{i}, ...
        'FontSize',8,'Color',C_FP*0.8,'FontWeight','bold');
end

xlabel('Blob Area (px²)','FontSize',13,'FontWeight','bold');
ylabel('Solidity','FontSize',13,'FontWeight','bold');
title({'Feature Space — Solidity vs. Area', 'Separability of True Nodules and False Positives'}, ...
    'FontSize',14,'FontWeight','bold');
legend('Location','southwest','FontSize',10,'Box','on');
xlim([0 6000]); ylim([0.0 1.05]);

exportgraphics(fig7, fullfile(OUT_DIR,'fig7_scatter_2d.png'),'Resolution',200);
fprintf('Saved fig7_scatter_2d.png\n');

%% ═══════════════════════════════════════════════════════════════════════════
%% FIGURE 8 — 3D Scatter: Eccentricity / MeanHU / Solidity
%% ═══════════════════════════════════════════════════════════════════════════
fig8 = figure('Position',[50 50 1000 700],'Color','white','Name','Scatter 3D');
ax8 = axes(fig8); hold on; box on; grid on;
set(ax8,'Color',[0.96 0.96 0.97],'GridColor',[0.55 0.55 0.55],'GridAlpha',0.4,'FontSize',11);

scatter3(gt_data(:,3), gt_data(:,4), gt_data(:,2), 130, C_GT, 'o',...
    'LineWidth',1.5,'MarkerFaceColor',C_GT,'MarkerFaceAlpha',0.75,...
    'DisplayName','True Nodule (Dice ≥ 0.50)');

scatter3(fp_data(:,3), fp_data(:,4), fp_data(:,2), 150, C_FP, 'x',...
    'LineWidth',2.5,'DisplayName','False Positive (Dice = 0)');

scatter3(partial_data(:,3), partial_data(:,4), partial_data(:,2), 130, C_PART, 'd',...
    'LineWidth',2,'MarkerFaceColor',C_PART,'MarkerFaceAlpha',0.65,...
    'DisplayName','Partial (0 < Dice < 0.50)');

% Decision plane: Eccentricity = 0.93
[emesh, hmesh] = meshgrid(linspace(0.3,1.0,10), linspace(-750,300,10));
smesh = 0.65 * ones(size(emesh));  % Solidity = sol_thr plane
surf(emesh*0+0.93, hmesh, smesh,...
    'FaceColor',[0.4 0.4 0.9],'FaceAlpha',0.12,'EdgeColor','none',...
    'DisplayName','ECC threshold = 0.93');

xlabel('Eccentricity','FontSize',12,'FontWeight','bold');
ylabel('Mean HU','FontSize',12,'FontWeight','bold');
zlabel('Solidity','FontSize',12,'FontWeight','bold');
title({'3D Feature Space — Eccentricity / MeanHU / Solidity',...
       'Algorithm decision boundaries shown as transparent planes'},...
    'FontSize',13,'FontWeight','bold');
legend('Location','northeast','FontSize',10,'Box','on');
xlim([0.3 1.0]); ylim([-750 300]); zlim([0.0 1.05]);
view([-42 28]);

% Rotate animation hint (static for PNG)
camlight('headlight'); lighting gouraud;

exportgraphics(fig8, fullfile(OUT_DIR,'fig8_scatter_3d.png'),'Resolution',200);
fprintf('Saved fig8_scatter_3d.png\n');

%% ═══════════════════════════════════════════════════════════════════════════
%% FIGURE 9 — Halo Score Distribution (GT vs FP discrimination)
%% ═══════════════════════════════════════════════════════════════════════════
fig9 = figure('Position',[50 50 900 520],'Color','white','Name','Halo Distribution');

ax9a = subplot(1,2,1); hold on; box on; grid on;
set(ax9a,'Color',C_BG,'FontSize',11);

edges_h = 0:0.05:1.05;
h_gt = histcounts(gt_data(:,5), edges_h, 'Normalization','probability');
h_fp = histcounts(fp_data(:,5), edges_h, 'Normalization','probability');

b1 = bar(edges_h(1:end-1)+0.025, h_gt, 1.0, 'FaceColor', C_GT, 'FaceAlpha',0.7,...
    'EdgeColor','none','DisplayName','True Nodule');
b2 = bar(edges_h(1:end-1)+0.025, h_fp, 0.6, 'FaceColor', C_FP, 'FaceAlpha',0.7,...
    'EdgeColor','none','DisplayName','False Positive');
xline(0.85, '--', 'ISOL threshold', 'Color',[0.2 0.2 0.8],'LineWidth',2,...
    'LabelHorizontalAlignment','left','FontSize',9);
xline(0.50, ':', 'halo > 0.50', 'Color',[0.5 0.5 0.5],'LineWidth',1.5,...
    'LabelHorizontalAlignment','left','FontSize',8);
xlabel('Halo Score','FontSize',12,'FontWeight','bold');
ylabel('Fraction','FontSize',12,'FontWeight','bold');
title('Halo Score Distribution','FontSize',13,'FontWeight','bold');
legend('Location','northwest','FontSize',10);
xlim([0 1.05]); ylim([0 0.55]);

% Solidity histogram
ax9b = subplot(1,2,2); hold on; box on; grid on;
set(ax9b,'Color',C_BG,'FontSize',11);
edges_s = 0.3:0.05:1.05;
s_gt = histcounts(gt_data(:,2), edges_s, 'Normalization','probability');
s_fp = histcounts(fp_data(:,2), edges_s, 'Normalization','probability');
bar(edges_s(1:end-1)+0.025, s_gt, 1.0, 'FaceColor',C_GT,'FaceAlpha',0.7,...
    'EdgeColor','none','DisplayName','True Nodule');
bar(edges_s(1:end-1)+0.025, s_fp, 0.6, 'FaceColor',C_FP,'FaceAlpha',0.7,...
    'EdgeColor','none','DisplayName','False Positive');
xline(0.65,'--','sol\_thr=0.65','Color',[0.2 0.2 0.8],'LineWidth',2,...
    'LabelHorizontalAlignment','left','FontSize',9);
xlabel('Solidity','FontSize',12,'FontWeight','bold');
ylabel('Fraction','FontSize',12,'FontWeight','bold');
title('Solidity Distribution','FontSize',13,'FontWeight','bold');
legend('Location','northwest','FontSize',10);
xlim([0.3 1.05]); ylim([0 0.55]);

sgtitle('Key Feature Distributions — True Nodules vs False Positives',...
    'FontSize',14,'FontWeight','bold');

exportgraphics(fig9, fullfile(OUT_DIR,'fig9_feature_distributions.png'),'Resolution',200);
fprintf('Saved fig9_feature_distributions.png\n');

close all;
fprintf('\nAll scatter figures done.\n');
