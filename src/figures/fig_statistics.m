%   Author: Angelo Mastrangelo

%  Produces:
%    fig1_dice_distribution.png  — Histogram + scatter/boxplot
%    fig2_ablation_study.png     — Ablation bar chart
%    fig3_failure_analysis.png   — Failure pie + margin bars

clear; close all; clc;

%% ─── Run 105 final results ───────────────────────────────────────────────
dice = [0.8683 0.5876 0.9397 0.8519 0.0000 0.8703 0.7094 0.7986 ...
        0.7939 0.4040 0.6421 0.8686 0.9015 0.5597 0.9174 0.0000 ...
        0.0000 0.0000 0.0000 0.6505 0.9220 0.9249 0.8411 0.8973 ...
        0.8874 0.0000 0.5674 0.3957 0.9150 0.9040 0.9169 0.0000 ...
        0.9049 0.7712 0.4558 0.8677 0.7089 0.9015 0.0000 0.8824 ...
        0.9057 0.0000];

pids = {'0013','0014','0016','0027','0037','0045','0052','0061', ...
        '0064','0067','0068','0076','0094','0109','0129','0133', ...
        '0152','0158','0163','0169','0171','0172','0175','0179', ...
        '0181','0187','0195','0196','0217','0219','0237','0246', ...
        '0252','0310','0387','0469','0543','0635','0655','0673', ...
        '0705','0939'};

n          = numel(dice);
n_zeros    = sum(dice == 0);
n_above06  = sum(dice > 0.60);
mu         = mean(dice);
med        = median(dice);
sd         = std(dice);

%% ─── Color palette ───────────────────────────────────────────────────────
C_HIGH = [0.133 0.545 0.133];   % green  — Dice > 0.60
C_MID  = [1.000 0.600 0.000];   % orange — 0 < Dice ≤ 0.60
C_ZERO = [0.800 0.100 0.100];   % red    — Dice = 0
C_MEAN = [0.200 0.400 0.800];   % blue   — mean
C_MED  = [0.550 0.150 0.750];   % purple — median
C_BG   = [0.965 0.965 0.965];   % panel bg

out = 'export_figures';
if ~exist(out,'dir'), mkdir(out); end

%% ══════════════════════════════════════════════════════════════════════════
%% FIGURE 1 — Dice Score Distribution (Histogram + Scatter/Box)
%% ══════════════════════════════════════════════════════════════════════════
fig1 = figure('Position',[50 50 1300 540],'Color','white','Name','Dice Distribution');

% ── Subplot 1: Histogram ─────────────────────────────────────────────────
ax1 = subplot(1,2,1); hold on; box on; set(ax1,'Color',C_BG);
edges   = -0.025 : 0.05 : 1.025;
centers = edges(1:end-1) + 0.025;

% Separate into 3 groups for coloring
for i = 1:numel(centers)
    lo = edges(i); hi = edges(i+1);
    cnt_z = sum(dice == 0   & dice >= lo & dice < hi + 1e-9*(hi==1));
    cnt_m = sum(dice >  0   & dice <= 0.6 & dice >= lo & dice < hi);
    cnt_h = sum(dice >  0.6 & dice >= lo & dice < hi);
    if cnt_z > 0, bar(centers(i), cnt_z, edges(2)-edges(1)-0.002, ...
            'FaceColor',C_ZERO,'EdgeColor','none'); end
    if cnt_m > 0, bar(centers(i), cnt_m, edges(2)-edges(1)-0.002, ...
            'FaceColor',C_MID,'EdgeColor','none'); end
    if cnt_h > 0, bar(centers(i), cnt_h, edges(2)-edges(1)-0.002, ...
            'FaceColor',C_HIGH,'EdgeColor','none'); end
end

ylim_top = 9;
plot([mu  mu ],[0 ylim_top],'--','Color',C_MEAN,'LineWidth',2.5);
plot([med med],[0 ylim_top],'--','Color',C_MED, 'LineWidth',2.5);
plot([0.60 0.60],[0 ylim_top],':','Color',[0.5 0.5 0.5],'LineWidth',1.8);

text(0.01, 7.6, sprintf('%d zeri\n(limite 2D)', n_zeros), ...
    'Color',C_ZERO,'FontSize',10,'FontWeight','bold');
text(0.75, 7.6, sprintf('%d/42 slices\n Dice>0.60', n_above06), ...
    'Color',C_HIGH,'FontSize',10,'FontWeight','bold');
text(mu+0.015, 6.5, sprintf('\\mu=%.4f', mu), ...
    'Color',C_MEAN,'FontSize',10,'FontWeight','bold');
text(med+0.015, 5.5, sprintf('\\theta_{50}=%.4f', med), ...
    'Color',C_MED,'FontSize',10,'FontWeight','bold');

xlabel('Dice Score','FontSize',13,'FontWeight','bold');
ylabel('N° Pazienti','FontSize',13,'FontWeight','bold');
title('Distribuzione del Dice Score','FontSize',14,'FontWeight','bold');
xlim([-0.05 1.05]); ylim([0 ylim_top]);
legend({'Dice = 0','0 < Dice \leq 0.60','Dice > 0.60','Mean','Median','Soglia 0.60'}, ...
    'Location','northwest','FontSize',9);
grid on; grid minor;

% ── Subplot 2: Strip chart + Boxplot ─────────────────────────────────────
ax2 = subplot(1,2,2); hold on; box on; set(ax2,'Color',C_BG);
rng(42);
jit = 0.18*(rand(1,n)-0.5);
for i = 1:n
    if     dice(i) == 0,   mc = C_ZERO;
    elseif dice(i) <= 0.6, mc = C_MID;
    else,                  mc = C_HIGH; end
    scatter(1+jit(i), dice(i), 70, mc, 'filled', ...
        'MarkerEdgeColor','none','MarkerFaceAlpha',0.85);
end

% Boxplot overlay (manual IQR box)
q1 = quantile(dice,0.25); q3 = quantile(dice,0.75);
rectangle('Position',[0.8 q1 0.4 (q3-q1)],'FaceColor',[0.9 0.9 0.9 0.5], ...
    'EdgeColor',[0.3 0.3 0.3],'LineWidth',2);
plot([0.8 1.2],[med med],'-','Color',[0.3 0.3 0.3],'LineWidth',3);
plot([1 1],[min(dice(dice>0)) q1],'k-','LineWidth',1.5);
plot([1 1],[q3 max(dice)],'k-','LineWidth',1.5);
plot([0.88 1.12],[min(dice(dice>0)) min(dice(dice>0))],'k-','LineWidth',1.5);
plot([0.88 1.12],[max(dice) max(dice)],'k-','LineWidth',1.5);
plot([0.8 1.2],[mu mu],'--','Color',C_MEAN,'LineWidth',2.5);
plot([0.65 1.35],[0.60 0.60],':','Color',[0.5 0.5 0.5],'LineWidth',1.8);

stats_str = sprintf('n = %d\n\\mu = %.4f\n\\theta_{50} = %.4f\n\\sigma = %.4f\nZeri = %d/%d\nDice>0.60 = %d/%d', ...
    n, mu, med, sd, n_zeros, n, n_above06, n);
text(1.38, 0.35, stats_str, 'FontSize',10,'VerticalAlignment','middle', ...
    'BackgroundColor',[0.95 0.95 1],'EdgeColor',[0.7 0.7 0.8],'Margin',4);

ylabel('Dice Score','FontSize',13,'FontWeight','bold');
title('Distribuzione per Paziente','FontSize',14,'FontWeight','bold');
set(ax2,'XTick',[],'YLim',[-0.05 1.05],'XLim',[0.5 1.9]);
grid on; grid minor;

sgtitle('Analisi Prestazioni — Run 105 (Segmentazione Noduli Polmonari)', ...
    'FontSize',15,'FontWeight','bold');
print(fig1, fullfile(out,'fig1_dice_distribution'), '-dpng','-r300');
fprintf('  Fig 1 salvata.\n');

%% ══════════════════════════════════════════════════════════════════════════
%% FIGURE 2 — Ablation Study
%% ══════════════════════════════════════════════════════════════════════════
fig2 = figure('Position',[50 50 1300 580],'Color','white','Name','Ablation Study');
ax_ab = axes('Parent',fig2); hold on; box on; set(ax_ab,'Color',C_BG);

run_labels = { ...
    sprintf('Baseline\n(Run 41)'), ...
    sprintf('+Halo Score\n(Run 50)'), ...
    sprintf('+2.5D Tube Filter\n+Shape (Run 88)'), ...
    sprintf('+ECC Gate\nCondiz. (Run 96)'), ...
    sprintf('+Branch GGO\nEnfisema (Run 101)'), ...
    sprintf('+Guard Solidity\nEnfisema (Run 105)')};
mean_vals = [0.3060, 0.3161, 0.5742, 0.5972, 0.6080, 0.6175];
zero_vals = [33,     33,     14,     12,     10,     9    ];
nr = numel(mean_vals);

% Color gradient red → green
reds   = linspace(0.85, 0.10, nr)';
greens = linspace(0.10, 0.55, nr)';
blues  = linspace(0.10, 0.15, nr)';
bar_colors = [reds greens blues];

% Background bands
yline(0.60,'--','Color',[0.8 0.1 0.1],'LineWidth',2.5,'Label','  Target 0.60', ...
    'LabelHorizontalAlignment','right','FontSize',11,'FontWeight','bold');

bh = bar(1:nr, mean_vals, 0.62, 'FaceColor','flat','EdgeColor','none');
for i = 1:nr, bh.CData(i,:) = bar_colors(i,:); end

% Labels on bars
for i = 1:nr
    text(i, mean_vals(i)+0.012, sprintf('%.4f', mean_vals(i)), ...
        'HorizontalAlignment','center','FontSize',12,'FontWeight','bold');
    if zero_vals(i) > 0
        text(i, 0.025, sprintf('%d\\bigcirc', zero_vals(i)), ...
            'HorizontalAlignment','center','FontSize',10,'Color','white','FontWeight','bold');
    end
end

% Delta annotations between bars
for i = 2:nr
    delta = mean_vals(i) - mean_vals(i-1);
    if delta > 0.002
        mid_y  = mean_vals(i-1) + delta/2;
        text(i-0.5, mid_y + 0.04, sprintf('+%.4f', delta), ...
            'Color',[0 0.45 0],'FontSize',9.5,'FontWeight','bold', ...
            'HorizontalAlignment','center');
        annotation('textarrow', [((i-1)/nr)*0.82+0.08 ((i)/nr)*0.82+0.08], ...
            [mid_y/0.72*0.75+0.10 mid_y/0.72*0.75+0.10], ...
            'Color',[0 0.45 0],'HeadSize',6,'LineWidth',1.2);
    end
end

% The breakthrough arrow for halo
annotation('textbox',[0.34 0.72 0.12 0.08],'String', ...
    sprintf('BREAKTHROUGH\nHalo Score\n+0.2581'), ...
    'FontSize',9,'FontWeight','bold','Color',[0.1 0.5 0.1], ...
    'BackgroundColor',[0.9 1 0.9],'EdgeColor',[0.1 0.5 0.1],'FitBoxToText','on');

set(ax_ab,'XTick',1:nr,'XTickLabel',run_labels,'YLim',[0 0.72],'FontSize',10);
xlabel('Modulo Algoritmo','FontSize',13,'FontWeight','bold');
ylabel('Mean Dice Score','FontSize',13,'FontWeight','bold');
title('Ablation Study — Evoluzione del Dice Score per Modulo','FontSize',14,'FontWeight','bold');
grid on;

print(fig2, fullfile(out,'fig2_ablation_study'), '-dpng','-r300');
fprintf('  Fig 2 salvata.\n');

%% ══════════════════════════════════════════════════════════════════════════
%% FIGURE 3 — Failure Analysis
%% ══════════════════════════════════════════════════════════════════════════
fig3 = figure('Position',[50 50 1300 560],'Color','white','Name','Failure Analysis');

% ── Left: Pie ────────────────────────────────────────────────────────────
ax3 = subplot(1,2,1);
pie_v = [4, 3, 1, 1];
pie_lb = {sprintf('Vaso sanguigno\ntrasversale (n=4)'), ...
          sprintf('Struttura pleurica\n/ parete (n=3)'), ...
          sprintf('GT morfologica-\nmente irregolare (n=1)'), ...
          sprintf('GGO sotto\nsoglia solidity (n=1)')};
pie_c = [0.85 0.18 0.18; 1.00 0.55 0.08; 0.55 0.10 0.60; 0.15 0.45 0.82];
explode_v = [1 0 0 1];

axes(ax3);
ph = pie(pie_v, explode_v);
for i = 1:length(pie_v)
    ph(2*i-1).FaceColor  = pie_c(i,:);
    ph(2*i-1).EdgeColor  = 'white';
    ph(2*i-1).LineWidth  = 2;
    ph(2*i).FontSize     = 9.5;
    ph(2*i).FontWeight   = 'bold';
    ph(2*i).String       = pie_lb{i};
end
title(sprintf('Cause Fisiche dei %d Zeri Rimanenti', n_zeros), ...
    'FontSize',13,'FontWeight','bold');

% ── Right: Horizontal bar — margin FP/GT ─────────────────────────────────
ax4 = subplot(1,2,2); hold on; box on; set(ax4,'Color',C_BG);

z_ids     = {'0037','0133','0152','0158','0163','0187','0246','0655','0939'};
z_margin  = [NaN,    2.8,   2.6,  1.7,   2.6,  1.8,   1.6,  NaN,   1.3 ];
z_cause   = [4,      1,     2,    1,     2,    2,     1,    3,     1   ];
z_notes   = {'Blob area 10759px, sol=0.129', ...
             'FP ISOL halo=1.00 2.8×GT', ...
             'FP solid A=806px tube≈0.92', ...
             'GT A=51px ov=38%', ...
             'FP juxta ECC=0.975 2.6×GT', ...
             'FP sol=0.827 juxta 1.8×GT', ...
             'FP ISOL halo=1.00 1.6×GT', ...
             'GT sol=0.406 sotto T1 thr', ...
             'FP quasi identico 1.3×GT'};

c_map = [pie_c(1,:); pie_c(2,:); pie_c(4,:); pie_c(3,:)];

for i = 1:numel(z_ids)
    val = z_margin(i);
    if isnan(val), val = 0.1; end
    bh2 = barh(i, val, 0.62, 'FaceColor', c_map(z_cause(i),:), ...
        'EdgeColor','none','FaceAlpha',0.85);
    if ~isnan(z_margin(i))
        text(z_margin(i)+0.04, i, sprintf('%.1f×', z_margin(i)), ...
            'FontSize',10,'FontWeight','bold','VerticalAlignment','middle');
    else
        text(0.12, i, 'Strutturale', 'FontSize',10,'FontWeight','bold', ...
            'Color','white','VerticalAlignment','middle');
    end
    text(-0.05, i, z_notes{i}, 'FontSize',7.5,'HorizontalAlignment','right', ...
        'Color',[0.35 0.35 0.35],'VerticalAlignment','middle');
end

plot([1 1],[0.4 numel(z_ids)+0.4],':k','LineWidth',1.5);
text(1.02, numel(z_ids)+0.55, 'Pareggio','FontSize',9,'Color',[0.4 0.4 0.4]);

set(ax4,'YTick',1:numel(z_ids), ...
    'YTickLabel', cellfun(@(x) ['LIDC-' x], z_ids,'UniformOutput',false), ...
    'FontSize',10,'XLim',[-1.8 3.4]);
xlabel('Score FP / Score GT  (volte)','FontSize',12,'FontWeight','bold');
title('Margine di Sconfitta per Zero','FontSize',13,'FontWeight','bold');
grid on; grid minor;

legend({'Vaso trasversale','Struttura pleurica','GGO borderline','GT irregolare'}, ...
    'Location','southeast','FontSize',9);

sgtitle('Analisi dei Fallimenti — 9 Zeri Rimanenti (Limite Fisico 2D)', ...
    'FontSize',14,'FontWeight','bold');
print(fig3, fullfile(out,'fig3_failure_analysis'), '-dpng','-r300');
fprintf('  Fig 3 salvata.\n');

fprintf('\n=== Tutte le figure statistiche salvate in: %s/ ===\n', out);
