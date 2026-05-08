%   Author: Angelo Mastrangelo

% Genera due grafici a barre:
%   1) metriche globali (includono anche i casi con Dice = 0)
%   2) metriche sui soli casi rilevati (Dice > 0)
%
% Richiede in workspace:
%   all_dice, all_prec, all_sens

clearvars -except all_dice all_prec all_sens all_labels
clc; close all;

OUT_DIR = 'export_figures_tesi';
if ~exist(OUT_DIR,'dir'), mkdir(OUT_DIR); end

% Impostazioni globali
set(0, 'DefaultAxesColor', 'w', ...
       'DefaultTextColor', 'k', ...
       'DefaultAxesXColor', 'k', ...
       'DefaultAxesYColor', 'k');

if ~exist('all_dice','var') || ~exist('all_prec','var') || ~exist('all_sens','var')
    error('Variabili all_dice / all_prec / all_sens non trovate nel workspace.');
end

% Nomi metriche
metrics_names = {'Dice Score', 'Precision', 'Sensitivity (Recall)'};

% Colori: Blu, Arancio, Verde
colors = [0.18 0.53 0.96;
          0.95 0.60 0.20;
          0.15 0.70 0.35];

%% ────────────────────────────────────────────────────────────────────────
% 1) METRICHE GLOBALI
%% ────────────────────────────────────────────────────────────────────────
mean_dice_global = mean(all_dice);
mean_prec_global = mean(all_prec);
mean_sens_global = mean(all_sens);

metrics_vals_global = [mean_dice_global, mean_prec_global, mean_sens_global];

fig1 = figure('Position',[100 100 900 650], 'Color','w', 'Name','Metriche Globali');
ax1 = axes('Parent',fig1); hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on');
set(ax1, 'GridColor', [0.85 0.85 0.85], ...
         'GridAlpha', 0.6, ...
         'FontSize', 13);

for i = 1:3
    bar(ax1, i, metrics_vals_global(i), 0.55, ...
        'FaceColor', colors(i,:), ...
        'EdgeColor', 'k', ...
        'LineWidth', 1.5);
    text(i, metrics_vals_global(i) + 0.03, sprintf('%.4f', metrics_vals_global(i)), ...
        'HorizontalAlignment', 'center', ...
        'FontSize', 15, ...
        'FontWeight', 'bold', ...
        'Color', 'k');
end

set(ax1, 'XTick', 1:3, 'XTickLabel', metrics_names, 'XTickLabelRotation', 0);
ylabel('Punteggio Medio (0.0 - 1.0)', 'FontSize', 15, 'FontWeight', 'bold');
title('Prestazioni Medie Globali', 'FontSize', 18, 'FontWeight', 'bold');
ylim([0 1.1]);
xlim([0.3 3.7]);

exportgraphics(fig1, fullfile(OUT_DIR, '04_metrics_barchart_global.png'), 'Resolution', 300);

%% ────────────────────────────────────────────────────────────────────────
% 2) METRICHE SUI CASI RILEVATI (Dice > 0)
%% ────────────────────────────────────────────────────────────────────────
valid_idx = (all_dice > 0);
if sum(valid_idx) == 0
    warning('Nessun caso con Dice > 0 trovato. Uso tutte le slice come fallback.');
    valid_idx = true(size(all_dice));
end

mean_dice_valid = mean(all_dice(valid_idx));
mean_prec_valid = mean(all_prec(valid_idx));
mean_sens_valid = mean(all_sens(valid_idx));

metrics_vals_valid = [mean_dice_valid, mean_prec_valid, mean_sens_valid];

fig2 = figure('Position',[120 120 900 650], 'Color','w', 'Name','Metriche Casi Rilevati');
ax2 = axes('Parent',fig2); hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on');
set(ax2, 'GridColor', [0.85 0.85 0.85], ...
         'GridAlpha', 0.6, ...
         'FontSize', 13);

for i = 1:3
    bar(ax2, i, metrics_vals_valid(i), 0.55, ...
        'FaceColor', colors(i,:), ...
        'EdgeColor', 'k', ...
        'LineWidth', 1.5);
    text(i, metrics_vals_valid(i) + 0.03, sprintf('%.4f', metrics_vals_valid(i)), ...
        'HorizontalAlignment', 'center', ...
        'FontSize', 15, ...
        'FontWeight', 'bold', ...
        'Color', 'k');
end

set(ax2, 'XTick', 1:3, 'XTickLabel', metrics_names, 'XTickLabelRotation', 0);
ylabel('Punteggio Medio (0.0 - 1.0)', 'FontSize', 15, 'FontWeight', 'bold');
title('Prestazioni Medie sui Pazienti Rilevati', 'FontSize', 18, 'FontWeight', 'bold');
ylim([0 1.1]);
xlim([0.3 3.7]);

exportgraphics(fig2, fullfile(OUT_DIR, '05_metrics_barchart_detected_only.png'), 'Resolution', 300);

%% ────────────────────────────────────────────────────────────────────────
% Stampa riepilogo console
%% ────────────────────────────────────────────────────────────────────────
fprintf('\n');
fprintf('═══════════════════════════════════════════════\n');
fprintf('  METRICHE GLOBALI\n');
fprintf('═══════════════════════════════════════════════\n');
fprintf('  Mean Dice        : %.4f\n', mean_dice_global);
fprintf('  Mean Precision   : %.4f\n', mean_prec_global);
fprintf('  Mean Sensitivity : %.4f\n', mean_sens_global);

fprintf('\n');
fprintf('═══════════════════════════════════════════════\n');
fprintf('  METRICHE SUI CASI RILEVATI (Dice > 0)\n');
fprintf('═══════════════════════════════════════════════\n');
fprintf('  Casi rilevati    : %d / %d\n', sum(valid_idx), length(all_dice));
fprintf('  Mean Dice        : %.4f\n', mean_dice_valid);
fprintf('  Mean Precision   : %.4f\n', mean_prec_valid);
fprintf('  Mean Sensitivity : %.4f\n', mean_sens_valid);

fprintf('\nGrafici salvati in: %s/\n', OUT_DIR);