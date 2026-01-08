clc;
clear all;
close all;

%% 
load('Energy_Statistics.mat'); 
% 变量: 
% <E^{m,n}>
 % E^{m,n}_{rms}
% <E_{total}>
% <E^{m,n}> / E^{m,n}_{rms}
 % <E^{m,n}> / <E_{total}>
% Ra

%行：(m,n) (1,1)-(1,2)-(1,3)-(2,1)-(2,2)-(2,3)...(3,3)
%列：Ra 1e8-1e12
%% 
plot_indices = [1, 2, 4, 5]; 
num_plots = length(plot_indices);

marker_styles = {'o', '^', 'd', 's'};

custom_colors = {
    [0 0.4470 0.7410], ... % Blue
    [0.8500 0.3250 0.0980], ... % Red/Orange
    [0.4660 0.6740 0.1880], ... % Green
    [0.4940 0.1840 0.5560]      % Purple
};

% 通用軸設定
x_limits = [5e7, 5e12];
x_ticks = [1e8, 1e9, 1e10, 1e11, 1e12];
line_width = 2;       % 線條/邊框寬度
marker_size = 10;     % 標記大小

fprintf('Plotting figures...\n');

%% 3. 繪製圖 1: 能量百分比 (<E^{m,n}> / <E_{total}>)

figure(1);
% 設定畫布大小 (寬, 高)
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 800, 600], 'Color', 'w');
hold on;

for k = 1:num_plots
    idx = plot_indices(k);
    
    % 使用 semilogx 繪製對數 X 軸
    % 'LineStyle', 'none' 表示不連線，只畫點
    h = semilogx(Ra_values, E_mn_percentage(idx, :), ...
        'LineStyle', 'none', ...
        'Marker', marker_styles{k}, ...
        'Color', custom_colors{k}, ...
        'MarkerFaceColor', 'none', ... % 空心標記 (若要實心改為 custom_colors{k})
        'MarkerEdgeColor', custom_colors{k}, ...
        'MarkerSize', marker_size, ...
        'LineWidth', line_width);
end

% 座標軸設定
xlim(x_limits);
ylim([0, 60]);
set(gca, 'XTick', x_ticks);
set(gca, 'XScale', 'log'); % 確保 X 軸為對數
box on; % 顯示方框
grid off; % 顯示網格 (可選)

% 字型與標籤設定
xlabel('\it{Ra}', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('\it{<E^{m,n}>} \it{/<E_{total}>} \rm{(%)}', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'LineWidth', 1.5); % 座標軸刻度字體大小與線寬

hold off;
% 儲存圖片
saveas(gcf, 'energyPercentageAvg_divide4.png');
fprintf('Saved energyPercentageAvg_divide4.png\n');


%% 4. 繪製圖 2: 與 RMS 的比值 (<E^{m,n}>/E^{m,n}_{rms})

figure(2);
set(gcf, 'Units', 'pixels', 'Position', [150, 150, 800, 600], 'Color', 'w');
hold on;

for k = 1:num_plots
    idx = plot_indices(k);
    
    h = semilogx(Ra_values, E_mn_over_E_rms(idx, :), ...
        'LineStyle', 'none', ...
        'Marker', marker_styles{k}, ...
        'Color', custom_colors{k}, ...
        'MarkerFaceColor', 'none', ...
        'MarkerEdgeColor', custom_colors{k}, ...
        'MarkerSize', marker_size, ...
        'LineWidth', line_width);
end

% 座標軸設定
xlim(x_limits);
% Y軸範圍根據數據自動調整，或者你可以手動設定，例如 ylim([0, 5]);
set(gca, 'XTick', x_ticks);
set(gca, 'XScale', 'log');
box on;
grid off;

% 字型與標籤設定
xlabel('\it{Ra}', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('\it{<E^{m,n}>/E^{m,n}_{rms}}', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'LineWidth', 1.5);

hold off;
% 儲存圖片
saveas(gcf, 'S4.png');
fprintf('Saved S4.png\n');

fprintf('All plotting tasks completed.\n');