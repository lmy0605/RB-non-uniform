clear; close all; clc;

%% 
Ra_list = [1e9, 1e10, 1e11, 1e12];
num_files = length(Ra_list);
all_msd = cell(1, num_files);
all_lag_steps = cell(1, num_files);
all_U_potential = cell(1, num_files);
all_r_potential = cell(1, num_files);
legend_entries = cell(1, num_files);
max_lag = 100;
num_bins = 10;

%% 
for i = 1:num_files
    Ra = Ra_list(i);
    fprintf('>>>>>> 正在处理 Ra = %.0e ...\n', Ra);
    
    filename = sprintf('LSC_1e%d.mat', log10(Ra));
    fprintf('正在加载文件: %s\n', filename);
    load(filename, 'LSC_center_x', 'LSC_center_y');

    % MSD
    fprintf('  计算MSD...\n');
    num_points = length(LSC_center_x);
    lag_steps = (1:max_lag)';
    msd = zeros(max_lag, 1);

    for m = 1:max_lag
        dx = LSC_center_x(1+m : end) - LSC_center_x(1 : end-m);
        dy = LSC_center_y(1+m : end) - LSC_center_y(1 : end-m);
        squared_displacements = dx.^2 + dy.^2;
        msd(m) = mean(squared_displacements);
    end
    
    all_msd{i} = msd;
    all_lag_steps{i} = lag_steps;
    
    % 有效势能
    fprintf('  计算有效势阱...\n');
    mean_x = mean(LSC_center_x);
    mean_y = mean(LSC_center_y);
    r = sqrt((LSC_center_x - mean_x).^2 + (LSC_center_y - mean_y).^2);
    
    [counts, edges] = histcounts(r, num_bins, 'Normalization', 'pdf');
    bin_centers = (edges(1:end-1) + edges(2:end)) / 2;

    prob_density = zeros(size(counts));
    valid_indices = bin_centers > 1e-9;
    prob_density(valid_indices) = counts(valid_indices);
    
    valid_potential_indices = prob_density > 1e-9;
    r_potential = bin_centers(valid_potential_indices);
    U_potential = -log(prob_density(valid_potential_indices));
    U_potential = U_potential - min(U_potential)+1; 
    
    all_U_potential{i} = U_potential;
    all_r_potential{i} = r_potential;
    
    legend_entries{i} = sprintf('Ra = 10^{%d}', log10(Ra));
    
    fprintf('处理完成.\n\n');
end


%% 
plot_matlab('x',all_lag_steps, ...
            'y',all_msd, ...
            'CanvasRatio', [10, 7], ...
            'YLim',[1e-4,5e-1],...
            'YTicks', [1e-3, 1e-1], ...
                 'LineStyle', '-', ...
                 'XLabel', '\tau', ...
                 'YLabel', '\rm{MSD}', ...
                 'LogScale', 'xy', ...
                 'Marker', {'s'}, ...
                 'MarkerSize', 8, ...
                 'MarkerInterval', 1,...
                 'MarkerFaceColor','w',...
                 'ShowLegend', 0, ...
                 'OutputFilename', 'MSD', ...
                 'Resolution', 600);

plot_matlab('x',all_r_potential, ...
            'y',all_U_potential, ...
            'CanvasRatio', [10, 7], ...
            'YLim',[7e-1,10],...
                 'LineStyle', '-', ...
                 'XLabel', '\it{r}', ...
                 'YLabel', '\it{U}\rm{(}\it{r}\rm{)}', ...
                 'LogScale', 'y', ...
                 'Marker', {'s'}, ...
                 'MarkerSize', 8, ...
                 'MarkerInterval', 1,...
                 'MarkerFaceColor','w',...
                 'ShowLegend', 0, ... 
                 'Legend', legend_entries, ...
                 'OutputFilename', 'U', ...
                 'Resolution', 600);
hold on;

full_legend_entries = legend_entries;

legend(full_legend_entries, 'Location', 'southeast', 'Interpreter', 'tex','Box', 'off', 'FontSize', 18);
print('U.png', '-dpng', '-r600');
