clc;
clear all;
close all;

Ra_values = [1e8, 3e8, 1e9, 3e9, 1e10, 3e10, 1e11, 3e11, 1e12];
Ra_str = {'1e8', '3e8', '1e9', '3e9', '1e10', '3e10', '1e11', '3e11', '1e12'};
num_Ra = length(Ra_values);

results = struct();
num_energy_files = 9;

%% 
fprintf('Processing energy percentage and loading data for each Ra...\n');

for i = 1:num_Ra
    results(i).summary = load(fullfile('../data', strcat(Ra_str{i}, '_inter'), 'Fourier_result_summary.dat'));
    
    energy_data_y = cell(1, num_energy_files);
    NuRe_data_y = cell(1, 4);
    error_data_y = cell(1, 2);
    
    for j = 1:num_energy_files
        energy_filename = fullfile('../data', strcat(Ra_str{i}, '_inter'), sprintf('energyPercentage_%d.dat', j));
        temp_data = load(energy_filename);
        if j == 1
            results(i).time = temp_data(:,1);
        end
        energy_data_y{j} = temp_data(:,2);
    end
        
    results(i).energy_percentages = energy_data_y;

    NuRe_filename = fullfile(strcat('../data/', Ra_str{i}, '_inter'), sprintf('NuRe_compare.dat'));
    temp_data_NuRe = load(NuRe_filename);
    for k=1:4
        NuRe_data_y{k} = temp_data_NuRe(k,:);
    end
    results(i).NuRe = NuRe_data_y;

    error_data_y{1} = abs(NuRe_data_y{2}-NuRe_data_y{1})./abs(NuRe_data_y{1})*100;
    error_data_y{2} = abs(NuRe_data_y{4}-NuRe_data_y{3})./abs(NuRe_data_y{3})*100;
    
    results(i).error = error_data_y;
end

fprintf('Done loading data.\n');

%% 
fprintf('Aggregating results and saving to a single MAT file...\n');

rows_count = size(results(1).summary, 1);
result_4 = zeros(rows_count, num_Ra);
result_5 = zeros(rows_count, num_Ra);
result_6 = zeros(rows_count, num_Ra);
result_7 = zeros(rows_count, num_Ra);

for i = 1:num_Ra
    result_4(:, i) = results(i).summary(:, 4); % <E^{m,n}>
    result_5(:, i) = results(i).summary(:, 5); % Percentage
    result_6(:, i) = results(i).summary(:, 6); % Ratio to RMS
    result_7(:, i) = results(i).summary(:, 7); % RMS value
end


% 1. <E^{m,n}>
E_mn_avg = result_4;            

% 2. E^{m,n}_{rms}
E_mn_rms = result_7;            

% 3. <E^{m,n}> / E^{m,n}_{rms}
E_mn_over_E_rms = result_6;     

% 4. <E^{m,n}> / <E_{total}> (Percentage)
E_mn_percentage = result_5;     

% 5. <E_{total}>
%  <E_{total}> = <E^{m,n}> / (Percentage / 100)
E_total_avg = result_4 ./ (result_5 ./ 100);

output_filename = 'Energy_Statistics.mat';

%行：(m,n) (1,1)-(1,2)-(1,3)-(2,1)-(2,2)-(2,3)...(3,3)
%列：Ra 1e8-1e12
save(output_filename, ...
     'E_mn_avg', ...           % <E^{m,n}>
     'E_mn_rms', ...           % E^{m,n}_{rms}
     'E_total_avg', ...        % <E_{total}>
     'E_mn_over_E_rms', ...    % <E^{m,n}> / E^{m,n}_{rms}
     'E_mn_percentage', ...    % <E^{m,n}> / <E_{total}>
     'Ra_values');             % Ra

fprintf('All tasks completed. Data saved to %s.\n', output_filename);