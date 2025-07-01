clear; close all; clc;

%% basic settings
fileNumStart=1501;
fileNumEnd=10000;
casename='1e9';
fileNumInterval=1;
fileSum=fileNumEnd-fileNumStart+1;
inputDir_Nu = strcat('C:\Users\user\OneDrive\Desktop\01-RB-0612\fig3-timeseries-0\',casename,'\timeseries_ReNu.plt');
inputDir_Nurest =strcat('C:\Users\user\OneDrive\Desktop\01-RB-0612\fig4_timeseries-stationary\',casename,'\timeseries_Nu.plt');
namebase = 'buoyancyCavity-';

nx=513;
ny=nx;
constA=2.1;
Rayleigh=1e9;
Prandtl=0.71;

params = calculateSystemParameters(nx,ny, Rayleigh, Prandtl,constA,'log.log');
viscosity=sqrt(Prandtl/Rayleigh);

%% Load Data
disp('Loading data from .plt files...');
[existing_data_cell, var_names] = read_tecplot_plt(inputDir_Nu);
NuWallAvg  = existing_data_cell{strcmp(var_names, 'Nu')}';
NuWallAvg  = NuWallAvg(fileNumStart:fileNumEnd);

[existing_data_cell, var_names] = read_tecplot_plt(inputDir_Nurest);
NuVolAvg  = existing_data_cell{strcmp(var_names, 'NuVolAvg')}';
NuEUAvg   = existing_data_cell{strcmp(var_names, 'NuEUAvg')}';
NuETAvg   = existing_data_cell{strcmp(var_names, 'NuETAvg')}';
disp('Data loading complete.');

%% --- log original statistics ---
disp('Calculating and logging original statistics...');

mean_NuWallAvg = mean(NuWallAvg);
std_NuWallAvg = std(NuWallAvg, 1);
mean_NuVolAvg = mean(NuVolAvg);
std_NuVolAvg = std(NuVolAvg, 1);
mean_NuEUAvg = mean(NuEUAvg);
std_NuEUAvg = std(NuEUAvg, 1);
mean_NuETAvg = mean(NuETAvg);
std_NuETAvg = std(NuETAvg, 1);
log_filename = strcat('statistics_log_',casename,'.txt');
fid = fopen(log_filename, 'w');
if fid == -1, error('Could not open log file for writing.'); end
fprintf(fid, 'Statistics Log for Nusselt Number Time Series\n');
fprintf(fid, '==================================================\n');
fprintf(fid, 'Date: %s\n', datestr(now));
fprintf(fid, 'Rayleigh Number: %e\n', Rayleigh);
fprintf(fid, 'Data points used: %d\n\n', length(NuWallAvg));
fprintf(fid, '--- Original Data Statistics (Population Std Dev) ---\n');
fprintf(fid, 'Variable      |          Mean |   Std Deviation\n');
fprintf(fid, '--------------------------------------------------\n');
fprintf(fid, 'NuWallAvg     | %13.6f | %15.6f\n', mean_NuWallAvg, std_NuWallAvg);
fprintf(fid, 'NuVolAvg      | %13.6f | %15.6f\n', mean_NuVolAvg, std_NuVolAvg);
fprintf(fid, 'NuEUAvg       | %13.6f | %15.6f\n', mean_NuEUAvg, std_NuEUAvg);
fprintf(fid, 'NuETAvg       | %13.6f | %15.6f\n', mean_NuETAvg, std_NuETAvg);
fprintf(fid, '--------------------------------------------------\n');
fclose(fid);
disp(['Statistics have been written to ', log_filename]);

%% Standardize the data
disp('Standardizing data for PDF calculation...');
NuWallAvg_std = (NuWallAvg - mean_NuWallAvg) / std_NuWallAvg;
NuVolAvg_std  = (NuVolAvg - mean_NuVolAvg) / std_NuVolAvg;
NuEUAvg_std   = (NuEUAvg - mean_NuEUAvg) / std_NuEUAvg;
NuETAvg_std   = (NuETAvg - mean_NuETAvg) / std_NuETAvg;

%% Calculate PDFs
disp('Calculating probability density functions (PDFs)...');
all_std_data_for_bins = [NuWallAvg_std(:); NuVolAvg_std(:); NuEUAvg_std(:); NuETAvg_std(:)];
min_val = min(all_std_data_for_bins);
max_val = max(all_std_data_for_bins);
num_bins = 50;
common_edges = linspace(min_val, max_val, num_bins + 1);
[pdf_values_wallavg, ~] = histcounts(NuWallAvg_std, common_edges, 'Normalization', 'pdf');
[pdf_values_NuVolAvg, ~] = histcounts(NuVolAvg_std, common_edges, 'Normalization', 'pdf');
[pdf_values_NuEUAvg, ~] = histcounts(NuEUAvg_std, common_edges, 'Normalization', 'pdf');
[pdf_values_NuETAvg, ~] = histcounts(NuETAvg_std, common_edges, 'Normalization', 'pdf');
bin_centers = common_edges(1:end-1) + diff(common_edges)/2;

%% Write PDF data to files
disp('Writing PDF data to .plt files...');
tec_file = liton_ordered_tec.TEC_FILE;
tec_file.FileName = strcat('PDF_Nu_',casename,'_hist');
tec_file.Variables = {'Nu_std_center','PDF_NuWallAvg','PDF_NuVolAvg','PDF_NuEUAvg','PDF_NuETAvg'};
tec_file.Zones = liton_ordered_tec.TEC_ZONE;
tec_file.Zones.Data = {bin_centers', pdf_values_wallavg', pdf_values_NuVolAvg', pdf_values_NuEUAvg', pdf_values_NuETAvg'};
tec_file = tec_file.write_plt();
disp('Empirical PDF data saved.');

%% --- Distribution Fitting Analysis ---
disp(newline);
disp('--- Starting Distribution Fitting Analysis for All Nu Types ---');

% --- 1. Open log file in APPEND mode ---
fid = fopen(log_filename, 'a');
if fid == -1, error('Could not open log file for appending.'); end

fprintf(fid, '--- Distribution Fitting Analysis ---\n');

all_data_to_fit = {NuWallAvg_std, NuVolAvg_std, NuEUAvg_std, NuETAvg_std};
all_pdf_values  = {pdf_values_wallavg, pdf_values_NuVolAvg, pdf_values_NuEUAvg, pdf_values_NuETAvg};
data_names      = {'NuWallAvg', 'NuVolAvg', 'NuEUAvg', 'NuETAvg'};

for i = 1:length(all_data_to_fit)

    fprintf('\n======================================================\n');
    fprintf('           ANALYZING: %s\n', data_names{i});
    fprintf('======================================================\n');

    fprintf(fid, '\n======================================================\n');
    fprintf(fid, '           ANALYZING: %s\n', data_names{i});
    fprintf(fid, '======================================================\n');
    
    % --- 2. Pass file ID (fid) to the analysis function ---
    analyze_and_plot_distribution( ...
        all_data_to_fit{i}, ...
        data_names{i}, ...
        bin_centers, ...
        all_pdf_values{i}, ...
        casename, ...
        fid ... 
    );
end

% --- 3. Close the log file ---
fclose(fid);
disp(['Distribution fitting results have been appended to ', log_filename]);
disp('All distribution fitting analyses are complete.');

%% --- Functions ---

function analyze_and_plot_distribution(data_std, data_name, bin_centers, pdf_empirical, casename, fid)
    % This function performs distribution fitting, plotting, and comparison
    % and writes the text output to the file identifier 'fid'.

    data_to_fit = data_std'; 

    % --- 1. Fit Distributions ---
    disp('Fitting Gaussian and GEV distributions...');
    pd_gaussian = fitdist(data_to_fit, 'Normal');
    pd_gev = fitdist(data_to_fit, 'Generalized Extreme Value');

    fprintf(fid, '\n--- Fit Results for %s ---\n', data_name);
    fprintf(fid, 'Gaussian Fit Parameters:\n');
    fprintf(fid, '    mu: %g,  sigma: %g\n', pd_gaussian.mu, pd_gaussian.sigma);
    fprintf(fid, 'GEV Fit Parameters:\n');
    fprintf(fid, '    k (shape): %g,  sigma (scale): %g,  mu (location): %g\n', pd_gev.k, pd_gev.sigma, pd_gev.mu);
    fprintf(fid, '-------------------------------------\n');

    % --- 2. Generate PDF curves from fitted distributions ---
    x_values = bin_centers'; 
    y_gaussian_fit = pdf(pd_gaussian, x_values);
    y_gev_fit = pdf(pd_gev, x_values);

    % --- 3. Visual Comparison: Plotting (no change needed here) ---
    figure('Visible', 'off');
    hold on;
    bar(bin_centers, pdf_empirical, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none', 'BarWidth', 1);
    plot(x_values, y_gaussian_fit, 'b-', 'LineWidth', 2);
    plot(x_values, y_gev_fit, 'r--', 'LineWidth', 2);
    title(sprintf('PDF of Standardized %s and Fitted Distributions', strrep(data_name, '_', '\_')));
    xlabel(sprintf('(%s - <%s>) / \\sigma', data_name, data_name));
    ylabel('Probability Density Function (PDF)');
    legend('Empirical PDF', 'Gaussian Fit', 'GEV Fit', 'Location', 'best');
    grid on;
    hold off;
    
    plot_filename = sprintf('PDF_Fit_Comparison_%s_%s.png', data_name, casename);
    saveas(gcf, plot_filename);
    disp(['Fit comparison plot saved as ', plot_filename]);
    fprintf(fid, 'Fit comparison plot saved as %s\n', plot_filename);
    close(gcf);

    % --- 4. Quantitative Comparison (Goodness-of-Fit) ---
    logL_gaussian = -pd_gaussian.NLogL;
    [aic_gaussian, bic_gaussian] = aicbic(logL_gaussian, pd_gaussian.NumParameters, length(data_to_fit));
    logL_gev = -pd_gev.NLogL;
    [aic_gev, bic_gev] = aicbic(logL_gev, pd_gev.NumParameters, length(data_to_fit));

    fprintf(fid, '\n--- Goodness-of-Fit Comparison for %s ---\n', data_name);
    fprintf(fid, 'Distribution | Log-Likelihood |   AIC    |   BIC    \n');
    fprintf(fid, '-----------------------------------------------------\n');
    fprintf(fid, 'Gaussian     | %14.4f | %8.2f | %8.2f \n', logL_gaussian, aic_gaussian, bic_gaussian);
    fprintf(fid, 'GEV          | %14.4f | %8.2f | %8.2f \n', logL_gev, aic_gev, bic_gev);
    fprintf(fid, '-----------------------------------------------------\n');

    if aic_gev < aic_gaussian
        fprintf(fid, 'Conclusion: GEV distribution provides a better fit (lower AIC).\n');
    else
        fprintf(fid, 'Conclusion: Gaussian distribution provides a better fit or is comparable (lower or equal AIC).\n');
    end

    % --- 5. Write fitted curves to a .plt file ---
    fit_filename = sprintf('PDF_Fits_%s_%s', data_name, casename);
    disp(['Writing fitted PDF curves to ', fit_filename, '.plt ...']);
    fprintf(fid, 'Fitted PDF data saved to %s.plt\n', fit_filename);
    
    fit_tec_file = liton_ordered_tec.TEC_FILE;
    fit_tec_file.FileName = fit_filename;
    fit_tec_file.Variables = {'Nu_std', 'PDF_Empirical', 'PDF_Gaussian_Fit', 'PDF_GEV_Fit'};
    fit_tec_file.Zones = liton_ordered_tec.TEC_ZONE;
    fit_tec_file.Zones.Data = {x_values, pdf_empirical', y_gaussian_fit, y_gev_fit};
    fit_tec_file.write_plt();
end


function [U, V ,T, rho] = readBinaryFile(file, nx, ny)
    % (The rest of the function code is identical to your original script)
    fid = fopen(file,'r');
    [~,~] = fread(fid,1,'int32');
    U = fread(fid,nx*ny,'float64');
    [~,~] = fread(fid,2,'int32');
    V = fread(fid,nx*ny,'float64');
    [~,~] = fread(fid,2,'int32');
    T = fread(fid,nx*ny,'float64');
    [~,~] = fread(fid,2,'int32');
    rho = fread(fid,nx*ny,'float64');
    [~,~] = fread(fid,1,'int32');
    fclose(fid);
end

function [U_dxWeighted,U_dyWeighted,U_areaWeighted]=nonUniformAverage(U,xGrid,yGrid)
    [nx, ny] = size(U);
    dx = zeros(1,nx); dy = zeros(1,ny);
    dx(1) = (xGrid(2) -0)/2; 
    for i = 2:nx, dx(i) = (xGrid(i+1) - xGrid(i-1))/2; end
    dy(1) = (yGrid(2) - 0)/2; 
    for j = 2:ny, dy(j) = (yGrid(j+1) - yGrid(j-1))/2; end
    U_dxWeighted = zeros(nx,ny); U_dyWeighted = zeros(nx,ny); U_areaWeighted = zeros(nx,ny);
    for j = 1:ny
        for i = 1:nx
            U_dxWeighted(i,j) = U(i,j) .* dx(i);      
            U_dyWeighted(i,j) = U(i,j) .* dy(j);      
            U_areaWeighted(i,j) = U(i,j) .* dx(i) .* dy(j);
        end
    end
end

function [deri_Ty_bottom,deri_Ty_top]=deri_Twall(T,T_bottom,T_top,length_LB,q,nx,ny)
    deri_Ty_bottom=zeros(nx,1); deri_Ty_top=zeros(nx,1);
    for i=1:1:nx
        deri_Ty_bottom(i)=(-4*q*(q+1)*T_bottom+(2*q+1)^2*T(i,1)-T(i,2))/(q*(2*q+1))*length_LB;
        deri_Ty_top(i)=(4*q*(q+1)*T_top-(2*q+1)^2*T(i,ny)+T(i,ny-1))/(q*(2*q+1))*length_LB;
    end
end

function meanU = cumulativeAverage(U)
    if ~isvector(U), error('input is not an array'); end
    U = U(:).'; 
    cumulative_sum = cumsum(U);
    divisors = 1:length(U);
    meanU = cumulative_sum ./ divisors;
end

function varPopU = cumulativePopulationVariance(U)
    if ~isvector(U), error('input is not an array'); end
    n = length(U);
    if n == 0, varPopU = []; return; end
    U_row = U(:).';
    cumulative_sum = cumsum(U_row);
    cumulative_sum_sq = cumsum(U_row.^2);
    t = 1:n;
    meanU = cumulative_sum ./ t;
    mean_U_squared = cumulative_sum_sq ./ t;
    varPopU = mean_U_squared - meanU.^2;
end
function [data_cell, var_names, file_info] = read_tecplot_plt(filename)
% read_tecplot_plt_binary Reads binary Tecplot .plt files created by 'liton_ordered_tec'.
%   This function is specifically designed to parse the binary format
%   written by the user's provided TEC_FILE class.
%
%   INPUT:
%       filename - A string containing the name of the binary .plt file.
%
%   OUTPUT:
%       data_cell - A 1xN cell array where N is the number of variables.
%                   Each cell contains a column vector of data for one variable.
%       var_names - A cell array of strings containing the variable names.
%       file_info - A struct containing metadata like Title, ZoneName, etc.

% --- 1. File Opening and Setup ---
fid = fopen(filename, 'rb'); % 'rb' for "read binary"
if fid == -1
    error('read_tecplot_plt_binary:CannotOpenFile', 'Cannot open file: %s', filename);
end
% Ensure the file is closed automatically when the function exits
cleanupObj = onCleanup(@() fclose(fid));

file_info = struct();

% --- 2. Read File Header ---
% i. Magic number, Version number
magic = fread(fid, 8, '*char')';
if ~strcmp(magic, '#!TDV112')
    error('read_tecplot_plt_binary:InvalidMagic', 'File is not a valid Tecplot binary file. Magic number is incorrect.');
end
file_info.Magic = magic;

% ii. Byte order check
fread(fid, 1, 'int32'); 

% iii. Title and variable names
file_info.FileType = fread(fid, 1, 'int32');
file_info.Title = read_null_terminated_string(fid);
num_vars = fread(fid, 1, 'int32');
var_names = cell(1, num_vars);
for i = 1:num_vars
    var_names{i} = read_null_terminated_string(fid);
end
file_info.Variables = var_names;

% --- 3. Read Zone Header ---
% This format only supports one zone, so we don't need a loop.
zone_marker = fread(fid, 1, 'float32');
if abs(zone_marker - 299.0) > 1e-6
    warning('read_tecplot_plt_binary:ZoneMarker', 'Zone marker is not 299.0 as expected.');
end

file_info.ZoneName = read_null_terminated_string(fid);
fread(fid, 1, 'int32'); % ParentZone
file_info.StrandId = fread(fid, 1, 'int32');
file_info.SolutionTime = fread(fid, 1, 'float64');
fread(fid, 1, 'int32'); % Not used
file_info.ZoneType = fread(fid, 1, 'int32');
fread(fid, 3, 'int32'); % VarLocation, RawFaceNeighbors, NumMiscFaceConnections

% Read IMax, JMax, KMax
zone_dims = fread(fid, 3, 'int32');
IMax = zone_dims(1);
JMax = zone_dims(2);
KMax = zone_dims(3);
num_points = IMax * JMax * KMax; % Total data points per variable

% Skip Auxiliary Data as it's not used in the user's case for reading
% But the writer supports it, so we must be able to skip it.
% The writer code writes "int32(1)" to indicate an aux pair, and "int32(0)" to end.
while true
    indicator = fread(fid, 1, 'int32');
    if indicator == 0 % No more aux data
        break;
    elseif indicator == 1 % Aux data pair follows
        read_null_terminated_string(fid); % Read and discard name
        fread(fid, 1, 'int32'); % Read and discard value format
        read_null_terminated_string(fid); % Read and discard value
    else
        % We've likely misread the stream, move back and assume end of aux
        fseek(fid, -4, 'cof');
        break;
    end
end


% --- 4. Read End of Header Marker ---
eoh_marker = fread(fid, 1, 'float32');
if abs(eoh_marker - 357.0) > 1e-6
    error('read_tecplot_plt_binary:EOHMarker', 'End of Header marker is not 357.0. Header parsing failed.');
end

% --- 5. Read Data Section ---
% i. Zone Data Header
data_zone_marker = fread(fid, 1, 'float32');
if abs(data_zone_marker - 299.0) > 1e-6
    warning('read_tecplot_plt_binary:DataZoneMarker', 'Data section zone marker is not 299.0.');
end

% ii. Read data types for each variable
% From the writer's gettype function: 1=single, 2=double, 3=int32, etc.
type_map = {'single', 'double', 'int32', 'int16', 'int8', 'bit1'};
data_types_int = zeros(1, num_vars);
data_types_str = cell(1, num_vars);
for i = 1:num_vars
    type_code = fread(fid, 1, 'int32');
    data_types_int(i) = type_code;
    data_types_str{i} = type_map{type_code};
end

% iii. Skip passive/shared variables info
fread(fid, 3, 'int32'); % Passive, Shared, Conn-Share

% iv. Skip Min/Max values (read 2 doubles for each variable)
fseek(fid, 2 * 8 * num_vars, 'cof');

% v. Read the actual zone data
disp(['Reading ', num2str(num_points), ' data points for each of the ', num2str(num_vars), ' variables.']);
data_cell = cell(1, num_vars);
for i = 1:num_vars
    % The '*' prefix tells fread to return data in its native class
    precision = ['*', data_types_str{i}];
    data_cell{i} = fread(fid, num_points, precision);
end

end


% --- Helper function to read null-terminated strings ---
function str = read_null_terminated_string(fid)
    % Reads characters until a null character (ASCII 0) is encountered.
    str_chars = [];
    while true
        char_val = fread(fid, 1, 'int32'); % The writer uses int32 for characters
        if isempty(char_val) || char_val == 0
            break;
        end
        str_chars(end+1) = char_val;
    end
    str = char(str_chars);
end
