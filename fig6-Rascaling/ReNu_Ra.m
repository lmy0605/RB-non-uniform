clear; close all; clc;

%% --- Configuration for Supplementary Data ---
fileNumStart = 1001;
fileNumEnd = 10000;
inputDir = 'C:\Users\user\OneDrive\Desktop\01-RB-0612\fig3-timeseries-0\1e8\timeseries_ReNu.plt';
inputDir2 = 'C:\Users\user\OneDrive\Desktop\01-RB-0612\fig4_timeseries-stationary\1e8\timeseries_Nu.plt';

[Revol_8, ~] = read_tecplot_plt(inputDir);
[Nu_8, ~] = read_tecplot_plt(inputDir2);
Re_8=mean(Revol_8{2}(fileNumStart:fileNumEnd));
Nuwall_8=mean(Revol_8{3}(fileNumStart:fileNumEnd));
%% --- Configuration for Supplementary Data ---
fileNumStart = 1501;
fileNumEnd = 10000;
inputDir = 'C:\Users\user\OneDrive\Desktop\01-RB-0612\fig3-timeseries-0\1e9\timeseries_ReNu.plt';
inputDir2 = 'C:\Users\user\OneDrive\Desktop\01-RB-0612\fig4_timeseries-stationary\1e9\timeseries_Nu.plt';

[Revol_9, ~] = read_tecplot_plt(inputDir);
[Nu_9, ~] = read_tecplot_plt(inputDir2);
Re_9=mean(Revol_9{2}(fileNumStart:fileNumEnd));
Nuwall_9=mean(Revol_9{3}(fileNumStart:fileNumEnd));
%% --- Configuration for Supplementary Data ---
fileNumStart = 2001;
fileNumEnd = 10000;
inputDir = 'C:\Users\user\OneDrive\Desktop\01-RB-0612\fig3-timeseries-0\1e10\timeseries_ReNu.plt';
inputDir2 = 'C:\Users\user\OneDrive\Desktop\01-RB-0612\fig4_timeseries-stationary\1e10\timeseries_Nu.plt';

[Revol_10, ~] = read_tecplot_plt(inputDir);
[Nu_10, ~] = read_tecplot_plt(inputDir2);
Re_10=mean(Revol_10{2}(fileNumStart:fileNumEnd));
Nuwall_10=mean(Revol_10{3}(fileNumStart:fileNumEnd));
%% --- Configuration for Supplementary Data ---
fileNumStart = 2001;
fileNumEnd = 5194;
inputDir = 'C:\Users\user\OneDrive\Desktop\01-RB-0612\fig3-timeseries-0\1e11\timeseries_ReNu_update.plt';
inputDir2 = 'C:\Users\user\OneDrive\Desktop\01-RB-0612\fig4_timeseries-stationary\1e11\timeseries_Nu_updated.plt.plt';

[Revol_11, ~] = read_tecplot_plt(inputDir);
[Nu_11, ~] = read_tecplot_plt(inputDir2);
Re_11=mean(Revol_11{2}(fileNumStart:fileNumEnd));
Nuwall_11=mean(Revol_11{3}(fileNumStart:fileNumEnd));
%% --- Configuration for Supplementary Data ---
fileNumStart = 2086;
fileNumEnd = 2609;
inputDir = 'C:\Users\user\OneDrive\Desktop\01-RB-0612\fig3-timeseries-0\1e12\timeseries_ReNu_update.plt';
inputDir2 = 'C:\Users\user\OneDrive\Desktop\01-RB-0612\fig4_timeseries-stationary\1e12\timeseries_Nu_updated.plt.plt';

[Revol_12, ~] = read_tecplot_plt(inputDir);
[Nu_12, ~] = read_tecplot_plt(inputDir2);
Re_12=mean(Revol_12{2}(fileNumStart:fileNumEnd));
Nuwall_12=mean(Revol_12{3}(fileNumStart:fileNumEnd));
%%
Ra=[1e8,1e9,1e10,1e11,1e12];
Re=[Re_8,Re_9,Re_10,Re_11,Re_12];
Nuwall=[Nuwall_8,Nuwall_9,Nuwall_10,Nuwall_11,Nuwall_12];
Nuvol=[25.2108,53.273,94.308,186.1131,394.1462];
Nueu=[24.7972,51.7653,87.0589,163.0573,247.715];
Nuet=[25.2346,52.8451,93.2809,181.778,369.672];

% %%
% Re_t8=Revol_8{2}(1001:10000)./Ra(1)^0.62;
% Re_t9=Revol_9{2}(1501:10000)./Ra(2)^0.62;
% Re_t10=Revol_10{2}(2001:10000)./Ra(3)^0.62;
% Re_t11=Revol_11{2}(2001:5194)./Ra(4)^0.62;
% Re_t12=Revol_12{2}(2086:2609)./Ra(5)^0.62;
% 
% Nuwall_t8=Revol_8{3}(1001:10000)./Ra(1)^0.3;
% Nuwall_t9=Revol_9{3}(1501:10000)./Ra(2)^0.3;
% Nuwall_t10=Revol_10{3}(2001:10000)./Ra(3)^0.3;
% Nuwall_t11=Revol_11{3}(2001:5194)./Ra(4)^0.3;
% Nuwall_t12=Revol_12{3}(2086:2609)./Ra(5)^0.3;
% 
% Nuvol_t8=Nu_8{2}./Ra(1)^0.29;
% Nuvol_t9=Nu_9{2}./Ra(2)^0.29;
% Nuvol_t10=Nu_10{2}./Ra(3)^0.29;
% Nuvol_t11=Nu_11{2}./Ra(4)^0.29;
% Nuvol_t12=Nu_12{2}./Ra(5)^0.29;

%% --- Power Law Fitting with Error Analysis ---
% We will perform a linear regression on the log-transformed data.
% The model is: log10(Y) = b1 + b2*log10(Ra)
% where b2 is the exponent 'alpha' and b1 is log10(C).
% The 'regress' function provides confidence intervals for the coefficients.

% Prepare data for fitting
data_to_fit = {Re, Nuwall, Nuvol, Nueu, Nuet};
data_names = {'Re', 'Nu_{wall}', 'Nu_{vol}', 'Nu_{eu}', 'Nu_{et}'};
results = struct(); % To store all fitting results
fit_lines_matrix = zeros(length(data_to_fit), length(Ra));
% Log-transform the independent variable (Ra) and create the design matrix X
% The first column of ones is for the intercept term (b1).
logRa = log10(Ra(:)); % Ensure it's a column vector
X = [ones(size(logRa)), logRa];

fprintf('--- Power Law Fitting Results (alpha with 95%% Confidence Interval) ---\n');

for i = 1:length(data_to_fit)
    % Get current dataset and log-transform it
    current_data = data_to_fit{i};
    logData = log10(current_data(:));
    
    % Perform linear regression
    % b: coefficients [intercept; slope]
    % bint: 95% confidence intervals for coefficients
    [b, bint] = regress(logData, X);
    
    % Extract fitting parameters
    intercept = b(1);
    alpha = b(2); % The slope is our exponent
    
    % The prefactor C = 10^intercept
    C = 10^intercept;
    
    % Calculate the error band for alpha from the confidence interval
    % error = (upper_bound - lower_bound) / 2
    alpha_error = (bint(2,2) - bint(2,1)) / 2;
    
    % Store results
    results(i).name = data_names{i};
    results(i).C = C;
    results(i).alpha = alpha;
    results(i).alpha_error = alpha_error;
    fit_lines_matrix(i, :) = results(i).C * Ra.^results(i).alpha;
    
    % Display results
    fprintf('Fit for %-8s: Y = %.4f * Ra^(%.4f ± %.4f)\n', ...
            results(i).name, results(i).C, results(i).alpha, results(i).alpha_error);
end
disp('---------------------------------------------------------------------');

%%
% --- Output cumulative mean and variance time series ---
tec_file = liton_ordered_tec.TEC_FILE;
tec_file.FileName = 'Ra_scaling';
tec_file.Variables = {'Ra','Re','Nuwal','Nuvol','Nueu','Nuet'};
tec_file.Zones = liton_ordered_tec.TEC_ZONE;
tec_file.Zones.Data = {Ra,Re,Nuwall,Nuvol,Nueu,Nuet};
tec_file = tec_file.write_plt();

tec_file = liton_ordered_tec.TEC_FILE;
tec_file.FileName = 'Ra_fitted';
tec_file.Variables = {'Ra','Re','Nuwal','Nuvol','Nueu','Nuet'};
tec_file.Zones = liton_ordered_tec.TEC_ZONE;
tec_file.Zones.Data = {Ra, ...
                       fit_lines_matrix(1, :), ...
                       fit_lines_matrix(2, :), ...
                       fit_lines_matrix(3, :), ...
                       fit_lines_matrix(4, :), ...
                       fit_lines_matrix(5, :)};
tec_file = tec_file.write_plt();

tec_file = liton_ordered_tec.TEC_FILE;
tec_file.FileName = 'zoomRa_scaling';
tec_file.Variables = {'Ra','Re','Nuwal','Nuvol','Nueu','Nuet'};
tec_file.Zones = liton_ordered_tec.TEC_ZONE;
tec_file.Zones.Data = {Ra,Re./Ra.^results(1).alpha,Nuwall./Ra.^results(2).alpha,Nuvol./Ra.^results(3).alpha,Nueu./Ra.^results(4).alpha,Nuet./Ra.^results(5).alpha};
tec_file = tec_file.write_plt();

Raz=linspace(1e6,1e13,5);
Rez=results(1).C*ones(1,5);
Nuwallz=results(2).C*ones(1,5);
Nuvolz=results(3).C*ones(1,5);
Nueuz=results(4).C*ones(1,5);
Nuetz=results(5).C*ones(1,5);
tec_file = liton_ordered_tec.TEC_FILE;
tec_file.FileName = 'zoomRa_fitted';
tec_file.Variables = {'Ra','Re','Nuwal','Nuvol','Nueu','Nuet'};
tec_file.Zones = liton_ordered_tec.TEC_ZONE;
tec_file.Zones.Data = {Raz, ...
                       Rez, ...
                       Nuwallz, ...
                       Nuvolz, ...
                       Nueuz, ...
                       Nuetz};
tec_file = tec_file.write_plt();

% tec_file = liton_ordered_tec.TEC_FILE;
% tec_file.FileName = 'Ra_t8';
% tec_file.Variables = {'t','ReRa8','NuwallRa8','NuvolRa8'};
% tec_file.Zones = liton_ordered_tec.TEC_ZONE;
% tec_file.Zones.Data = {Nu_8{1},Re_t8,Nuwall_t8,Nuvol_t8};
% tec_file = tec_file.write_plt();
% 
% tec_file = liton_ordered_tec.TEC_FILE;
% tec_file.FileName = 'Ra_t9';
% tec_file.Variables = {'t','ReRa9','NuwallRa9','NuvolRa9'};
% tec_file.Zones = liton_ordered_tec.TEC_ZONE;
% tec_file.Zones.Data = {Nu_9{1},Re_t9,Nuwall_t9,Nuvol_t9};
% tec_file = tec_file.write_plt();
% 
% tec_file = liton_ordered_tec.TEC_FILE;
% tec_file.FileName = 'Ra_t10';
% tec_file.Variables = {'t','ReRa10','NuwallRa10','NuvolRa10'};
% tec_file.Zones = liton_ordered_tec.TEC_ZONE;
% tec_file.Zones.Data = {Nu_10{1},Re_t10,Nuwall_t10,Nuvol_t10};
% tec_file = tec_file.write_plt();
% 
% tec_file = liton_ordered_tec.TEC_FILE;
% tec_file.FileName = 'Ra_t11';
% tec_file.Variables = {'t','ReRa11','NuwallRa11','NuvolRa11'};
% tec_file.Zones = liton_ordered_tec.TEC_ZONE;
% tec_file.Zones.Data = {Nu_11{1},Re_t11,Nuwall_t11,Nuvol_t11};
% tec_file = tec_file.write_plt();
% 
% tec_file = liton_ordered_tec.TEC_FILE;
% tec_file.FileName = 'Ra_t12';
% tec_file.Variables = {'t','ReRa12','NuwallRa12','NuvolRa12'};
% tec_file.Zones = liton_ordered_tec.TEC_ZONE;
% tec_file.Zones.Data = {Nu_12{1},Re_t12,Nuwall_t12,Nuvol_t12};
% tec_file = tec_file.write_plt();

%%
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