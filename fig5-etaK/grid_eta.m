clear; close all; clc;

%% basic settings
casename='1e12';
inputDir_etaK = strcat(casename,'\etaK_timeAvg.plt'); % please rename data folder as "binFile"
inputDir_etaU = strcat(casename,'\etaU_timeAvg.plt');

nx=4097;
ny=nx;
constA=3.1;
Rayleigh=1e12;
Prandtl=0.71;

params = calculateSystemParameters(nx,ny, Rayleigh, Prandtl,constA,'log.log');
viscosity=sqrt(Prandtl/Rayleigh);
%%
[existing_data_cell_etaK, var_names_etaK] = read_tecplot_plt(inputDir_etaK);
[existing_data_cell_etaU, var_names_etaU] = read_tecplot_plt(inputDir_etaU);
etaK=reshape(existing_data_cell_etaK{3},nx,ny);
etaU=reshape(existing_data_cell_etaU{3},nx,ny);
Cx=reshape(existing_data_cell_etaU{1},nx,ny);
Cy=reshape(existing_data_cell_etaU{2},nx,ny);

[deltax,deltay,deltaxy]=calculate_node_area_weights(params.xGrid,params.yGrid);
[deltaX,deltaY]=ndgrid(deltax,deltay);

grid_resolution = max(deltaX, deltaY); 
grid_resolution2 = sqrt(deltaxy);

resolution_ratio_etaK = grid_resolution ./ etaK;
resolution_ratio_etaU = grid_resolution ./ etaU;
resolution_ratio_etaK2 = grid_resolution2 ./ etaK;
resolution_ratio_etaU2 = grid_resolution2 ./ etaU;
aspect_ratio=max(deltaX, deltaY)./min(deltaX, deltaY);

%%
tec_file = liton_ordered_tec.TEC_FILE;
tec_file.FileName = strcat('resolution_ratio_etaK_',casename);
tec_file.Variables = {'x','y','resolution_ratio_etaK','aspect_ratio'};
tec_file.Zones = liton_ordered_tec.TEC_ZONE;
tec_file.Zones.Data = {Cx,Cy,resolution_ratio_etaK,aspect_ratio};
tec_file = tec_file.write_plt();

tec_file = liton_ordered_tec.TEC_FILE;
tec_file.FileName = strcat('resolution_ratio_etaU_',casename);
tec_file.Variables = {'x','y','resolution_ratio_etaU','aspect_ratio'};
tec_file.Zones = liton_ordered_tec.TEC_ZONE;
tec_file.Zones.Data = {Cx,Cy,resolution_ratio_etaU,aspect_ratio};
tec_file = tec_file.write_plt();

tec_file = liton_ordered_tec.TEC_FILE;
tec_file.FileName = strcat('2resolution_ratio_etaK_',casename);
tec_file.Variables = {'x','y','resolution_ratio_etaK','aspect_ratio'};
tec_file.Zones = liton_ordered_tec.TEC_ZONE;
tec_file.Zones.Data = {Cx,Cy,resolution_ratio_etaK2,aspect_ratio};
tec_file = tec_file.write_plt();

tec_file = liton_ordered_tec.TEC_FILE;
tec_file.FileName = strcat('2resolution_ratio_etaU_',casename);
tec_file.Variables = {'x','y','resolution_ratio_etaU','aspect_ratio'};
tec_file.Zones = liton_ordered_tec.TEC_ZONE;
tec_file.Zones.Data = {Cx,Cy,resolution_ratio_etaU2,aspect_ratio};
tec_file = tec_file.write_plt();
%% --- Functions ---
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
    % (The rest of the function code is identical to your original script)
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
    % (The rest of the function code is identical to your original script)
    deri_Ty_bottom=zeros(nx,1); deri_Ty_top=zeros(nx,1);
    for i=1:1:nx
        deri_Ty_bottom(i)=(-4*q*(q+1)*T_bottom+(2*q+1)^2*T(i,1)-T(i,2))/(q*(2*q+1))*length_LB;
        deri_Ty_top(i)=(4*q*(q+1)*T_top-(2*q+1)^2*T(i,ny)+T(i,ny-1))/(q*(2*q+1))*length_LB;
    end
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

function [dx_contrib,dy_contrib,node_area_weights] = calculate_node_area_weights(x_node_coords, y_node_coords)
    % x_node_coords: 1D array of x-coordinates of nodes (length nx)
    % y_node_coords: 1D array of y-coordinates of nodes (length ny)
    % Returns: An nx x ny matrix of area weights for each node.

    nx_nodes = length(x_node_coords)-1;
    ny_nodes = length(y_node_coords)-1;

    dx_contrib = zeros(nx_nodes, 1); % Column vector for dx contributions
    if nx_nodes == 1
        warning('calculate_node_area_weights: Only one X node. Assuming unit contribution or check logic.');
        dx_contrib(1) = 1; % Or handle as an error, or require domain width
    else
        % Contribution from the first node
        dx_contrib(1) = (x_node_coords(2) - 0) / 2;
        % Contribution from internal nodes
        for i = 2:(nx_nodes)
            dx_contrib(i) = (x_node_coords(i+1) - x_node_coords(i-1)) / 2;
        end
    end

    dy_contrib = zeros(ny_nodes, 1); % Column vector for dy contributions
    if ny_nodes == 1
        warning('calculate_node_area_weights: Only one Y node. Assuming unit contribution or check logic.');
        dy_contrib(1) = 1; % Or handle as an error, or require domain height
    else
        dy_contrib(1) = (y_node_coords(2) - 0) / 2;
        for i = 2:(ny_nodes)
            dy_contrib(i) = (y_node_coords(i+1) - y_node_coords(i-1)) / 2;
        end
    end
    
    % Ensure no negative contributions (e.g., if grid coords are not monotonic)
    if any(dx_contrib < 0)
        error('Negative dx_contrib calculated. Check x_node_coords order and values.');
    end
    if any(dy_contrib < 0)
        error('Negative dy_contrib calculated. Check y_node_coords order and values.');
    end

    % Create the 2D area weights matrix
    % T is (nx, ny), so weights should be (nx, ny)
    % ndgrid(dx_contrib, dy_contrib) will produce DX(nx,ny), DY(nx,ny)
    % where DX(:,j) = dx_contrib and DY(i,:) = dy_contrib' (if dy_contrib is row)
    % or DX has columns as dx_contrib, DY has rows as dy_contrib
    % We need area_weights(i,j) = dx_contrib(i) * dy_contrib(j)
    
    % dx_contrib is (nx x 1), dy_contrib is (ny x 1)
    % We want node_area_weights(ix, iy) = dx_contrib(ix) * dy_contrib(iy)
    % This can be achieved by outer product: dx_contrib * dy_contrib'
    node_area_weights = dx_contrib * dy_contrib'; % Results in (nx x ny) matrix
    
    % If you prefer ndgrid:
    % [DX_C, DY_C] = ndgrid(dx_contrib, dy_contrib);
    % node_area_weights = DX_C .* DY_C; % This also gives nx x ny
end