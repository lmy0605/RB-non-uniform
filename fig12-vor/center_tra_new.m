clear; close all; clc;

%% ==================== Basic Settings ====================
fileNumStart=1501; % This is not used in a loop, but kept for context
fileNumEnd=10000;
fileNumInterval=1;
fileSum=fileNumEnd-fileNumStart+1;
inputDir = '/nfsdata4/AXu/RB-non-uniform/Ra1e9-mesh513/binFile-1-10000/'; 
namebase = 'buoyancyCavity-';
casename='1e9'; 

% --- Check if input directory exists before proceeding ---
disp('Verifying input directory...');
if ~isfolder(inputDir)
    error('Input directory not found: %s', inputDir);
end
disp(['Input directory found: ', inputDir]);

nx=513;
ny=nx;
constA=2.1;
Rayleigh=1e9;
Prandtl=0.71;


params = calculateSystemParameters(nx,ny, Rayleigh, Prandtl,constA,'log.log');
viscosity=sqrt(Prandtl/Rayleigh);

%% ==================== Pre-loop Setup ====================
% --- Pre-calculate constant values BEFORE the loop to improve performance ---
area_weights = calculate_node_area_weights(params.xGrid, params.yGrid);
[X_coords, Y_coords] = ndgrid(params.xGrid(1:nx), params.yGrid(1:ny));
% Define the sub-region for direction check
idx_range = floor(0.2*nx):floor(0.8*nx);
area_weights_direction = area_weights(idx_range, idx_range);

% --- Initialize storage for LSC center trajectory ---
lsc_centers_over_time = zeros(fileSum, 2);

%% ==================== Main Time Loop ====================
disp(['Starting to process ', num2str(fileSum), ' files...']);
% To use parfor, ensure analyze_lsc_ipt and other functions are accessible
% and handle the sliced variable 'lsc_centers_over_time' correctly.
for fileNum = fileNumStart:fileNumInterval:fileNumEnd
    t = fileNum - fileNumStart + 1;
    if(mod(fileNum,100)==0)
        disp(['Current data file is ', [namebase, num2str(fileNum),'.bin']]);
    end

    [U,V,~,~] = readBinaryFile(fullfile(inputDir, [namebase, num2str(fileNum),'.bin']),nx,ny);
    U = reshape(U,nx,ny);
    V = reshape(V,nx,ny);
    U = U/params.velocityUnit;
    V = V/params.velocityUnit;

    [~,~,VX,VY] = GRAD1(U,V,params.dx,params.dy); % We only need VX and UY for vor_z
    [~,UY,~,~] = GRAD1(U,V,params.dx,params.dy);
    vor_z = (VX-UY);

    % --- Step 1: Automatically Identify LSC Direction (in central region) ---
    vor_z_direction = vor_z(idx_range, idx_range);
    
    pos_mask = vor_z_direction > 0;
    neg_mask = vor_z_direction < 0;
    
    positive_enstrophy = sum(vor_z_direction(pos_mask).^2 .* area_weights_direction(pos_mask));
    negative_enstrophy = sum(vor_z_direction(neg_mask).^2 .* area_weights_direction(neg_mask));

    is_LSC_clockwise = negative_enstrophy > positive_enstrophy;
    
    % --- Step 2: Determine and Apply Threshold ---
    threshold_magnitude = 2;
    if is_LSC_clockwise
        binary_mask = vor_z < -threshold_magnitude;
    else
        binary_mask = vor_z > threshold_magnitude;
    end
    
    % --- Step 3: Analyze Vortices and find the LSC center ---
    lsc_center = analyze_lsc_ipt(binary_mask, X_coords, Y_coords, area_weights);
    
    % Store the result for this time step
    lsc_centers_over_time(t, :) = lsc_center;
end

%% ==================== Save and Plot Trajectory ====================
% Filter out any invalid points (NaNs) before saving
trajectory_data = lsc_centers_over_time(~any(isnan(lsc_centers_over_time), 2), :);

if isempty(trajectory_data)
    warning('No valid LSC centers were found in any of the files.');
else
    disp('Saving LSC trajectory to Tecplot file...');
    tec_file = liton_ordered_tec.TEC_FILE;
    tec_file.FileName = 'lsc_trajectory';
    tec_file.Variables = { 'X_center', 'Y_center'};
    tec_file.Zones = liton_ordered_tec.TEC_ZONE;
    % Corrected data order: X is column 1, Y is column 2
    tec_file.Zones.Data = {trajectory_data(:,1), trajectory_data(:,2)};
    tec_file.write_plt();
    disp('LSC trajectory saved successfully to lsc_trajectory.plt');
end

%% ======================== HELPER FUNCTIONS ========================

function lsc_center = analyze_lsc_ipt(binary_mask, X_coords, Y_coords, area_weights)
    % This function finds the LSC by identifying the largest vortex within a central region.
    % It returns the LSC centroid or [NaN, NaN] if not found.
    
    lsc_center = [NaN, NaN]; % Default output

    if ~any(binary_mask, 'all'); return; end

    % 1. Label all vortices in the full domain
    [labeled_mask, num_vortices] = bwlabel(binary_mask, 8);
    if num_vortices == 0; return; end

    % 2. Get pixel indices for each vortex
    stats = regionprops(labeled_mask, 'PixelIdxList');

    % 3. Calculate physical properties for ALL vortices
    physical_centroids = zeros(num_vortices, 2);
    physical_areas = zeros(num_vortices, 1);
    
    for k = 1:num_vortices
        indices_k = stats(k).PixelIdxList;
        weights_k = area_weights(indices_k);
        total_area = sum(weights_k);
        
        physical_areas(k) = total_area;
        
        if total_area > 0
            moment_x = sum(X_coords(indices_k) .* weights_k);
            moment_y = sum(Y_coords(indices_k) .* weights_k);
            physical_centroids(k, :) = [moment_x / total_area, moment_y / total_area];
        else
            physical_centroids(k, :) = [NaN, NaN];
        end
    end

    % 4. Filter to find which of these vortices are in the central region
    is_in_center = physical_centroids(:, 1) > 0.25 & physical_centroids(:, 1) < 0.75 & ...
                   physical_centroids(:, 2) > 0.25 & physical_centroids(:, 2) < 0.75;
    
    central_vortex_labels = find(is_in_center);

    % 5. Among the central vortices, find the one with the largest area
    if ~isempty(central_vortex_labels)
        areas_of_central_vortices = physical_areas(central_vortex_labels);
        [~, max_idx] = max(areas_of_central_vortices);
        
        % The label of the LSC is the 'max_idx'-th element of our filtered list
        lsc_label = central_vortex_labels(max_idx);
        
        % The final LSC center is the centroid corresponding to this label
        lsc_center = physical_centroids(lsc_label, :);
    end
end

function [U, V ,T, rho] = readBinaryFile(file, nx, ny)
fid = fopen(file,'r');

[~,~] = fread(fid,1,'int32'); % Read one tag at the beginning
U = fread(fid,nx*ny,'float64');

[~,~] = fread(fid,2,'int32'); % Read two tags...
V = fread(fid,nx*ny,'float64');

[~,~] = fread(fid,2,'int32'); % Read two tags...
T = fread(fid,nx*ny,'float64');

[~,~] = fread(fid,2,'int32'); % Read two tags...
rho = fread(fid,nx*ny,'float64');

[~,~] = fread(fid,1,'int32'); % Read one tag at the end

fclose(fid);
end

%nonUniform Weight
function [U_dxWeighted,U_dyWeighted,U_areaWeighted]=nonUniformAverage(U,xGrid,yGrid)
[nx, ny] = size(U);
dx = zeros(1,nx);
dy = zeros(1,ny);

dx(1) = (xGrid(2) -0)/2; 
for i = 2:nx
    dx(i) = (xGrid(i+1) - xGrid(i-1))/2;
end

dy(1) = (yGrid(2) - 0)/2; 
for j = 2:ny
    dy(j) = (yGrid(j+1) - yGrid(j-1))/2; 
end

U_dxWeighted = zeros(nx,ny);
U_dyWeighted = zeros(nx,ny);
U_areaWeighted = zeros(nx,ny);

for j = 1:ny
    for i = 1:nx
        U_dxWeighted(i,j) = U(i,j) .* dx(i);      
        U_dyWeighted(i,j) = U(i,j) .* dy(j);      
        U_areaWeighted(i,j) = U(i,j) .* dx(i) .* dy(j);
    end
end
end

%wall derivation
function [deri_Ty_bottom,deri_Ty_top]=deri_Twall(T,T_bottom,T_top,length_LB,q,nx,ny)
deri_Ty_bottom=zeros(nx,1);
deri_Ty_top=zeros(nx,1);
for i=1:1:nx
    deri_Ty_bottom(i)=(-4*q*(q+1)*T_bottom+(2*q+1)^2*T(i,1)-T(i,2))/(q*(2*q+1))*length_LB;
    deri_Ty_top(i)=(4*q*(q+1)*T_top-(2*q+1)^2*T(i,ny)+T(i,ny-1))/(q*(2*q+1))*length_LB;
end
end

function [GRAD_UX,GRAD_UY,GRAD_VX,GRAD_VY]=GRAD1(U,V,dx,dy)
%center dx=(dx(i)^2*(f(i+1)-f(i))+dx(i+1)^2*(f(i)-f(i-1)))/(dx(i)*dx(i+1)*(dx(i)+dx(i+1)))
%forward dx=(-(2*dx(i+1)*dx(i+2)+dx(i+2)^2)*f(i)+(dx(i+1)+dx(i+2))^2*f(i+1)-dx(i+1)^2*f(i+2))/(dx(i+2)*dx(i+1)*(dx(i+2)+dx(i+1)))
%backward dx=((2*dx(i)*dx(i-1)+dx(i-1)^2)*f(i)-(dx(i)+dx(i-1))^2*f(i-1)+dx(i)^2*f(i-2))/(dx(i)*dx(i-1)*(dx(i)+dx(i-1)))
    [X,Y]=size(U);
    GRAD_UX= zeros(X, Y);
    GRAD_UY= zeros(X, Y);
    GRAD_VX= zeros(X, Y);
    GRAD_VY= zeros(X, Y);
    for j = 1:1:Y
        for i = 1:1:X
            if i==1
                GRAD_UX(i,j) = (-(2*dx(i+1)*dx(i+2)+dx(i+2)^2)*U(i,j) + (dx(i+1)+dx(i+2))^2*U(i+1,j) - dx(i+1)^2*U(i+2,j)) / (dx(i+2)*dx(i+1)*(dx(i+2)+dx(i+1)));
                GRAD_VX(i,j) = (-(2*dx(i+1)*dx(i+2)+dx(i+2)^2)*V(i,j) + (dx(i+1)+dx(i+2))^2*V(i+1,j) - dx(i+1)^2*V(i+2,j)) / (dx(i+2)*dx(i+1)*(dx(i+2)+dx(i+1)));
            elseif i==X
                GRAD_UX(i,j)=((2*dx(i)*dx(i-1)+dx(i-1)^2)*U(i,j)-(dx(i)+dx(i-1))^2*U(i-1,j)+dx(i)^2*U(i-2,j))/(dx(i)*dx(i-1)*(dx(i)+dx(i-1)));
                GRAD_VX(i,j)=((2*dx(i)*dx(i-1)+dx(i-1)^2)*V(i,j)-(dx(i)+dx(i-1))^2*V(i-1,j)+dx(i)^2*V(i-2,j))/(dx(i)*dx(i-1)*(dx(i)+dx(i-1)));
            else
                GRAD_UX(i,j)=(dx(i)^2*(U(i+1,j)-U(i,j))+dx(i+1)^2*(U(i,j)-U(i-1,j)))/(dx(i)*dx(i+1)*(dx(i)+dx(i+1)));
                GRAD_VX(i,j)=(dx(i)^2*(V(i+1,j)-V(i,j))+dx(i+1)^2*(V(i,j)-V(i-1,j)))/(dx(i)*dx(i+1)*(dx(i)+dx(i+1)));
            end
            
            if j==1
                 GRAD_UY(i,j) = (-(2*dy(j+1)*dy(j+2)+dy(j+2)^2)*U(i,j) + (dy(j+1)+dy(j+2))^2*U(i,j+1) - dy(j+1)^2*U(i,j+2)) / (dy(j+2)*dy(j+1)*(dy(j+2)+dy(j+1)));
                GRAD_VY(i,j) = (-(2*dy(j+1)*dy(j+2)+dy(j+2)^2)*V(i,j) + (dy(j+1)+dy(j+2))^2*V(i,j+1) - dy(j+1)^2*V(i,j+2)) / (dy(j+2)*dy(j+1)*(dy(j+2)+dy(j+1)));
            elseif j==Y
                GRAD_UY(i,j)=((2*dy(j)*dy(j-1)+dy(j-1)^2)*U(i,j)-(dy(j)+dy(j-1))^2*U(i,j-1)+dy(j)^2*U(i,j-2))/(dy(j)*dy(j-1)*(dy(j)+dy(j-1)));
                GRAD_VY(i,j)=((2*dy(j)*dy(j-1)+dy(j-1)^2)*V(i,j)-(dy(j)+dy(j-1))^2*V(i,j-1)+dy(j)^2*V(i,j-2))/(dy(j)*dy(j-1)*(dy(j)+dy(j-1)));
            else
                GRAD_UY(i,j)=(dy(j)^2*(U(i,j+1)-U(i,j))+dy(j+1)^2*(U(i,j)-U(i,j-1)))/(dy(j)*dy(j+1)*(dy(j)+dy(j+1)));
                GRAD_VY(i,j)=(dy(j)^2*(V(i,j+1)-V(i,j))+dy(j+1)^2*(V(i,j)-V(i,j-1)))/(dy(j)*dy(j+1)*(dy(j)+dy(j+1)));
            end
        end
    end
end

% --- calculate_node_area_weights (Defined above, ensure it's here or accessible) ---
% function [node_area_weights] = calculate_node_area_weights(x_node_coords, y_node_coords) ... end
function [node_area_weights] = calculate_node_area_weights(x_node_coords, y_node_coords)
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

% function [labeled_mask, label_count] = custom_labeling(binary_mask)
%     [rows, cols] = size(binary_mask);
%     visited = false(rows, cols);
%     labeled_mask = zeros(rows, cols);
%     label_count = 0;
% 
%     for r = 1:rows
%         for c = 1:cols
%             if binary_mask(r, c) && ~visited(r, c)
%                 label_count = label_count + 1;
%                 
%                 stack = [r, c];
%                 visited(r, c) = true;
%                 labeled_mask(r, c) = label_count;
%                 
%                 while ~isempty(stack)
%                     curr = stack(end, :);
%                     stack(end, :) = [];
%                     curr_r = curr(1);
%                     curr_c = curr(2);
%                     
%                     for dr = -1:1
%                         for dc = -1:1
%                             if dr == 0 && dc == 0; continue; end
%                             nr = curr_r + dr;
%                             nc = curr_c + dc;
%                             
%                             if nr >= 1 && nr <= rows && nc >= 1 && nc <= cols && ...
%                                binary_mask(nr, nc) && ~visited(nr, nc)
%                                 visited(nr, nc) = true;
%                                 labeled_mask(nr, nc) = label_count;
%                                 stack(end+1, :) = [nr, nc];
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end

function [stats, label_count] = custom_region_analysis(binary_mask, xGrid, yGrid, area_weights)
    % This version uses pre-calculated area_weights for physical calculations.
    
    [rows, cols] = size(binary_mask);
    visited = false(rows, cols);
    stats = struct('Area', {}, 'Centroid', {});
    label_count = 0;

    for r = 1:rows
        for c = 1:cols
            if binary_mask(r, c) && ~visited(r, c)
                label_count = label_count + 1;
                
                physical_area = 0;
                moment_x = 0; % Sum of xi * Ai
                moment_y = 0; % Sum of yi * Ai
                
                stack = [r, c];
                visited(r, c) = true;
                
                while ~isempty(stack)
                    curr = stack(end, :);
                    stack(end, :) = [];
                    curr_r = curr(1);
                    curr_c = curr(2);
                    
                    area_cell = area_weights(curr_r, curr_c);
                    x_cell = xGrid(curr_c);
                    y_cell = yGrid(curr_r);
                    
                    physical_area = physical_area + area_cell;
                    moment_x = moment_x + x_cell * area_cell;
                    moment_y = moment_y + y_cell * area_cell;
                    
                    % 8-connectivity DFS search
                    for dr = -1:1
                        for dc = -1:1
                            if dr == 0 && dc == 0; continue; end
                            nr = curr_r + dr;
                            nc = curr_c + dc;
                            
                            if nr >= 1 && nr <= rows && nc >= 1 && nc <= cols && ...
                               binary_mask(nr, nc) && ~visited(nr, nc)
                                visited(nr, nc) = true;
                                stack(end+1, :) = [nr, nc];
                            end
                        end
                    end
                end
                
                stats(label_count).Area = physical_area;
                if physical_area > 0
                    stats(label_count).Centroid = [moment_x / physical_area, moment_y / physical_area];
                else
                    stats(label_count).Centroid = [NaN, NaN];
                end
            end
        end
    end
end