clear; close all; clc;

%% basic settings
fileNumStart=2001; % This is not used in a loop, but kept for context
fileNumEnd=10000;
fileNumInterval=1;
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
%% ==================== Calculation of instantaneous fields ====================
% We will process a single snapshot
fileNum = fileNumEnd;
bin_filename = fullfile(inputDir, [namebase, num2str(fileNum),'.bin']);
disp(['Loading snapshot file: ', [namebase, num2str(fileNum),'.bin']]);

% Check if the binary file exists
if ~isfile(bin_filename)
    error('Binary file not found: %s', bin_filename);
end

[U,V,T,~] = readBinaryFile(bin_filename,nx,ny);
U = reshape(U,nx,ny);
V = reshape(V,nx,ny);
T = reshape(T,nx,ny);

% Non-dimensionalize
U_nd=U/params.velocityUnit;
V_nd=V/params.velocityUnit;

% Calculate buoyancy-related vorticity
[vor_buoyancy,~,~,~]=GRAD1(T,T,params.dx,params.dy);

disp('Instantaneous fields calculated.');

%% ==================== Step 2: Determine and Apply Thresholds for Both Signs ====================

% Define the absolute value of the threshold
threshold_magnitude = 0.1; % <--- Adjust the threshold magnitude here

% Create a binary mask for positive (counter-clockwise) vortices
binary_mask_pos = vor_buoyancy > threshold_magnitude;

% Create a binary mask for negative (clockwise) vortices
binary_mask_neg = vor_buoyancy < -threshold_magnitude;

% Pre-calculate grid properties needed for analysis
[Cx, Cy] = ndgrid(params.xGrid(1:nx), params.yGrid(1:ny));
area_weights = calculate_node_area_weights(params.xGrid, params.yGrid);

%% ==================== Analysis of Positive Vortices ====================
disp(' '); % Add a space for readability
disp('--- Analyzing Positive Vortices (Counter-Clockwise) ---');

% Step 1 & 2: Perform initial analysis to find area and centroid of each region
[stats_pos, num_vortices_pos] = custom_region_analysis(binary_mask_pos, params.xGrid(1:nx), params.yGrid(1:ny), area_weights);
disp(['Found ', num2str(num_vortices_pos), ' distinct positive vortices.']);

relabeled_mask_pos = zeros(size(binary_mask_pos));
lsc_center_pos = [NaN, NaN]; % Default value

if num_vortices_pos > 0
    % Step 3: Extract areas and sort them in descending order
    areas_pos = [stats_pos.Area];
    [~, sort_indices_pos] = sort(areas_pos, 'descend');

    % Step 4: Relabel the mask according to size
    [labeled_mask_pos, ~] = custom_labeling(binary_mask_pos);
    
    for i = 1:num_vortices_pos
        old_label = sort_indices_pos(i); % The original label of the i-th largest vortex
        new_label = i;                  % The new, sorted label
        relabeled_mask_pos(labeled_mask_pos == old_label) = new_label;
    end
    disp('Positive vortices sorted and relabeled by size.');

    % Extract the center of the largest positive vortex
    lsc_original_label_pos = sort_indices_pos(1);
    lsc_center_pos = stats_pos(lsc_original_label_pos).Centroid;
    disp(['Largest positive vortex center at (X, Y) = (', ...
          num2str(lsc_center_pos(1)), ', ', num2str(lsc_center_pos(2)), ')']);

%     % --- Visualization for Positive Vortices ---
%     figure;
%     custom_colormap_pos = [1 1 1; jet(num_vortices_pos)];
%     contourf(Cx, Cy, relabeled_mask_pos', 0:num_vortices_pos, 'LineStyle', 'none');
%     
%     hold on;
%     plot(lsc_center_pos(1), lsc_center_pos(2), 'p', ...
%         'MarkerEdgeColor','k', 'MarkerFaceColor','y', 'MarkerSize',15);
%     hold off;
%     
%     axis equal tight;
%     colormap(custom_colormap_pos);
%     c = colorbar('Ticks', (1:num_vortices_pos) - 0.5, 'TickLabels', 1:num_vortices_pos);
%     c.Label.String = 'Vortex Label (Sorted by Size)';
%     title({'Sorted Positive (CCW) Vortices', 'with Largest Center Marked'});
%     xlabel('x'); ylabel('y');
end

%% ==================== Analysis of Negative Vortices ====================
disp(' ');
disp('--- Analyzing Negative Vortices (Clockwise) ---');

% Step 1 & 2: Perform initial analysis
[stats_neg, num_vortices_neg] = custom_region_analysis(binary_mask_neg, params.xGrid(1:nx), params.yGrid(1:ny), area_weights);
disp(['Found ', num2str(num_vortices_neg), ' distinct negative vortices.']);

relabeled_mask_neg = zeros(size(binary_mask_neg));
lsc_center_neg = [NaN, NaN]; % Default value

if num_vortices_neg > 0
    % Step 3: Extract areas and sort
    areas_neg = [stats_neg.Area];
    [~, sort_indices_neg] = sort(areas_neg, 'descend');

    % Step 4: Relabel the mask
    [labeled_mask_neg, ~] = custom_labeling(binary_mask_neg);
    
    for i = 1:num_vortices_neg
        old_label = sort_indices_neg(i);
        new_label = i;
        relabeled_mask_neg(labeled_mask_neg == old_label) = new_label;
    end
    disp('Negative vortices sorted and relabeled by size.');

    % Extract the center of the largest negative vortex
    lsc_original_label_neg = sort_indices_neg(1);
    lsc_center_neg = stats_neg(lsc_original_label_neg).Centroid;
    disp(['Largest negative vortex center at (X, Y) = (', ...
          num2str(lsc_center_neg(1)), ', ', num2str(lsc_center_neg(2)), ')']);
          
%     % --- Visualization for Negative Vortices ---
%     figure;
%     custom_colormap_neg = [1 1 1; jet(num_vortices_neg)];
%     contourf(Cx, Cy, relabeled_mask_neg', 0:num_vortices_neg, 'LineStyle', 'none');
%     
%     hold on;
%     plot(lsc_center_neg(1), lsc_center_neg(2), 'p', ...
%         'MarkerEdgeColor','k', 'MarkerFaceColor','y', 'MarkerSize',15);
%     hold off;
%     
%     axis equal tight;
%     colormap(custom_colormap_neg);
%     c = colorbar('Ticks', (1:num_vortices_neg) - 0.5, 'TickLabels', 1:num_vortices_neg);
%     c.Label.String = 'Vortex Label (Sorted by Size)';
%     title({'Sorted Negative (CW) Vortices', 'with Largest Center Marked'});
%     xlabel('x'); ylabel('y');
end


%% ==================== Save All Results ====================
disp(' ');
disp('Saving all result files...');

% --- Save the original continuous field (unchanged) ---
tec_file = liton_ordered_tec.TEC_FILE;
tec_file.FileName = 'vor_buoyancy';
tec_file.Variables = {'x','y','vor_buoyancy'};
tec_file.Zones = liton_ordered_tec.TEC_ZONE;
tec_file.Zones.Data = {Cx,Cy,vor_buoyancy};
tec_file.write_plt();

% --- Save Positive Vortex Results ---
if num_vortices_pos > 0
    % Save binary mask for positive vortices
    tec_file = liton_ordered_tec.TEC_FILE;
    tec_file.FileName = 'vor_buoyancy_inden_pos';
    tec_file.Variables = {'x','y','vor_buoyancy_inden_pos'};
    tec_file.Zones = liton_ordered_tec.TEC_ZONE;
    tec_file.Zones.Data = {Cx,Cy,double(binary_mask_pos)};
    tec_file.write_plt();
    
    % Save sorted labels for positive vortices
    tec_file = liton_ordered_tec.TEC_FILE;
    tec_file.FileName = 'vortexBuoyancy_sorted_labels_pos';
    tec_file.Variables = {'x','y','label_pos'};
    tec_file.Zones = liton_ordered_tec.TEC_ZONE;
    tec_file.Zones.Data = {Cx,Cy,relabeled_mask_pos};
    tec_file.write_plt();

    % Save center coordinates of the largest positive vortex
    fid = fopen('lsc_center_coords_pos.txt', 'w');
    fprintf(fid, '%.8f, %.8f\n', lsc_center_pos(1), lsc_center_pos(2));
    fclose(fid);
    disp('Positive vortex files saved.');
end

% --- Save Negative Vortex Results ---
if num_vortices_neg > 0
    % Save binary mask for negative vortices
    tec_file = liton_ordered_tec.TEC_FILE;
    tec_file.FileName = 'vor_buoyancy_inden_neg';
    tec_file.Variables = {'x','y','vor_buoyancy_inden_neg'};
    tec_file.Zones = liton_ordered_tec.TEC_ZONE;
    tec_file.Zones.Data = {Cx,Cy,double(binary_mask_neg)};
    tec_file.write_plt();
    
    % Save sorted labels for negative vortices
    tec_file = liton_ordered_tec.TEC_FILE;
    tec_file.FileName = 'vortexBuoyancy_sorted_labels_neg';
    tec_file.Variables = {'x','y','label_neg'};
    tec_file.Zones = liton_ordered_tec.TEC_ZONE;
    tec_file.Zones.Data = {Cx,Cy,relabeled_mask_neg};
    tec_file.write_plt();

    % Save center coordinates of the largest negative vortex
    fid = fopen('lsc_center_coords_neg.txt', 'w');
    fprintf(fid, '%.8f, %.8f\n', lsc_center_neg(1), lsc_center_neg(2));
    fclose(fid);
    disp('Negative vortex files saved.');
end

disp('All files saved successfully.');

%%
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
function [node_area_weights] = calculate_node_area_weights(x_node_coords, y_node_coords)
    nx_nodes = length(x_node_coords)-1;
    ny_nodes = length(y_node_coords)-1;

    dx_contrib = zeros(nx_nodes, 1);
    if nx_nodes == 1
        dx_contrib(1) = 1;
    else
        dx_contrib(1) = (x_node_coords(2) - 0) / 2;
        for i = 2:(nx_nodes)
            dx_contrib(i) = (x_node_coords(i+1) - x_node_coords(i-1)) / 2;
        end
    end

    dy_contrib = zeros(ny_nodes, 1);
    if ny_nodes == 1
        dy_contrib(1) = 1;
    else
        dy_contrib(1) = (y_node_coords(2) - 0) / 2;
        for i = 2:(ny_nodes)
            dy_contrib(i) = (y_node_coords(i+1) - y_node_coords(i-1)) / 2;
        end
    end
    
    if any(dx_contrib < 0); error('Negative dx_contrib calculated.'); end
    if any(dy_contrib < 0); error('Negative dy_contrib calculated.'); end

    node_area_weights = dx_contrib * dy_contrib';
end

function [labeled_mask, label_count] = custom_labeling(binary_mask)
    [rows, cols] = size(binary_mask);
    visited = false(rows, cols);
    labeled_mask = zeros(rows, cols);
    label_count = 0;

    for r = 1:rows
        for c = 1:cols
            if binary_mask(r, c) && ~visited(r, c)
                label_count = label_count + 1;
                
                stack = [r, c];
                visited(r, c) = true;
                labeled_mask(r, c) = label_count;
                
                while ~isempty(stack)
                    curr = stack(end, :);
                    stack(end, :) = [];
                    curr_r = curr(1);
                    curr_c = curr(2);
                    
                    for dr = -1:1
                        for dc = -1:1
                            if dr == 0 && dc == 0; continue; end
                            nr = curr_r + dr;
                            nc = curr_c + dc;
                            
                            if nr >= 1 && nr <= rows && nc >= 1 && nc <= cols && ...
                               binary_mask(nr, nc) && ~visited(nr, nc)
                                visited(nr, nc) = true;
                                labeled_mask(nr, nc) = label_count;
                                stack(end+1, :) = [nr, nc];
                            end
                        end
                    end
                end
            end
        end
    end
end

function [stats, label_count] = custom_region_analysis(binary_mask, xGrid, yGrid, area_weights)
    [rows, cols] = size(binary_mask);
    visited = false(rows, cols);
    stats = struct('Area', {}, 'Centroid', {});
    label_count = 0;

    for r = 1:rows
        for c = 1:cols
            if binary_mask(r, c) && ~visited(r, c)
                label_count = label_count + 1;
                
                physical_area = 0;
                moment_x = 0;
                moment_y = 0;
                
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