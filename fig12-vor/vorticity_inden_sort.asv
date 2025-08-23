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
% We will process a single snapshot, as in your script
fileNum = 2110; % Process the last file
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

% Calculate velocity gradients and derived quantities
[UX,UY,VX,VY]=GRAD1(U_nd,V_nd,params.dx,params.dy);
vor_z=(VX-UY);

% % NOTE: Your original code passed T twice to GRAD1. Assuming you want dT/dx and dT/dy
% [vor_buoyancy,~,~,~]=GRAD1(T,T,params.dx,params.dy);
% 
% enstrophy_buoyancy = vor_z .* vor_buoyancy;

disp('Instantaneous fields calculated.');

%% ==================== Step 1: Automatically Identify LSC Direction ====================

% 将涡量场分为正部和负部
positive_vorticity = vor_z(vor_z > 0);
negative_vorticity = vor_z(vor_z < 0);

sum_positive_vor = length(positive_vorticity);
sum_negative_vor_abs = length(negative_vorticity);

% 比较哪个总和更大
if sum_positive_vor > sum_negative_vor_abs
    is_LSC_clockwise = false; % 逆时针 (正涡量) 占主导
    disp('LSC direction identified automatically: Counter-Clockwise (Positive Vorticity)');
else
    is_LSC_clockwise = true; % 顺时针 (负涡量) 占主导
    disp('LSC direction identified automatically: Clockwise (Negative Vorticity)');
end

%% ==================== Step 2: Determine and Apply Threshold ====================

% --- 交互式地确定阈值 ---
% 为了方便调整，把阈值的绝对值定义为一个变量
threshold_magnitude = 2; % <--- 在这里调整阈值的绝对值

if is_LSC_clockwise
    % 如果LSC是顺时针 (蓝色, vor_z < 0)
    threshold = -threshold_magnitude; 
    binary_mask = vor_z < threshold; % 小于阈值的点为1
else
    % 如果LSC是逆时针 (红色, vor_z > 0)
    threshold = threshold_magnitude;
    binary_mask = vor_z > threshold; % 大于阈值的点为1
end

%%
try
    [labeled_mask, num_vortices] = bwlabel(binary_mask, 8); % 8-connectivity
    disp(['Found ', num2str(num_vortices), ' distinct vortices.']);
catch ME
    if strcmp(ME.identifier, 'MATLAB:UndefinedFunction')
        disp('Warning: Image Processing Toolbox not found. `bwlabel` is unavailable.');
        disp('Will plot the simple binary mask without labeling vortices.');
        labeled_mask = double(binary_mask); % 如果没有工具箱，直接用0和1的矩阵
        num_vortices = -1; % 表示未进行标记
    else
        rethrow(ME);
    end
end

%% ==================== Vortex Labeling and Sorting (Toolbox-Free version) ====================
disp('Performing vortex analysis (Toolbox-Free)...');
area_weights=calculate_node_area_weights(params.xGrid,params.yGrid);
% Step 1 & 2: 使用自定义函数进行初始分析，得到面积和质心
[stats, num_vortices] = custom_region_analysis(binary_mask, params.xGrid(1:nx),params.yGrid(1:ny), area_weights);
disp(['Found ', num_vortices, ' distinct vortices.']);

if num_vortices > 0
    % Step 3: 提取面积并排序
    areas = [stats.Area];
    [sorted_areas, sort_indices] = sort(areas, 'descend');

    % Step 4: 重新标记。这一步比较复杂，需要再次遍历图像
    % 我们先创建一个映射关系：旧标签 -> 新标签
    % 注意：我们的自定义函数没有输出labeled_mask，所以需要重新生成
    [labeled_mask, ~] = custom_labeling(binary_mask); % 需要一个新的辅助函数
    
    relabeled_mask = zeros(size(labeled_mask));
    for i = 1:num_vortices
        old_label = sort_indices(i); % 第i大的涡旋的旧标签
        new_label = i; % 新的、有序的标签
        relabeled_mask(labeled_mask == old_label) = new_label;
    end

    disp('Vortices sorted and relabeled by size.');

      % --- Extract the center of the largest vortex (LSC) ---
    lsc_original_label = sort_indices(1);
    lsc_center = stats(lsc_original_label).Centroid;
    disp(['Largest vortex (LSC) center found at (X, Y) = (', ...
    num2str(lsc_center(1)), ', ', num2str(lsc_center(2)), ')']);
end
% %% ==================== 5. Visualization with LSC Center Marked ====================
% figure;
% if num_vortices > 0
%     [Cx, Cy] = ndgrid(params.xGrid(1:nx), params.yGrid(1:ny));
%     custom_colormap = [1 1 1; jet(num_vortices)];
%     contourf(Cx, Cy, relabeled_mask', 0:num_vortices, 'LineStyle', 'none');
%     
%     hold on;
%     % Mark the LSC center with a star
%     plot(lsc_center(1), lsc_center(2), 'p', ...
%         'MarkerEdgeColor','k', 'MarkerFaceColor','y', 'MarkerSize',15);
%     hold off;
%     
%     axis equal tight;
%     colormap(custom_colormap);
%     c = colorbar('Ticks', (1:num_vortices) - 0.5, 'TickLabels', 1:num_vortices);
%     c.Label.String = 'Vortex Label (Sorted by Size)';
%     title({'Sorted and Labeled Vortices', 'with LSC Center Marked'});
%     xlabel('x'); ylabel('y');
% end
%% ========== MODIFIED: Step 4: Create and Save Visualization without Displaying ==========
disp('Generating and saving visualization as PNG...');

% Create a figure, but keep it invisible
fig = figure('Visible', 'off'); 

if num_vortices > 0
    [Cx, Cy] = ndgrid(params.xGrid(1:nx), params.yGrid(1:ny));
    
    % Define a colormap: first color is white (for background), rest are for vortices
    custom_colormap = [1 1 1; jet(max(1, num_vortices))]; 
    
    % Use contourf for filled regions. The levels ensure each label gets a color.
    contourf(Cx, Cy, relabeled_mask, 0:num_vortices, 'LineStyle', 'none');
    
    hold on;
    % Mark the LSC center with a prominent star
    plot(lsc_center(2), lsc_center(1), 'p', ...
        'MarkerEdgeColor','k', 'MarkerFaceColor','y', 'MarkerSize',15, 'LineWidth', 1.5);
    hold off;
    
    axis equal tight;
    colormap(custom_colormap);
    c = colorbar('Ticks', (1:num_vortices) - 0.5, 'TickLabels', 1:num_vortices);
    c.Label.String = 'Vortex Label (Sorted by Size)';
    title({'Sorted and Labeled Vortices', ['Snapshot: ', namebase, num2str(fileNum)], ...
           ['LSC Center: (', sprintf('%.3f', lsc_center(1)), ', ', sprintf('%.3f', lsc_center(2)), ')']});
    xlabel('x'); ylabel('y');
else
    % If no vortices, just show an empty box to indicate this
    axis([0 1 0 1]);
    title({'No Vortices Found', ['Snapshot: ', namebase, num2str(fileNum)]});
    xlabel('x'); ylabel('y');
    set(gca, 'XTick', [], 'YTick', []);
end

% Save the figure as a PNG file
output_png_filename = 'vortex_visualization.png';
saveas(fig, output_png_filename);
disp(['Visualization saved to: ', output_png_filename]);

% Close the invisible figure to free up memory
close(fig);
%% ==================== Save the original continuous fields ====================
[Cx, Cy] = ndgrid(params.xGrid(1:nx), params.yGrid(1:ny)); % Match grid size to data

disp('Saving original continuous fields...');
tec_file = liton_ordered_tec.TEC_FILE;
tec_file.FileName = 'vor_z_2110';
tec_file.Variables = {'x','y','vor_z'};
tec_file.Zones = liton_ordered_tec.TEC_ZONE;
tec_file.Zones.Data = {Cx,Cy,vor_z};
tec_file.write_plt();

tec_file = liton_ordered_tec.TEC_FILE;
tec_file.FileName = 'vor_z_inden_2110';
tec_file.Variables = {'x','y','vor_z_inden'};
tec_file.Zones = liton_ordered_tec.TEC_ZONE;
tec_file.Zones.Data = {Cx,Cy,double(binary_mask)};
tec_file.write_plt();
% 
tec_file = liton_ordered_tec.TEC_FILE;
tec_file.FileName = 'vortex_sorted_labels_2110';
tec_file.Variables = {'x','y','label'};
tec_file.Zones = liton_ordered_tec.TEC_ZONE;
tec_file.Zones.Data = {Cx,Cy,relabeled_mask};
tec_file.write_plt();

% --- Save LSC center coordinates to a text file ---
if ~isnan(lsc_center(1))
    fid = fopen('lsc_center_coords.txt', 'w');
    fprintf(fid, '%.8f, %.8f\n', lsc_center(1), lsc_center(2));
    fclose(fid);
    disp('LSC center coordinates saved to lsc_center_coords.txt');
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