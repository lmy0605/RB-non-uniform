clear; close all; clc;

% --- Configuration ---
Ra_vals     = [1e8,  1e9,   1e10,   1e11,   1e12];
nx_vals     = [257,  513,   1025,   2049,   4097];
constA_vals = [1.5,  2.1,   2.5,    2.9,    3.1];
Prandtl     = 0.71;
baseDir     = 'C:\Users\user\OneDrive\Desktop\01-RB-0612\fig14-tke'; % 基础目录
p           = zeros(1, length(Ra_vals));
d           = zeros(1, length(Ra_vals));

% --- Main Loop ---
for i = 1:length(Ra_vals)
    Rayleigh = Ra_vals(i);
    nx = nx_vals(i);
    ny = nx;
    constA = constA_vals(i);
    
    exponent = log10(Rayleigh);
    folderName = sprintf('1e%d', exponent); 
    inputDir = fullfile(baseDir, folderName, 'production.plt');
    
    fprintf('Processing Ra = %g, nx = %d...\n', Rayleigh, nx);
    
    params = calculateSystemParameters(nx, ny, Rayleigh, Prandtl, constA, 'log.log');
    
    if ~exist(inputDir, 'file')
        warning('File not found, skipping: %s', inputDir);
        p(i) = NaN; 
        continue;
    end
    [production, ~] = read_tecplot_plt(inputDir);

    production_res = reshape(production{5}, nx, ny);
    
    [~, ~, production_w] = nonUniformAverage(production_res, params.xGrid, params.yGrid);

    p(i) = sum(production_w(:)) / params.length0^2;
    
    %---------------------------diaasipation-------------------------------------------------
    inputDir2 = fullfile(baseDir, folderName, 'dissipation.plt');
   
    if ~exist(inputDir, 'file')
        warning('File not found, skipping: %s', inputDir2);
        p(i) = NaN; 
        continue;
    end
    [dissipation, ~] = read_tecplot_plt(inputDir2);

    dissipation_res = reshape(dissipation{4}, nx, ny);
    
    [~, ~, dissipation_w] = nonUniformAverage(dissipation_res, params.xGrid, params.yGrid);

    d(i) = sum(dissipation_w(:)) / params.length0^2;
end

% --- Plotting ---
disp('Plotting results...');
plot_matlab('x', Ra_vals, 'y', p, ...
            'CanvasRatio', [10, 7], ...
            'LineWidth', 1.5, ...
            'LineStyle',  '-', 'LineColor','color1',...
            'Marker', 'o', ...
            'MarkerInterval', 1, ...
            'MarkerSize', 8, ...
            'MarkerFaceColor', 'w', ...
            'XLabel', '\it{Ra}', ...
            'YLabel', '\it{Production}', 'YLim', [1e-4 2e-3], ...
            'LogScale', 'xy', ...
            'OutputFilename', 'p_ra.png', ...
            'Resolution', 600);
plot_matlab('x', Ra_vals, 'y', d, ...
            'CanvasRatio', [10, 7], ...
            'LineWidth', 1.5, ...
            'LineStyle',  '-', 'LineColor','color2',...
            'Marker', 's', ...
            'MarkerInterval',1,...
            'MarkerSize', 8, ...
            'MarkerFaceColor', 'w', ...
            'XLabel', '\it{Ra}', ...
            'YLabel', '\it{Dissipation}', ...
            'LogScale','xy','YLim',[1e-4 2e-3],...
            'OutputFilename', 'd_ra.png', ...
            'Resolution', 600);
%%
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