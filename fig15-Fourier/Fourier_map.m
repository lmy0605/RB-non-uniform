clear; close all; clc;

%% basic settings
fileNumStart=2001;
fileNumEnd=10000;
fileNumInterval=1;
fileSum=fileNumEnd-fileNumStart+1;
inputDir = '/nfsdata4/AXu/RB-non-uniform/Ra1e8-mesh257/binFile-1-10000/'; 
outputDir=strcat(pwd,'/');
namebase = 'buoyancyCavity-';

% --- Check if input directory exists before proceeding ---
disp('Verifying input directory...');
if ~isfolder(inputDir)
    error('Input directory not found: %s', inputDir);
end
disp(['Input directory found: ', inputDir]);

nx=257;
ny=nx;
constA=1.5;
Rayleigh=1e8;
Prandtl=0.71;

params = calculateSystemParameters(nx,ny, Rayleigh, Prandtl,constA,'log.log');
viscosity=sqrt(Prandtl/Rayleigh);

rankNum = 3;
outputModeIllustration = 0;

exportFigure=1;
%%
nx_old=257;
ny_old=nx_old;
[x_old, y_old] = generate_grid(nx_old, ny_old, 0, constA);

nx_new=257;
ny_new=nx_new;
[x_new, y_new] = generate_grid(nx_new, ny_new, 1,1);
%%
[Cx, Cy] = ndgrid(x_new, y_new);
% [Cx, Cy] = ndgrid(params.xGrid(1:end-1), params.yGrid(1:end-1));
modeNum = 0;
for m = 1:1:rankNum
    for n = 1:1:rankNum
        modeNum = modeNum+1;
        uBasis = 2.0*sin(m*pi*Cx).*cos(n*pi*Cy);
        vBasis = -2.0*cos(m*pi*Cx).*sin(n*pi*Cy);
        if (outputModeIllustration == 1)
            exportFourier(params.xGrid(1:end-1), params.yGrid(1:end-1),uBasis,vBasis,m,n,outputDir);
        end
        uHat(:,modeNum) = reshape(uBasis,[],1);
        vHat(:,modeNum) = reshape(vBasis,[],1);
    end
end
%%
disp(["Begin to get Fourier coefficient >>>>>>>>>"])
[~,~,deltaxy]=calculate_node_area_weights(params.xGrid,params.yGrid);
deltaxy=reshape(deltaxy,[],1);

NuVolAvg_old = zeros(1,fileSum);
NuVolAvg_new = zeros(1,fileSum);
ReAvg_old = zeros(1,fileSum);
ReAvg_new = zeros(1,fileSum);

for fileNum = fileNumStart:fileNumInterval:fileNumEnd
    t=fileNum-fileNumStart+1;
    if(mod(fileNum,100)==0)
        disp(['Current data file is ', [namebase, num2str(fileNum),'.bin']]);
    end
    [U,V,T,~] = readBinaryFile(fullfile(inputDir, [namebase, num2str(fileNum),'.bin']),nx,ny);

    U = reshape(U,nx_old,ny_old);
    V = reshape(V,nx_old,ny_old);
    T = reshape(T,nx_old,ny_old);
    data_old = {U, V, T};

    data_new=conservative_remap_2d(x_old, y_old, data_old, x_new, y_new);

    U_new=reshape(data_new{1},nx_new,ny_new);
    V_new=reshape(data_new{2},nx_new,ny_new);
    T_new=reshape(data_new{3},nx_new,ny_new);

    U_new=U_new/params.velocityUnit;
    V_new=V_new/params.velocityUnit;
    U_old=U/params.velocityUnit;
    V_old=V/params.velocityUnit;
    T_old=T;

    [~,~,NuVol_old] =nonUniformAverage(V_old.*T_old,params.xGrid,params.yGrid);
    NuVolAvg_old(1,t)=sqrt(Prandtl*Rayleigh)*sum(NuVol_old(:))/params.length0.^2+1;
    NuVol_new=V_new.*T_new;
    NuVolAvg_new(1,t)=sqrt(Prandtl*Rayleigh)*mean(NuVol_new(:))+1;

    [~,~,Re_old] =nonUniformAverage(V_old.^2+U_old.^2,params.xGrid,params.yGrid);
    ReAvg_old(1,t)=sqrt(Rayleigh/Prandtl)*sqrt(sum(Re_old(:))/params.length0.^2);
    Re_new=V_new.^2+U_new.^2;
    ReAvg_new(1,t)=sqrt(Rayleigh/Prandtl)*sqrt(mean(Re_new(:)));

    U_new = reshape(U_new,[],1);
    V_new = reshape(V_new,[],1);

    for modeNum=1:rankNum^2
        Ax(:,fileNum-fileNumStart+1) = U_new'*uHat;
        Ay(:,fileNum-fileNumStart+1) = V_new'*vHat;
    end
end

disp(['NuVolAvg_old = ', num2str(mean(NuVolAvg_old(:)))]);
disp(['NuVolAvg_new = ', num2str(mean(NuVolAvg_new(:)))]);
disp(['ReAvg_old = ', num2str(mean(ReAvg_old(:)))]);
disp(['ReAvg_new = ', num2str(mean(ReAvg_new(:)))]);
%%
%Evaluate contributions from each mode
energyCoeff = sqrt(Ax.^2 + Ay.^2);
energyAll = sum(energyCoeff,1);  %E_total(t)
energyPercentage = energyCoeff./energyAll * 100;  %E^{m,n}(t)/E_total(t)
energyRMS=std(energyCoeff, 1,2);  %E_rms^{m,n}

for modeNum = 1:rankNum^2
    output(:,1) = linspace(1, size(Ax, 2), size(Ax, 2));
    output(:,2) = Ax(modeNum, :);
    save([outputDir,'Fourier_coeff_Ax_',num2str(modeNum),'.dat'], 'output', '-ascii')
    
    output(:,1) = linspace(1, size(Ay, 2), size(Ay, 2));
    output(:,2) = Ay(modeNum, :);
    save([outputDir,'Fourier_coeff_Ay_',num2str(modeNum),'.dat'], 'output', '-ascii')
    
    output(:,1) = linspace(1, size(energyPercentage, 2), size(energyPercentage, 2));
    output(:,2) = energyPercentage(modeNum, :);
    save([outputDir,'energyPercentage_',num2str(modeNum),'.dat'], 'output', '-ascii')
end

for modeNum = 1:rankNum^2
    result(modeNum, 1) = mean(abs(Ax(modeNum, :)));
    result(modeNum, 2) = mean(abs(Ay(modeNum, :)));
    result(modeNum, 3) = mean(energyCoeff(modeNum, :));
    result(modeNum, 4) = mean(energyPercentage(modeNum, :));  %<E^{m,n}(t)/E_total(t)>
end
for modeNum = 1:rankNum^2
    result(modeNum, 5) = result(modeNum,3)/mean(energyAll)*100.0;  %<E^{m,n}(t)>/<E_total(t)>
    result(modeNum, 6) = result(modeNum,3)/energyRMS(modeNum); %S11=<E^(m,n)>/E^{m,n}_rms
    result(modeNum, 7) = energyRMS(modeNum)/mean(energyAll); %E^{m,n}_rms/<E_total>
end
save(['Fourier_result_summary.dat'], 'result', '-ascii')

save(['NuRe_compare.dat'],'NuVolAvg_old','NuVolAvg_new',"ReAvg_old","ReAvg_new",'-ascii');
%%
if exportFigure==1
    energyPercentage_cell = num2cell(energyPercentage, 2);

    plot_matlab('x', linspace(0,fileSum,fileSum)', 'y', energyPercentage_cell, ...
           'CanvasRatio', [10, 7], ...
            'LineWidth', 1.2, ...
            'LineColor', {'color1','color9','color8','color6','color7','color5','color3','color4','color2'}, ...
            'YLim', [0, 80], ...
            'XLim', [0, fileSum], ... 
            'XLabel', '\it{t}', ...
            'YLabel', '\it{E^{m,n}} \rm{(}\it{t}\rm{)}\it{/E_{total}} \rm{(}\it{t}\rm{)} (%)', ...
            'ShowLegend', 0, ...
            'OutputFilename', sprintf('energyPercentage_1e%d.png', log10(Rayleigh)), ...
            'Resolution', 600);

    plot_matlab('x', linspace(0,fileSum,fileSum)', 'y', {NuVolAvg_old,NuVolAvg_new}, ...
           'CanvasRatio', [10, 7], ...
            'LineWidth', 1.2, ...
            'XLim', [0, fileSum], ... 
            'XLabel', '\it{t}', ...
            'YLabel', '\it{Nu_{vol}}', ...
            'ShowLegend', 0, ...
            'OutputFilename', sprintf('Nu_vol_1e%d.png', log10(Rayleigh)), ...
            'Resolution', 600);

      plot_matlab('x', linspace(0,fileSum,fileSum)', 'y', {ReAvg_old,ReAvg_new}, ...
           'CanvasRatio', [10, 7], ...
            'LineWidth', 1.2, ...
            'XLim', [0, fileSum], ... 
            'XLabel', '\it{t}', ...
            'YLabel', '\it{Re}', ...
            'ShowLegend', 0, ...
            'OutputFilename', sprintf('Re_1e%d.png', log10(Rayleigh)), ...
            'Resolution', 600);
end
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

function exportFourier(gridNumX, gridNumY, U, V, M, N, outputDir)
X = gridNumX;
Y = gridNumY;
fileID=fopen([outputDir,'FourierBasis_', num2str(M),num2str(N),'.dat'],'w');
fprintf(fileID, 'Title="RB-flow"\n');
fprintf(fileID, 'Variables=X, Y, U, V \n');
fprintf(fileID, 'zone i= %d j = %d f = point \n', size(X,2), size(Y,2));
for j=1:1:size(Y,2)
    for i=1:1:size(X,2)
        fprintf(fileID, '%13.6f %13.6f %13.6f %13.6f\n', X(i), Y(j), U(i,j), V(i,j));
    end
end
fclose(fileID);
end

function [dx_contrib,dy_contrib,node_area_weights] = calculate_node_area_weights(x_node_coords, y_node_coords)
    % x_node_coords: 1D array of x-coordinates of nodes (length nx)
    % y_node_coords: 1D array of y-coordinates of nodes (length ny)

    nx_nodes = length(x_node_coords)-1;
    ny_nodes = length(y_node_coords)-1;

    dx_contrib = zeros(nx_nodes, 1); 
    if nx_nodes == 1
        warning('Only one X node.');
        dx_contrib(1) = 1; 
    else
        dx_contrib(1) = (x_node_coords(2) - 0) / 2;
        for i = 2:(nx_nodes)
            dx_contrib(i) = (x_node_coords(i+1) - x_node_coords(i-1)) / 2;
        end
    end

    dy_contrib = zeros(ny_nodes, 1);
    if ny_nodes == 1
        warning('Only one Y node.');
        dy_contrib(1) = 1; 
    else
        dy_contrib(1) = (y_node_coords(2) - 0) / 2;
        for i = 2:(ny_nodes)
            dy_contrib(i) = (y_node_coords(i+1) - y_node_coords(i-1)) / 2;
        end
    end
    
    if any(dx_contrib < 0)
        error('Negative dx_contrib.');
    end
    if any(dy_contrib < 0)
        error('Negative dy_contrib.');
    end

    node_area_weights = dx_contrib * dy_contrib'; 
end
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

function [x, y] = generate_grid(nx, ny, is_uniform, stretch_factor)
    % --- Generate X-coordinates ---
    if is_uniform
        % Generate a uniform grid from 0 to 1, with points at cell centers
        x = ( (1:nx) - 0.5 ) / nx;
    else
        % Generate a non-uniform grid stretched towards the center and boundaries
        % Create a normalized coordinate from 0 to 1
        xi = (1:nx) / (nx + 1);
        % Apply the error function based stretching formula
        x = 0.5 * (erf(stretch_factor * (xi - 0.5)) / erf(0.5 * stretch_factor) + 1);
    end

    % --- Generate Y-coordinates ---
    % Assumes the same logic is used for the y-direction
    if is_uniform
        y = ( (1:ny) - 0.5 ) / ny;
    else
        xi = (1:ny) / (ny + 1);
        y = 0.5 * (erf(stretch_factor * (xi - 0.5)) / erf(0.5 * stretch_factor) + 1);
    end

end