clear; close all; clc;

%% basic settings
fileNumStart=1001;
fileNumEnd=10000;
fileNumInterval=1;
fileSum=fileNumEnd-fileNumStart+1;
inputDir = '/nfsdata4/AXu/RB-non-uniform/Ra1e8-mesh257/binFile-1-10000/'; % please rename data folder as "binFile"
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
[Cx, Cy] = ndgrid(params.xGrid(1:end-1), params.yGrid(1:end-1));
modeNum = 0;
for m = 1:1:rankNum
    for n = 1:1:rankNum
        modeNum = modeNum+1;
        uBasis = 2.0*sin(m*pi*Cx).*cos(n*pi*Cy);
        vBasis = -2.0*cos(m*pi*Cx).*sin(n*pi*Cy);
        if (outputModeIllustration == 1)
            exportFourier(nx, ny,uBasis,vBasis,m,n,outputDir);
        end
        uHat(:,modeNum) = reshape(uBasis,[],1);
        vHat(:,modeNum) = reshape(vBasis,[],1);
    end
end
%%
disp(["Begin to get Fourier coefficient >>>>>>>>>"])
for fileNum = fileNumStart:fileNumInterval:fileNumEnd
    t=fileNum-fileNumStart+1;
    if(mod(fileNum,100)==0)
        disp(['Current data file is ', [namebase, num2str(fileNum),'.bin']]);
    end
    [U,V,~,~] = readBinaryFile(fullfile(inputDir, [namebase, num2str(fileNum),'.bin']),nx,ny);
%     U = reshape(U,nx,ny);
%     V = reshape(V,nx,ny);

    U=U/params.velocityUnit;
    V=V/params.velocityUnit;

    Ax(:,fileNum-fileNumStart+1) = U'*uHat;
    Ay(:,fileNum-fileNumStart+1) = V'*vHat;
end

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
    result(modeNum, 5) = result(modeNum,3)/sum(result(:,3))*100.0; 
    result(modeNum, 6) = result(modeNum,3)/energyRMS(modeNum); %S11=<E^(m,n)>/E^{m,n}_rms
    result(modeNum, 7) = energyRMS(modeNum)/mean(energyAll); %E^{m,n}_rms/<E_total>
end
save(['Fourier_result_summary.dat'], 'result', '-ascii')

%%
if exportFigure==1
    energyPercentage_cell = num2cell(energyPercentage, 2);

    plot_matlab('x', linspace(0,fileSum,fileSum)', 'y', energyPercentage_cell, ...
            'CanvasRatio', [10, 7], ...
            'LineWidth', 1.5, ...
            'XLabel', '\it{t} ', ...
            'YLabel', '\it{E^{m,n}} \rm{(}\it{t}\rm{)} \it{/E_{total}} \rm{(}\it{t}\rm{)}', ...
            'ShowLegend',0,...
            'OutputFilename', 'energyPercentage.png', ...
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
X = linspace(1,gridNumX,gridNumX);
Y = linspace(1,gridNumY,gridNumY);
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