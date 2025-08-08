clear; close all; clc;

%% Basic settings
% --- Define first data set ---
fileNumStart1=1;
fileNumEnd1=2000;
inputDir1 = '/nfsdata4/AXu/RB-non-uniform/Ra1e11-mesh2049/binFile-1-2000/';
timeStep1 = 0.5; % Time interval per file (in free-fall time)

% --- Define second data set ---
fileNumStart2=2001;
fileNumEnd2=5839;
inputDir2 = '/nfsdata4/AXu/RB-non-uniform/Ra1e11-mesh2049/binFile-2001-5839/';
timeStep2 = 1.0; % Time interval per file (in free-fall time)
timeStart2 = fileNumEnd1 * timeStep1; % Physical start time for the second data set

fileNumStart3 = 5840;
fileNumEnd3 =11001;
inputDir3 = '/nfsdata4/AXu/RB-non-uniform/Ra1e11-mesh2049/binFile-5840-11001/';
timeStep3 = 1.0; % Time interval per file (in free-fall time)
% Calculate start time for the third set based on the end of the second set
timeStart3 = timeStart2 + (fileNumEnd2 - fileNumStart2 + 1) * timeStep2;

% --- Calculate total file count ---
fileSum1 = fileNumEnd1 - fileNumStart1 + 1;
fileSum2 = fileNumEnd2 - fileNumStart2 + 1;
fileSum3 = fileNumEnd3 - fileNumStart3 + 1;
totalFileCount = fileSum1 + fileSum2 + fileSum3;

namebase = 'buoyancyCavity-';

% --- Check if input directory exists before proceeding ---
disp('Verifying input directory...');
if ~isfolder(inputDir3)
    error('Input directory not found: %s', inputDir3);
end
disp(['Input directory found: ', inputDir3]);

% --- Physical and grid parameters ---
nx=2049;
ny=nx;
constA=2.9;
Rayleigh=1e11; % Update Rayleigh number
Prandtl=0.71;

params = calculateSystemParameters(nx,ny, Rayleigh, Prandtl,constA,'log.log');

%% Calculation
% --- Pre-allocate memory ---
ReAvg = zeros(1,totalFileCount);
NuWallAvg = zeros(1,totalFileCount);
physicalTime = zeros(1,totalFileCount); % To store the actual physical time

deri_Ty_bottom = zeros(nx, totalFileCount);
deri_Ty_top = zeros(nx, totalFileCount);

t_idx = 0; % Global index for the time series

% --- Loop to process the first set of files ---
disp(['Processing first set of files from: ', inputDir1]);
for fileNum = fileNumStart1:fileNumEnd1
    t_idx = t_idx + 1;
    if(mod(fileNum,100)==0)
        disp(['Current data file is ', [namebase, num2str(fileNum),'.bin']]);
    end
    
    % Calculate physical time
    physicalTime(t_idx) = fileNum * timeStep1;

    % Read and process data
    [U,V,T,~] = readBinaryFile(fullfile(inputDir1, [namebase, num2str(fileNum),'.bin']),nx,ny);
    U = reshape(U,nx,ny);
    V = reshape(V,nx,ny);
    T = reshape(T,nx,ny);

    U=U/params.velocityUnit;
    V=V/params.velocityUnit;

    [~,~,Re] =nonUniformAverage(V.^2+U.^2,params.xGrid,params.yGrid);
    [deri_Ty_bottom(:,t_idx),deri_Ty_top(:,t_idx)] =deri_Twall(T,params.Thot,params.Tcold,params.length_LB,params.gridRatioQ,nx,ny);
    [deri_Ty_bottom_xWeighted,~,~]=nonUniformAverage(deri_Ty_bottom(:,t_idx),params.xGrid,params.yGrid);
    [deri_Ty_top_xWeighted,~,~]=nonUniformAverage(deri_Ty_top(:,t_idx),params.xGrid,params.yGrid);

    ReAvg(1,t_idx)=sqrt(Rayleigh/Prandtl)*sqrt(sum(Re(:))/params.length0.^2);
    NuWallAvg(1,t_idx)=-0.5*(sum(deri_Ty_bottom_xWeighted(:))+sum(deri_Ty_top_xWeighted(:)))/params.length0;
end

% --- Loop to process the second set of files ---
disp(['Processing second set of files from: ', inputDir2]);
for fileNum = fileNumStart2:fileNumEnd2
    t_idx = t_idx + 1;
    if(mod(fileNum,100)==0)
        disp(['Current data file is ', [namebase, num2str(fileNum),'.bin']]);
    end
    
    % Calculate physical time
    physicalTime(t_idx) = timeStart2 + (fileNum - fileNumStart2 + 1) * timeStep2;

    % Read and process data
    [U,V,T,~] = readBinaryFile(fullfile(inputDir2, [namebase, num2str(fileNum),'.bin']),nx,ny);
    U = reshape(U,nx,ny);
    V = reshape(V,nx,ny);
    T = reshape(T,nx,ny);

    U=U/params.velocityUnit;
    V=V/params.velocityUnit;

    [~,~,Re] =nonUniformAverage(V.^2+U.^2,params.xGrid,params.yGrid);
    [deri_Ty_bottom(:,t_idx),deri_Ty_top(:,t_idx)] =deri_Twall(T,params.Thot,params.Tcold,params.length_LB,params.gridRatioQ,nx,ny);
    [deri_Ty_bottom_xWeighted,~,~]=nonUniformAverage(deri_Ty_bottom(:,t_idx),params.xGrid,params.yGrid);
    [deri_Ty_top_xWeighted,~,~]=nonUniformAverage(deri_Ty_top(:,t_idx),params.xGrid,params.yGrid);

    ReAvg(1,t_idx)=sqrt(Rayleigh/Prandtl)*sqrt(sum(Re(:))/params.length0.^2);
    NuWallAvg(1,t_idx)=-0.5*(sum(deri_Ty_bottom_xWeighted(:))+sum(deri_Ty_top_xWeighted(:)))/params.length0;
end

disp(['Processing third set of files from: ', inputDir3]);
for fileNum = fileNumStart3:fileNumEnd3
    t_idx = t_idx + 1;
    if(mod(fileNum,100)==0)
        disp(['Current data file is ', [namebase, num2str(fileNum),'.bin']]);
    end
    
    % Calculate physical time
    physicalTime(t_idx) = timeStart3 + (fileNum - fileNumStart3 + 1) * timeStep3;

    % Read and process data
    [U,V,T,~] = readBinaryFile(fullfile(inputDir3, [namebase, num2str(fileNum),'.bin']),nx,ny);
    U = reshape(U,nx,ny);
    V = reshape(V,nx,ny);
    T = reshape(T,nx,ny);

    U=U/params.velocityUnit;
    V=V/params.velocityUnit;

    [~,~,Re] =nonUniformAverage(V.^2+U.^2,params.xGrid,params.yGrid);
    [deri_Ty_bottom(:,t_idx),deri_Ty_top(:,t_idx)] =deri_Twall(T,params.Thot,params.Tcold,params.length_LB,params.gridRatioQ,nx,ny);
    [deri_Ty_bottom_xWeighted,~,~]=nonUniformAverage(deri_Ty_bottom(:,t_idx),params.xGrid,params.yGrid);
    [deri_Ty_top_xWeighted,~,~]=nonUniformAverage(deri_Ty_top(:,t_idx),params.xGrid,params.yGrid);

    ReAvg(1,t_idx)=sqrt(Rayleigh/Prandtl)*sqrt(sum(Re(:))/params.length0.^2);
    NuWallAvg(1,t_idx)=-0.5*(sum(deri_Ty_bottom_xWeighted(:))+sum(deri_Ty_top_xWeighted(:)))/params.length0;
end
%% Post-processing and Output
disp(['Total files processed: ', num2str(t_idx)]);
disp(['ReAvg = ', num2str(mean(ReAvg(:)))])
disp(['NuWallAvg = ', num2str(mean(NuWallAvg(:)))])

meanReAvg_cumsum = cumulativeAverage(ReAvg);
meanNuWallAvg_cumsum = cumulativeAverage(NuWallAvg);

varReAvg_cumsum = cumulativePopulationVariance(ReAvg);
varNuWallAvg_cumsum = cumulativePopulationVariance(NuWallAvg);

% --- Output instantaneous time series ---
tec_file = liton_ordered_tec.TEC_FILE;
tec_file.FileName = 'timeseries_ReNu';
tec_file.Variables = {'PhysicalTime','Re','Nu'}; % Use 'PhysicalTime' as the variable name
tec_file.Zones = liton_ordered_tec.TEC_ZONE;
% Use the physicalTime array as the first variable (time axis)
tec_file.Zones.Data = {physicalTime,ReAvg,NuWallAvg}; 
tec_file = tec_file.write_plt();

% --- Output cumulative mean and variance time series ---
tec_file = liton_ordered_tec.TEC_FILE;
tec_file.FileName = 'timeseries_ReNu_mean';
tec_file.Variables = {'PhysicalTime','Re_mean','Nu_mean','Re_var','Nu_var'};
tec_file.Zones = liton_ordered_tec.TEC_ZONE;
% Use the physicalTime array as the first variable (time axis)
tec_file.Zones.Data = {physicalTime,meanReAvg_cumsum,meanNuWallAvg_cumsum,varReAvg_cumsum,varNuWallAvg_cumsum};
tec_file = tec_file.write_plt();

disp('All processing and file writing complete.');

%% Functions (unchanged from original)
function [U, V ,T, rho] = readBinaryFile(file, nx, ny)
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

function [deri_Ty_bottom,deri_Ty_top]=deri_Twall(T,T_bottom,T_top,length_LB,q,nx,ny)
deri_Ty_bottom=zeros(nx,1);
deri_Ty_top=zeros(nx,1);
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