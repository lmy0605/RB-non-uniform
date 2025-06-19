clear; close all; clc;

%% basic settings
fileNumStart=1501;
fileNumEnd=10000;
fileNumInterval=1;
fileSum=fileNumEnd-fileNumStart+1;
inputDir = '/nfsdata4/AXu/RB-non-uniform/Ra1e9-mesh513/binFile-1-10000/'; % please rename data folder as "binFile"
namebase = 'buoyancyCavity-';

nx=513;
ny=nx;
constA=2.1;
Rayleigh=1e9;
Prandtl=0.71;

params = calculateSystemParameters(nx,ny, Rayleigh, Prandtl,constA,'log.log');

%% calculation
[Cx, Cy] = ndgrid(params.xGrid(1:end-1), params.yGrid(1:end-1));

%
fileNum=fileNumEnd;
[U,V,T,rho] = readBinaryFile(fullfile(inputDir, [namebase, num2str(fileNum),'.bin']),nx,ny);
U = reshape(U,nx,ny);
V = reshape(V,nx,ny);
T = reshape(T,nx,ny);
rho = reshape(rho,nx,ny);

U=U/params.velocityUnit;
V=V/params.velocityUnit;

u_module=(U.^2+V.^2).^0.5;

tec_file = liton_ordered_tec.TEC_FILE;
tec_file.FileName = 'inst_uvT';
tec_file.Variables = {'X','Y','U','V','T','Umodule'};
tec_file.Zones = liton_ordered_tec.TEC_ZONE;
tec_file.Zones.Data = {Cx,Cy,U,V,T,u_module};
tec_file = tec_file.write_plt();

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
