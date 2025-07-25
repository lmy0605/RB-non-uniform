clear; close all; clc;

%% basic settings
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

%%
Uavg = zeros(nx, ny);
Vavg = zeros(nx, ny);
Tavg = zeros(nx, ny);

for fileNum = fileNumStart:fileNumInterval:fileNumEnd
    t=fileNum-fileNumStart+1;
    if(mod(fileNum,100)==0)
        disp(['Current data fsile is ', [namebase, num2str(fileNum),'.bin']]);
    end
    [U,V,T,~] = readBinaryFile(fullfile(inputDir, [namebase, num2str(fileNum),'.bin']),nx,ny);
    U = reshape(U,nx,ny);
    V = reshape(V,nx,ny);
    T = reshape(T,nx,ny);

    U=U/params.velocityUnit;
    V=V/params.velocityUnit;
    T=(T-params.Tref)/(params.Thot-params.Tcold);

    Uavg = Uavg+U;
    Vavg = Vavg+V;
    Tavg = Tavg+T;
end

Uavg = Uavg/(fileSum);
Vavg = Vavg/(fileSum);
Tavg = Tavg/(fileSum);

%Reynold
UU_PRIME = zeros(nx, ny);
UV_PRIME = zeros(nx, ny);
VV_PRIME = zeros(nx, ny);
VT_PRIME = zeros(nx, ny);

GRAD_UX_SUM= zeros(nx, ny);
GRAD_UY_SUM= zeros(nx, ny);
GRAD_VX_SUM= zeros(nx, ny);
GRAD_VY_SUM= zeros(nx, ny);

PSEUDO11=zeros(nx, ny);
PSEUDO12=zeros(nx, ny);
PSEUDO21=zeros(nx, ny);
PSEUDO22=zeros(nx, ny);

true11=zeros(nx, ny);
true12=zeros(nx, ny);
true22=zeros(nx, ny);

diff11= zeros(nx, ny);
diff12= zeros(nx, ny);
diff22= zeros(nx, ny);

diff211= zeros(nx, ny);
diff212= zeros(nx, ny);
diff222= zeros(nx, ny);

for fileNum = fileNumStart:fileNumInterval:fileNumEnd
    t=fileNum-fileNumStart+1;
    if(mod(fileNum,100)==0)
        disp(['Current data file is ', [namebase, num2str(fileNum),'.bin']]);
    end
    [U,V,T,~] = readBinaryFile(fullfile(inputDir, [namebase, num2str(fileNum),'.bin']),nx,ny);
    U = reshape(U,nx,ny);
    V = reshape(V,nx,ny);
    T = reshape(T,nx,ny);

    U=U/params.velocityUnit;
    V=V/params.velocityUnit;
    T=(T-params.Tref)/(params.Thot-params.Tcold);

    UU_PRIME = UU_PRIME+(U-Uavg).^2;
    UV_PRIME = UV_PRIME+(V-Vavg).*(U-Uavg);
    VV_PRIME = VV_PRIME+(V-Vavg).^2;
    VT_PRIME = VT_PRIME+(V-Vavg).*(T-Tavg);

    [UX,UY,VX,VY]=GRAD1(U,V,params.dx,params.dy);
            
    GRAD_UX_SUM= GRAD_UX_SUM+UX;
    GRAD_UY_SUM= GRAD_UY_SUM+UY;
    GRAD_VX_SUM= GRAD_VX_SUM+VX;
    GRAD_VY_SUM= GRAD_VY_SUM+VY;

    [UX_prime,UY_prime,VX_prime,VY_prime]=GRAD1(U-Uavg,V-Vavg,params.dx,params.dy);

    PSEUDO11=PSEUDO11+UX_prime.^2;
    PSEUDO12=PSEUDO12+UY_prime.^2;
    PSEUDO21=PSEUDO21+VX_prime.^2;
    PSEUDO22=PSEUDO22+VY_prime.^2;

    true11=true11+UX_prime.^2;
    true12=true12+0.5*(VX_prime+UY_prime).^2;
    true22=true22+VY_prime.^2;

end

UU_PRIME = UU_PRIME/fileSum;
UV_PRIME = UV_PRIME/fileSum;
VV_PRIME = VV_PRIME/fileSum;
VT_PRIME = VT_PRIME/fileSum;

% GRAD_UX_SUM= GRAD_UX_SUM/fileSum*params.length_LB;
% GRAD_UY_SUM= GRAD_UY_SUM/fileSum*params.length_LB;
% GRAD_VX_SUM= GRAD_VX_SUM/fileSum*params.length_LB;
% GRAD_VY_SUM= GRAD_VY_SUM/fileSum*params.length_LB;
% 
% PSEUDO11=PSEUDO11/fileSum*params.length_LB.^2;
% PSEUDO12=PSEUDO12/fileSum*params.length_LB.^2;
% PSEUDO21=PSEUDO21/fileSum*params.length_LB.^2;
% PSEUDO22=PSEUDO22/fileSum*params.length_LB.^2;
% 
% true11=true11/fileSum*params.length_LB.^2;
% true12=true12/fileSum*params.length_LB.^2;
% true22=true22/fileSum*params.length_LB.^2;
% 
% diff11=diff11/fileSum*params.length_LB.^2;
% diff22=diff22/fileSum*params.length_LB.^2;
% diff12=diff12/fileSum*params.length_LB.^2;
% 
% diff211=diff211/fileSum*params.length_LB.^2;
% diff222=diff222/fileSum*params.length_LB.^2;
% diff212=diff212/fileSum*params.length_LB.^2;

GRAD_UX_SUM= GRAD_UX_SUM/fileSum;
GRAD_UY_SUM= GRAD_UY_SUM/fileSum;
GRAD_VX_SUM= GRAD_VX_SUM/fileSum;
GRAD_VY_SUM= GRAD_VY_SUM/fileSum;

PSEUDO11=PSEUDO11/fileSum;
PSEUDO12=PSEUDO12/fileSum;
PSEUDO21=PSEUDO21/fileSum;
PSEUDO22=PSEUDO22/fileSum;

true11=true11/fileSum;
true12=true12/fileSum;
true22=true22/fileSum;

shearP=-UU_PRIME.*GRAD_UX_SUM-UV_PRIME.*GRAD_UY_SUM-UV_PRIME.*GRAD_VX_SUM-VV_PRIME.*GRAD_VY_SUM;
buoyancyP=VT_PRIME;
Pruduction=shearP+buoyancyP;

pseudo_dissipation=sqrt(Prandtl/Rayleigh)*(PSEUDO11+PSEUDO22+PSEUDO21+PSEUDO12);
true_dissipation=sqrt(Prandtl/Rayleigh)*(true22+true12+true11)*2;

[~,~,pavg]=nonUniformAverage(Pruduction,params.xGrid,params.yGrid);
pavg=sum(pavg(:))/params.length0.^2;
disp(pavg);
%%
[Cx, Cy] = ndgrid(params.xGrid(1:end-1), params.yGrid(1:end-1));

%production
tec_file = liton_ordered_tec.TEC_FILE;
tec_file.FileName = 'production';
tec_file.Variables = {'X','Y','s','b','p'};
tec_file.Zones = liton_ordered_tec.TEC_ZONE;
tec_file.Zones.Data = {Cx,Cy,shearP,buoyancyP,Pruduction};
tec_file = tec_file.write_plt();

%dissipation
tec_file = liton_ordered_tec.TEC_FILE;
tec_file.FileName = 'dissipation';
tec_file.Variables = {'X','Y','pseudo','true'};
tec_file.Zones = liton_ordered_tec.TEC_ZONE;
tec_file.Zones.Data = {Cx,Cy,pseudo_dissipation,true_dissipation};
tec_file = tec_file.write_plt();

%---------------------------------------------
% read explanation here: https://zhuanlan.zhihu.com/p/62792201
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