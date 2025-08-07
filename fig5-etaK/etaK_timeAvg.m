clear; close all; clc;

%% basic settings
fileNumStart=2086;
fileNumEnd=2609;
fileNumInterval=1;
fileSum=fileNumEnd-fileNumStart+1;
inputDir = '/nfsdata2/AXu/RB-non-uniform/Ra1e12-mesh4097/binFile-1971-2609/'; % please rename data folder as "binFile"
namebase = 'buoyancyCavity-';

% --- Check if input directory exists before proceeding ---
disp('Verifying input directory...');
if ~isfolder(inputDir)
    error('Input directory not found: %s', inputDir);
end
disp(['Input directory found: ', inputDir]);

nx=4097;
ny=nx;
constA=3.1;
Rayleigh=1e12;
Prandtl=0.71;

params = calculateSystemParameters(nx,ny, Rayleigh, Prandtl,constA,'log.log');
viscosity=sqrt(Prandtl/Rayleigh);
%% calculation
% Kolmogorov scale using  kinetic energy dissipation rates
UAvg = zeros(nx,ny);
VAvg = zeros(nx,ny);
EUavg = zeros(nx,ny);

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

    UAvg=UAvg+U;
    VAvg=VAvg+V;

    [UX,UY,VX,VY]=GRAD1(U,V,params.dx,params.dy);
    EUavg=EUavg+(2.*(UX).^2+2.*(VY).^2+(VX+UY).^2)*viscosity;
end

UAvg=UAvg/fileSum;
VAvg=VAvg/fileSum;
EUavg=EUavg/fileSum;
etaUAvg=(viscosity.^3./EUavg).^0.25;

% Kolmogorov scale using tke dissipation rates
UPrime = zeros(nx,ny);
VPrime = zeros(nx,ny);
dissipation=zeros(nx,ny);

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

    UPrime=U-UAvg;
    VPrime=V-VAvg;

    [UX_prime,UY_prime,VX_prime,VY_prime]=GRAD1(UPrime,VPrime,params.dx,params.dy);
    dissipation=dissipation+(2*UX_prime.^2+2*VY_prime.^2+(VX_prime+UY_prime).^2)*viscosity;
end
dissipation=dissipation/fileSum;
etaKAvg=(viscosity.^3./dissipationAvg).^0.25;
%%
[Cx, Cy] = ndgrid(params.xGrid(1:end-1), params.yGrid(1:end-1));

%
tec_file = liton_ordered_tec.TEC_FILE;
tec_file.FileName = 'etaU_timeAvg';
tec_file.Variables = {'x','y','eta'};
tec_file.Zones = liton_ordered_tec.TEC_ZONE;
tec_file.Zones.Data = {Cx,Cy,etaUAvg};
tec_file = tec_file.write_plt();

tec_file = liton_ordered_tec.TEC_FILE;
tec_file.FileName = 'etaK_timeAvg';
tec_file.Variables = {'x','y','eta'};
tec_file.Zones = liton_ordered_tec.TEC_ZONE;
tec_file.Zones.Data = {Cx,Cy,etaKAvg};
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

function meanU = cumulativeAverage(U)
    % CUMULATIVEAVERAGE 
    if ~isvector(U), error('input is not an array'); end
    U = U(:).'; 
    cumulative_sum = cumsum(U);
    divisors = 1:length(U);
    meanU = cumulative_sum ./ divisors;
end

function varPopU = cumulativePopulationVariance(U)
    % CUMULATIVEPOPULATIONVARIANCE 
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

function [GRAD_UX,GRAD_UY,GRAD_VX,GRAD_VY]=GRAD1(U,V,dx,dy)
%center df=(dx(i)^2*(f(i+1)-f(i))+dx(i+1)^2*(f(i)-f(i-1)))/(dx(i)*dx(i+1)*(dx(i)+dx(i+1)))
%forward df=(-(2*dx(i+1)*dx(i+2)+dx(i+2)^2)*f(i)+(dx(i+1)+dx(i+2))^2*f(i+1)-dx(i+1)^2*f(i+2))/(dx(i+2)*dx(i+1)*(dx(i+2)+dx(i+1)))
%backward df=((2*dx(i)*dx(i-1)+dx(i-1)^2)*f(i)-(dx(i)+dx(i-1))^2*f(i-1)+dx(i)^2*f(i-2))/(dx(i)*dx(i-1)*(dx(i)+dx(i-1)))
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