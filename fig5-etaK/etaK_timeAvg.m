clear; close all; clc;

%% basic settings
fileNumStart=2001;
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
etaKAvg=(viscosity.^3./dissipation).^0.25;

%%
[deltax,deltay,deltaxy]=calculate_node_area_weights(params.xGrid,params.yGrid);
[deltaX,deltaY]=ndgrid(deltax,deltay);

grid_resolution = max(deltaX, deltaY); 
grid_resolution2 = sqrt(deltaxy);

resolution_ratio_etaK = grid_resolution ./ etaKAvg;
resolution_ratio_etaU = grid_resolution ./ etaUAvg;
resolution_ratio_etaK2 = grid_resolution2 ./ etaKAvg;
resolution_ratio_etaU2 = grid_resolution2 ./ etaUAvg;
aspect_ratio=max(deltaX, deltaY)./min(deltaX, deltaY);
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

tec_file = liton_ordered_tec.TEC_FILE;
tec_file.FileName = strcat('resolution_ratio_etaK_',casename);
tec_file.Variables = {'x','y','resolution_ratio_etaK','aspect_ratio'};
tec_file.Zones = liton_ordered_tec.TEC_ZONE;
tec_file.Zones.Data = {Cx,Cy,resolution_ratio_etaK,aspect_ratio};
tec_file = tec_file.write_plt();

tec_file = liton_ordered_tec.TEC_FILE;
tec_file.FileName = strcat('resolution_ratio_etaU_',casename);
tec_file.Variables = {'x','y','resolution_ratio_etaU','aspect_ratio'};
tec_file.Zones = liton_ordered_tec.TEC_ZONE;
tec_file.Zones.Data = {Cx,Cy,resolution_ratio_etaU,aspect_ratio};
tec_file = tec_file.write_plt();

tec_file = liton_ordered_tec.TEC_FILE;
tec_file.FileName = strcat('2resolution_ratio_etaK_',casename);
tec_file.Variables = {'x','y','resolution_ratio_etaK','aspect_ratio'};
tec_file.Zones = liton_ordered_tec.TEC_ZONE;
tec_file.Zones.Data = {Cx,Cy,resolution_ratio_etaK2,aspect_ratio};
tec_file = tec_file.write_plt();

tec_file = liton_ordered_tec.TEC_FILE;
tec_file.FileName = strcat('2resolution_ratio_etaU_',casename);
tec_file.Variables = {'x','y','resolution_ratio_etaU','aspect_ratio'};
tec_file.Zones = liton_ordered_tec.TEC_ZONE;
tec_file.Zones.Data = {Cx,Cy,resolution_ratio_etaU2,aspect_ratio};
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
function [dx_contrib,dy_contrib,node_area_weights] = calculate_node_area_weights(x_node_coords, y_node_coords)
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

    node_area_weights = dx_contrib * dy_contrib'; % Results in (nx x ny) matrix
end