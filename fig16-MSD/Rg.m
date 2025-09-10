clear; close all; clc;
%%
inputDir = 'C:\Users\user\OneDrive\Desktop\01-RB-0612\fig12-vor\1e8\lsc_trajectory.plt';
nx=257;
ny=nx;
constA=1.5;
Rayleigh=1e8;
Prandtl=0.71;
params = calculateSystemParameters(nx,ny, Rayleigh, Prandtl,constA,'log.log');

[LSC_center, ~] = read_tecplot_plt(inputDir);
LSC_center_x=LSC_center{1};
LSC_center_y=LSC_center{2};

save('LSC_1e8.mat', ...
     'LSC_center_x', ...
     'LSC_center_y');

LSC_center_xmean=mean(LSC_center_x);
LSC_center_ymean=mean(LSC_center_y);
LSC_center_x=(LSC_center_x-LSC_center_xmean).^2;
LSC_center_y=(LSC_center_y-LSC_center_ymean).^2;

Rg_1e8=sqrt((sum(LSC_center_x(:))+sum(LSC_center_y(:)))/length(LSC_center_y))/params.length0;

%% --- Configuration for Supplementary Data ---
inputDir = 'C:\Users\user\OneDrive\Desktop\01-RB-0612\fig12-vor\1e9\lsc_trajectory.plt';
nx=513;
ny=nx;
constA=2.1;
Rayleigh=1e9;
Prandtl=0.71;
params = calculateSystemParameters(nx,ny, Rayleigh, Prandtl,constA,'log.log');

[LSC_center, ~] = read_tecplot_plt(inputDir);
LSC_center_x=LSC_center{1};
LSC_center_y=LSC_center{2};
LSC_center_x=LSC_center_x(1001:9000);
LSC_center_y=LSC_center_y(1001:9000);

save('LSC_1e9.mat', ...
     'LSC_center_x', ...
     'LSC_center_y');

LSC_center_xmean=mean(LSC_center_x);
LSC_center_ymean=mean(LSC_center_y);
LSC_center_x=(LSC_center_x-LSC_center_xmean).^2;
LSC_center_y=(LSC_center_y-LSC_center_ymean).^2;

Rg_1e9=sqrt((sum(LSC_center_x(:))+sum(LSC_center_y(:)))/length(LSC_center_y))/params.length0;
%%
inputDir = 'C:\Users\user\OneDrive\Desktop\01-RB-0612\fig12-vor\1e10\lsc_trajectory.plt';
nx=1025;
ny=nx;
constA=2.5;
Rayleigh=1e10;
Prandtl=0.71;
params = calculateSystemParameters(nx,ny, Rayleigh, Prandtl,constA,'log.log');

[LSC_center, ~] = read_tecplot_plt(inputDir);
LSC_center_x=LSC_center{1};
LSC_center_y=LSC_center{2};

save('LSC_1e10.mat', ...
     'LSC_center_x', ...
     'LSC_center_y');

LSC_center_xmean=mean(LSC_center_x);
LSC_center_ymean=mean(LSC_center_y);
LSC_center_x=(LSC_center_x-LSC_center_xmean).^2;
LSC_center_y=(LSC_center_y-LSC_center_ymean).^2;

Rg_1e10=sqrt((sum(LSC_center_x(:))+sum(LSC_center_y(:)))/length(LSC_center_y))/params.length0;
%%
inputDir = 'C:\Users\user\OneDrive\Desktop\01-RB-0612\fig12-vor\1e11\lsc_trajectory.plt';
nx=2049;
ny=nx;
constA=2.9;
Rayleigh=1e11;
Prandtl=0.71;
params = calculateSystemParameters(nx,ny, Rayleigh, Prandtl,constA,'log.log');

[LSC_center, ~] = read_tecplot_plt(inputDir);
LSC_center_x=LSC_center{1};
LSC_center_y=LSC_center{2};

save('LSC_1e11.mat', ...
     'LSC_center_x', ...
     'LSC_center_y');

LSC_center_xmean=mean(LSC_center_x);
LSC_center_ymean=mean(LSC_center_y);
LSC_center_x=(LSC_center_x-LSC_center_xmean).^2;
LSC_center_y=(LSC_center_y-LSC_center_ymean).^2;

Rg_1e11=sqrt((sum(LSC_center_x(:))+sum(LSC_center_y(:)))/length(LSC_center_y))/params.length0;
%%
inputDir = 'C:\Users\user\OneDrive\Desktop\01-RB-0612\fig12-vor\1e12\lsc_trajectory.plt';
nx=4097;
ny=nx;
constA=3.1;
Rayleigh=1e12;
Prandtl=0.71;
params = calculateSystemParameters(nx,ny, Rayleigh, Prandtl,constA,'log.log');

[LSC_center, ~] = read_tecplot_plt(inputDir);
LSC_center_x=LSC_center{1};
LSC_center_y=LSC_center{2};

save('LSC_1e12.mat', ...
     'LSC_center_x', ...
     'LSC_center_y');

LSC_center_xmean=mean(LSC_center_x);
LSC_center_ymean=mean(LSC_center_y);
LSC_center_x=(LSC_center_x-LSC_center_xmean).^2;
LSC_center_y=(LSC_center_y-LSC_center_ymean).^2;

Rg_1e12=sqrt((sum(LSC_center_x(:))+sum(LSC_center_y(:)))/length(LSC_center_y))/params.length0;
%%
Ra_list=[1e8,1e9,1e10,1e11,1e12];
Rg_list=[Rg_1e8,Rg_1e9,Rg_1e10,Rg_1e11,Rg_1e12];
plot_matlab(...
        'x', Ra_list, ...
        'y', Rg_list, ...
        'CanvasRatio', [10, 7], ...
        'Height', 540, ...
        'Units', 'points', ...
        'XLim',[5e7,2e12],...
        'YLim',[0,0.2],...
        'XTicks', [1e8,1e10, 1e12], ...
        'XLabel', '\it{Ra}', ...
        'YLabel', '\it{R_{g} /H}', ...
        'LogScale', 'x', ...          
        'LineColor', 'color1', ... 
        'LineStyle', 'none', ...          
        'Marker', {'o', 's'}, ...       
        'MarkerSize', 20, ...
        'MarkerInterval',1,...
        'ShowLegend', 0, ...
        'OutputFilename', 'Rg_Ra', ...
        'Resolution', 600);