% clear; close all; clc;
% format long
function params = calculateSystemParameters(nx,ny, Rayleigh, Prandtl,constA, logFileName)
% CALCULATESYSTEMPARAMETERS 计算系统参数
%   params = calculateSystemParameters(nx,ny, Rayleigh, Prandtl, constA, logFileName)
%   input:
%       nx          - x gridNumber
%       ny          - y gridNumber
%       Rayleigh
%       Prandtl 
%       constA    
%       logFileName - (Optional) Name of the log file to write parameters to.
%   output:
%       params      

% --- MODIFICATION START: Handle optional log file argument ---
if nargin < 6
    logFileName = ''; % If no log file is specified, set to empty
end

% Check if we should write to a file or display to console
writeToFile = ~isempty(logFileName);
fileID = -1; % Initialize file ID

if writeToFile
    fileID = fopen(logFileName, 'a'); % Open file for writing, 'w' overwrites existing file
    if fileID == -1
        error('Cannot open log file for writing: %s', logFileName);
    end
end
% --- MODIFICATION END ---

format long;

%% basic settings
% (The commented-out code remains the same)
nxHalf=(nx-1)/2+1;
nyHalf=(nx-1)/2+1;

Mach=0.1;

rho0 = 1.0;
Thot=1.0;
Tcold=0.0;
Tref=0.0d0;

%xGrid yGrid
denom = erf(0.5 * constA);  
i = 1:nx+1;  
xi = i / (nx+1);  
xGrid = 0.5 * (erf(constA * (xi - 0.5)) ./ denom + 1);
xGrid0=0.5 * (erf(constA * (0 - 0.5)) ./ denom + 1);

j = 1:ny+1;
yj = j / (ny+1);
yGrid = 0.5 * (erf(constA * (yj - 0.5)) ./ denom + 1);
yGrid0 = 0.5 * (erf(constA * (0 - 0.5)) ./ denom + 1);
%% 0 system parameters  '0' represents dimensionless
length0 = xGrid(nx);  
viscosity0 = Mach * length0 * sqrt(Prandtl) / sqrt(3.0*Rayleigh);

dx(2:nx) = xGrid(2:end-1) - xGrid(1:end-2);
dy(2:ny) = yGrid(2:end-1) - yGrid(1:end-2);
dx(1)=xGrid(1)-0;
dy(1)=yGrid(1)-0;
dy0 = yGrid(1)-yGrid0;
gridRatioQ = dx(2) / dx(1);
%% LB system parameters
length_LB = 1.0/ dy0 -1.0;  
viscosity_LB = Mach * length_LB * sqrt(Prandtl) / sqrt(3.0 * Rayleigh);
tauf = viscosity_LB * 3.0d0 + 0.5d0;
diffusivity_LB = viscosity_LB / Prandtl;
paraA = 20.0d0*sqrt(3.0d0) * diffusivity_LB - 4.0d0;
gBeta1 = Rayleigh * viscosity_LB * diffusivity_LB / length_LB;
gBeta  = gBeta1 / length_LB / length_LB;
time_LB = sqrt(length_LB / gBeta);
U_LB = sqrt(gBeta * length_LB);

timeUnit = time_LB; 
velocityUnit = U_LB;

% --- MODIFICATION START: Replaced disp with fprintf for logging ---
% This section now either prints to the file (if fileID is valid)
% or to the command window (if fileID is -1, which we can achieve
% by using fileID=1 for standard output).
output_target = 1; % Default to command window (stdout)
if writeToFile
    output_target = fileID; % If writing to file, use the file's ID
    fprintf(output_target, '--- System Parameters Log ---\n');
    fprintf(output_target, 'Date: %s\n', datestr(now));
    fprintf(output_target, '-----------------------------\n\n');
end

% The logic to either display or write to file
if writeToFile
    % Write to file using fprintf with the file ID
    fprintf(fileID, 'xGrid0 = %.15g\n', xGrid0);
    fprintf(fileID, 'xGrid(nx+1) = %.15g\n', xGrid(nx+1));
    fprintf(fileID, 'xGrid(1) = %.15g\n', xGrid(1));
    fprintf(fileID, 'xGrid(nx) = %.15g\n', xGrid(nx));
    fprintf(fileID, 'xGrid(1)+xGrid(nx) = %.15g\n\n', xGrid(1)+xGrid(nx));
    
    fprintf(fileID, 'length0 = %.15g\n', length0);
    fprintf(fileID, 'viscosity0 = %.15g\n\n', viscosity0);
    
    fprintf(fileID, 'dy0 = %.15g\n', dy0);
    fprintf(fileID, 'gridRatioQ = %.15g\n\n', gridRatioQ);
    
    fprintf(fileID, 'length_LB = %.15g\n', length_LB);
    fprintf(fileID, 'characteristic length (cell height) = %.15g l.u.\n', length_LB);
    fprintf(fileID, 'viscosity_LB = %.15g l.u.^2/t.s.\n', viscosity_LB);
    fprintf(fileID, 'tauf = %.15g\n', tauf);
    fprintf(fileID, 'diffusivity_LB = %.15g l.u.^2/t.s.\n', diffusivity_LB);
    fprintf(fileID, 'paraA = %.15g; should be between -4 and 1\n', paraA);
    fprintf(fileID, 'characteristic time = %.15g t.s.\n', time_LB);
    fprintf(fileID, 'characteristic velocity = %.15g l.u./t.s.\n', U_LB);

else
    % Original behavior: display to command window
    disp(['xGrid0 = ', num2str(xGrid0,15)])
    disp(['xGrid(nx+1) = ', num2str(xGrid(nx+1),15)])
    disp(['xGrid(1) = ', num2str(xGrid(1),15)])
    disp(['xGrid(nx) = ', num2str(xGrid(nx),15)])
    disp(['xGrid(1)+xGrid(nx) = ', num2str(xGrid(1)+xGrid(nx),15)])
    
    disp(' ')
    disp(['length0 = ', num2str(length0,15)])
    disp(['viscosity0 = ', num2str(viscosity0,15)])
    
    disp(' ')
    disp(['dy0 = ', num2str(dy0,15)])
    disp(['gridRatioQ = ', num2str(gridRatioQ,15)])
    
    disp(' ')
    disp(['length_LB = ', num2str(length_LB,15)])
    disp(['characteristic length (cell height) =', num2str(length_LB,15), 'l.u.'])
    disp(['viscosity_LB = ', num2str(viscosity_LB,15), 'l.u.^2/t.s.'])
    disp(['tauf = ', num2str(tauf,15)])
    disp(['diffusivity_LB = ', num2str(diffusivity_LB,15),'l.u.^2/t.s.'])
    disp(['paraA = ', num2str(paraA,15),'; should be between -4 and 1'])
    disp(['characteristic time =', num2str(time_LB,15), 't.s.'])
    disp(['characteristic velocity =', num2str(U_LB,15), 'l.u./t.s.'])
end

if writeToFile
    fclose(fileID); % Close the file
    disp(['Parameters successfully written to ', logFileName]);
end
% --- MODIFICATION END ---

%% save
params = struct(...
    'nxHalf',nxHalf,...
    'nyHalf',nyHalf,...
    'rho0',rho0,...
    'Thot',Thot,...
    'Tcold',Tcold,...
    'Tref',Tref,...
    'nx', nx, ...
    'ny', ny, ...
    'xGrid', xGrid, ...
    'yGrid', yGrid, ...
    'xGrid0', xGrid0, ...
    'yGrid0', yGrid0, ...
    'dx', dx, ...
    'dy', dy, ...
    'dy0', dy0, ...
    'length0', length0, ...
    'viscosity0', viscosity0, ...
    'length_LB', length_LB, ...
    'viscosity_LB', viscosity_LB, ...
    'tauf', tauf, ...
    'diffusivity_LB', diffusivity_LB, ...
    'paraA', paraA, ...
    'gBeta', gBeta, ...
    'time_LB', time_LB, ...
    'U_LB', U_LB, ...
    'timeUnit', timeUnit, ...
    'velocityUnit', velocityUnit, ...
    'gridRatioQ', gridRatioQ);
end