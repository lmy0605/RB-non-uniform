clear; close all; clc;

%% 1. 参数设置 (Parameter Configuration)
% 将所有变化的参数存入数组，保持索引对应
Ra_list     = [1e8,   3e8,   1e9,   3e9,   1e10,  3e10,  1e11,  3e11,  1e12];
nx_list     = [257,   385,   513,   769,   1025,  1537,  2049,  3073,  4097];
constA_list = [1.5,   1.8,   2.1,   2.3,   2.5,   2.7,   2.9,   3.0,   3.1];

% 文件名后缀 (字符串形式，用于匹配文件名)
suffixes    = {'1e8', '3e8', '1e9', '3e9', '1e10','3e10','1e11','3e11','1e12'};

% 常量
Prandtl = 0.71;

% 预分配结果数组
Rg_list = zeros(size(Ra_list));

%% 2. 循环处理 (Loop Processing)
fprintf('开始处理数据...\n');

for i = 1:length(Ra_list)
    % --- 获取当前循环的参数 ---
    Rayleigh = Ra_list(i);
    nx = nx_list(i);
    ny = nx;
    constA = constA_list(i);
    suffix = suffixes{i};
    
    % --- 计算系统参数 ---
    % 注意：确保 calculateSystemParameters 函数在路径中
    params = calculateSystemParameters(nx, ny, Rayleigh, Prandtl, constA, 'log.log');
    
    % --- 构建文件路径并加载数据 ---
    filename = sprintf('../data/LSC_%s.mat', suffix);
    
    if isfile(filename)
        fprintf('正在读取: %s (Ra=%.1e, nx=%d)\n', filename, Rayleigh, nx);
        
        data = load(filename, 'LSC_center_x', 'LSC_center_y');
        
        % 获取数据
        LSC_x = data.LSC_center_x;
        LSC_y = data.LSC_center_y;
        
        % --- 计算 Rg (Gyration Radius) ---
        % 1. 计算均值
        LSC_xmean = mean(LSC_x);
        LSC_ymean = mean(LSC_y);
        
        % 2. 计算去均值后的平方 (Squared deviation)
        LSC_x_sq = (LSC_x - LSC_xmean).^2;
        LSC_y_sq = (LSC_y - LSC_ymean).^2;
        
        % 3. 计算 Rg 并归一化
        % 公式: sqrt( mean(dx^2 + dy^2) ) / length0
        mean_sq_dist = (sum(LSC_x_sq(:)) + sum(LSC_y_sq(:))) / length(LSC_y_sq);
        Rg_val = sqrt(mean_sq_dist) / params.length0;
        
        % 存储结果
        Rg_list(i) = Rg_val;
    else
        warning('文件 %s 不存在，跳过该数据点。', filename);
        Rg_list(i) = NaN;
    end
end
