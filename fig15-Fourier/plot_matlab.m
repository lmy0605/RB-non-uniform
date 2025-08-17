function plot_matlab(varargin)
% PLOT_MATLAB - A standardized, high-resolution plotting function with multiple data sources.
% 语法:
%   plot_matlab(x_data, y_data, 'ParameterName', ParameterValue, ...)
%   plot_matlab('DataSource', 'file', 'Filename', 'data.plt', ...)
%
% 详细说明:
%   此函数是一个可配置的绘图引擎，它支持两种主要的数据源：
%
%   1. 'variable' (默认): 直接从MATLAB工作区的变量（如向量、矩阵）绘图。
%   2. 'file': 从外部数据文件（支持Tecplot二进制.plt文件和简单的ASCII文本文件）
%      读取数据并绘图。
%
% 主要功能:
%   - 自动数据处理: 轻松处理单个或多个Y数据系列。
%   - 文件读取: 能够直接从Tecplot .plt文件或文本文件读取数据。
%   - 数据选择: 从文件中绘图时，可通过列索引、变量名或甚至
%     数学表达式（例如 'YVariable', 'Velocity*0.5'）来选择X和Y数据。
%   - 外观: 通过参数-值对，可以控制几乎所有绘图元素，
%     包括线型、颜色、标记、坐标轴范围、标签等。
%   - 预设样式: 内置了一套默认的颜色和样式，以确保图表的美观和一致性。
%     可通过 'LineColor', 'color1' 等方式调用。
%   - 输出: 将生成的图表保存为指定分辨率的PNG图像。
%% 1. 输入预处理器
if nargin > 0 && isnumeric(varargin{1})
    if nargin < 2 || ~(isnumeric(varargin{2}) || iscell(varargin{2})), error('When the first argument is numeric (x-data), the second argument must be y-data (numeric or cell).'); end
    processed_varargin = [{'x', varargin{1}, 'y', varargin{2}}, varargin(3:end)];
else, processed_varargin = varargin; end

%% 2. 输入解析器
p = inputParser;
p.KeepUnmatched = true; p.PartialMatching = false;
validateXInput = @(x) isnumeric(x) || iscell(x);
validateYInput = @(y) isnumeric(y) || iscell(y);
validateLabels = @(c) iscell(c) && all(cellfun(@(s) isstring(s) || ischar(s), c));
validateInterval = @(v) (isnumeric(v) && all(v > 0) && all(mod(v,1)==0)) || (iscell(v) && all(cellfun(@(x) isnumeric(x) && isscalar(x) && x > 0 && mod(x,1) == 0, v)));
validateLogScale = @(s) any(validatestring(s, {'none', 'x', 'y', 'xy', 'yx'}));
validateLogical = @(b) islogical(b) || (isnumeric(b) && isscalar(b) && (b==0 || b==1));
validateFilename = @(f) (ischar(f) || isstring(f)) || (iscell(f) && all(cellfun(@(s) ischar(s) || isstring(s), f)));
validateVarSpec = @(v) isempty(v) || (ischar(v) || isstring(v)) || (isnumeric(v) && isvector(v)) || iscell(v);
p.addParameter('DataSource', 'variable', @(s) any(validatestring(s, {'variable', 'file'})));
% p.addParameter('x', [], @isnumeric);
p.addParameter('x', [], validateXInput);
p.addParameter('y', [], validateYInput);
p.addParameter('LegendLabels', {}, validateLabels);
p.addParameter('ShowLegend', true, validateLogical);
p.addParameter('MarkerInterval', 5, validateInterval);
p.addParameter('LogScale', 'none', validateLogScale);
p.addParameter('Filename', '', validateFilename);
p.addParameter('XVariable', [], validateVarSpec);
p.addParameter('YVariable', [], validateVarSpec);
p.addParameter('CanvasRatio', [10, 7], @(r) isnumeric(r) && isvector(r) && numel(r) == 2 && all(r > 0));
p.addParameter('XLim', [], @(lim) isnumeric(lim) && isvector(lim) && numel(lim) == 2);
p.addParameter('YLim', [], @(lim) isnumeric(lim) && isvector(lim) && numel(lim) == 2);
p.addParameter('XTickInterval', [], @(s) isnumeric(s) && isscalar(s) && s > 0);
p.addParameter('YTickInterval', [], @(s) isnumeric(s) && isscalar(s) && s > 0);
p.addParameter('XLabel', '', @(s) isstring(s) || ischar(s));
p.addParameter('YLabel', '', @(s) isstring(s) || ischar(s));
p.addParameter('OutputFilename', 'plot_output.png', @(s) isstring(s) || ischar(s));
p.addParameter('Resolution', 600, @(s) isnumeric(s) && isscalar(s) && s > 0);
p.addParameter('LabelInterpreter', 'tex', @(s) any(validatestring(s, {'tex', 'latex', 'none'})));
p.addParameter('AxisLayer', 'top', @(s) any(validatestring(s, {'top', 'bottom'})));
p.addParameter('TickDirection', 'out', @(s) any(validatestring(s, {'in', 'out'})));
p.addParameter('Padding', false, validateLogical);
p.addParameter('XTicks', [], @isnumeric);

p.parse(processed_varargin{:});
params = p.Results;
unmatched_styles = p.Unmatched;

if ~isempty(params.Filename) && strcmp(params.DataSource, 'variable')
    params.DataSource = 'file';
end

%% 3. 数据获取与验证
x_data_cell = {}; y_data_cell = {}; y_labels_cell = {};
x_label_default = 'x'; y_label_default = '';

if strcmpi(params.DataSource, 'variable')
    if isempty(params.x) || isempty(params.y), error('Data source is "variable", but x and/or y data are missing.'); end
    
    % --- 标准化 Y 数据 ---
    % 无论输入是什么，都将Y转换为一个元胞数组 temp_y_cell
    if isnumeric(params.y)
        if isvector(params.y)
            temp_y_cell = {params.y}; % 单个向量
        else
            temp_y_cell = num2cell(params.y, 1); % 矩阵的每一列是一条曲线
        end
    elseif iscell(params.y)
        temp_y_cell = params.y; % 已经是元胞数组
    end
    num_curves = numel(temp_y_cell);

    % --- 标准化 X 数据，使其与 Y 匹配 ---
    % 无论输入是什么，都将X转换为一个包含 num_curves 个元素的元胞数组 temp_x_cell
    if isnumeric(params.x)
        % 如果 X 是单个数值数组，则为每条 Y 曲线复制它
        temp_x_cell = repmat({params.x(:)}, 1, num_curves);
    elseif iscell(params.x)
        if numel(params.x) == 1
            % 如果 X 是单元素的元胞，也为每条 Y 曲线复制它
            temp_x_cell = repmat({params.x{1}(:)}, 1, num_curves);
        elseif numel(params.x) ~= num_curves
            % 如果 X 是多元素的元胞，其数量必须与 Y 匹配
            error('If x is a cell array with multiple elements, its size must match the number of y-curves.');
        else
            % X 和 Y 一一对应
            temp_x_cell = params.x;
        end
    end
    
    % --- 验证并填充最终的数据元胞 ---
    for k = 1:num_curves
        x_vec = temp_x_cell{k}(:);
        y_vec = temp_y_cell{k}(:);
        
        if ~isnumeric(x_vec) || ~isvector(x_vec), error('Each x-data series must be a numeric vector.'); end
        if ~isnumeric(y_vec) || ~isvector(y_vec), error('Each y-data series must be a numeric vector.'); end
        if length(y_vec) ~= length(x_vec), error('Each y-vector must have the same length as its corresponding x-vector.'); end
        
        x_data_cell{end+1} = x_vec; 
        y_data_cell{end+1} = y_vec;
    end
    y_labels_cell = arrayfun(@(n) sprintf('Curve %d', n), 1:num_curves, 'UniformOutput', false);
else % strcmpi(params.DataSource, 'file')
    % --- Data from one or more files ---
    if isempty(params.Filename), error('DataSource is ''file'', but ''Filename'' was not specified.'); end
    if ischar(params.Filename) || isstring(params.Filename), file_list = {params.Filename}; else, file_list = params.Filename; end

    for file_idx = 1:numel(file_list)
        current_filename = file_list{file_idx};
        if ~exist(current_filename, 'file'), error('File does not exist: %s', current_filename); end

        is_ascii_fallback = false;
        try, [data_cell, var_names] = read_tecplot_plt(current_filename);
        catch, try, ascii_data = readmatrix(current_filename, 'FileType', 'text'); if isempty(ascii_data), error('ASCII file is empty or contains no numeric data.'); end, [~, num_cols] = size(ascii_data); var_names = arrayfun(@(n) sprintf('Var%d', n), 1:num_cols, 'UniformOutput', false); data_cell = num2cell(ascii_data, 1); is_ascii_fallback = true; catch ME, error('plot_matlab:FileReadError', 'Failed to read file "%s" as either Tecplot binary or simple ASCII. Error: %s', current_filename, ME.message); end, end
        num_vars = numel(var_names);

        % --- Process X variable (can be name or expression) ---
        current_x_var_spec = params.XVariable; if iscell(params.XVariable), if numel(params.XVariable) < file_idx, error('Not enough XVariable entries for the number of files.'); end, current_x_var_spec = params.XVariable{file_idx}; end
        if isempty(current_x_var_spec), current_x_var_spec = 1; end % Default to first column

        x_data_for_file = []; x_label_for_file = '';
        if isnumeric(current_x_var_spec)
            if current_x_var_spec > num_vars, error('XVariable index %d is out of bounds for file %s.', current_x_var_spec, current_filename); end
            x_idx = current_x_var_spec;
            x_data_for_file = data_cell{x_idx}; x_label_for_file = var_names{x_idx};
        else % It's a string, can be a name or an expression
            x_idx = find(strcmpi(var_names, current_x_var_spec), 1);
            if ~isempty(x_idx)
                x_data_for_file = data_cell{x_idx}; x_label_for_file = var_names{x_idx};
            else
                try
                    [x_data_for_file, x_label_for_file] = evaluate_expression(current_x_var_spec, data_cell, var_names);
                catch ME
                    error('Failed to evaluate XVariable expression "%s" for file %s. Error: %s', current_x_var_spec, current_filename, ME.message);
                end
            end
        end
        if file_idx == 1, x_label_default = x_label_for_file; end
        
        % --- Process Y variable(s) (can be name, index, or expression) ---
        current_y_var_spec = params.YVariable; if iscell(params.YVariable), if numel(params.YVariable) < file_idx, error('Not enough YVariable entries for the number of files.'); end, current_y_var_spec = params.YVariable{file_idx}; end
        
        y_specs_to_process = {};
        if isempty(current_y_var_spec), y_specs_to_process = setdiff(1:num_vars, x_idx);
        elseif iscell(current_y_var_spec), y_specs_to_process = current_y_var_spec;
        else, y_specs_to_process = {current_y_var_spec}; end
        
        for i = 1:numel(y_specs_to_process)
            spec = y_specs_to_process{i};
            y_data = []; y_label = '';
            
            if isnumeric(spec)
                if spec > num_vars, warning('YVariable index %d is out of bounds for file %s. Skipping.', spec, current_filename); continue; end
                y_data = data_cell{spec}; y_label = var_names{spec};
            elseif ischar(spec) || isstring(spec)
                y_idx = find(strcmpi(var_names, spec), 1);
                if ~isempty(y_idx)
                    y_data = data_cell{y_idx}; y_label = var_names{y_idx};
                else
                    try
                        [y_data, y_label] = evaluate_expression(spec, data_cell, var_names);
                    catch ME
                        warning('Could not evaluate YVariable expression "%s" for file %s. Skipping. Error: %s', spec, current_filename, ME.message);
                        continue;
                    end
                end
            else
                warning('Invalid YVariable specification. Skipping.'); continue;
            end
            
            x_data_cell{end+1} = x_data_for_file;
            y_data_cell{end+1} = y_data;
            y_labels_cell{end+1} = y_label;
        end
    end
end


% Final check and application of user-provided legend labels
if ~isempty(params.LegendLabels)
    if numel(params.LegendLabels) ~= numel(y_data_cell)
        error('The number of provided ''LegendLabels'' (%d) must match the total number of curves being plotted (%d).', numel(params.LegendLabels), numel(y_data_cell));
    end
    y_labels_cell = params.LegendLabels;
end

% Set final axis labels
x_label = params.XLabel; if isempty(x_label), x_label = x_label_default; end
y_label = params.YLabel; if isempty(y_label), y_label = y_label_default; end


%% 4. 创建图窗和坐标轴
fig = figure('Visible', 'on');
pos = get(fig, 'Position');
set(fig, 'Position', [pos(1), pos(2), pos(3), pos(3) * (params.CanvasRatio(2) / params.CanvasRatio(1))]);
ax = gca;

%% 5. 绘制图形
hold(ax, 'on');
default_colors = {[200, 36, 35]/255;[52, 128, 184]/255;[255, 190, 122]/255;[173, 211, 226]/255;[250, 136, 120]/255;[141, 206, 200]/255;[155, 191, 138]/255;[130, 175, 218]/255;[247, 144, 89]/255;[231, 219, 211]/255;[194, 189, 222]/255};
default_styles.LineWidth = 1.5; default_styles.MarkerSize = 8; 
%default_styles.MarkerFaceColor = 'w';
style_map = { 'LineColor', 'Color'; 'LineStyle', 'LineStyle'; 'LineWidth', 'LineWidth'; 'Marker', 'Marker'; 'MarkerSize', 'MarkerSize'; 'MarkerEdgeColor', 'MarkerEdgeColor'; 'MarkerFaceColor', 'MarkerFaceColor' };
for i = 1:numel(y_data_cell)
    current_plot_styles = default_styles;
    
    % --- Style processing loop ---
    for j = 1:size(style_map, 1)
        user_prop_name = style_map{j, 1};
        if isfield(unmatched_styles, user_prop_name)
            style_val = unmatched_styles.(user_prop_name);
            
            if iscell(style_val)
                raw_val = style_val{mod(i-1, numel(style_val)) + 1};
            else
                raw_val = style_val;
            end
            resolved_val = raw_val;
            
            % <-- Simplified logic to only handle 'colorN' -->
            if strcmpi(user_prop_name, 'LineColor') && (ischar(raw_val) || isstring(raw_val))
                tokens = regexp(lower(char(raw_val)), '^color(\d+)$', 'tokens');
                if ~isempty(tokens)
                    idx = str2double(tokens{1}{1});
                    if idx > 0 && idx <= numel(default_colors)
                        resolved_val = default_colors{idx};
                    end
                end
            end
            current_plot_styles.(user_prop_name) = resolved_val;
        end
    end
    
    % Apply default color if not specified by user
    if ~isfield(current_plot_styles, 'LineColor')
        current_plot_styles.LineColor = default_colors{mod(i-1, numel(default_colors)) + 1};
    end

    if ~isfield(current_plot_styles, 'MarkerEdgeColor')
        current_plot_styles.MarkerEdgeColor = current_plot_styles.LineColor;
    end
    
    if ~isfield(current_plot_styles, 'MarkerFaceColor')
    current_plot_styles.MarkerFaceColor = current_plot_styles.MarkerEdgeColor;
    end
    % Build plot arguments
    plot_args = {};
    for j = 1:size(style_map, 1)
        user_prop_name = style_map{j, 1};
        matlab_prop_name = style_map{j, 2};
        if isfield(current_plot_styles, user_prop_name)
            plot_args = [plot_args, {matlab_prop_name, current_plot_styles.(user_prop_name)}];
        end
    end
    
    % MarkerInterval logic and plot command
    current_interval = 1;
    if iscell(params.MarkerInterval), current_interval = params.MarkerInterval{mod(i-1, numel(params.MarkerInterval)) + 1};
    elseif numel(params.MarkerInterval) > 1, current_interval = params.MarkerInterval(mod(i-1, numel(params.MarkerInterval)) + 1);
    else, current_interval = params.MarkerInterval; end
    if current_interval > 1
        marker_indices = 1:current_interval:numel(x_data_cell{i});
        plot_args = [plot_args, {'MarkerIndices', marker_indices}];
    end
    plot(ax, x_data_cell{i}, y_data_cell{i}, plot_args{:});
end
% --- 将循环顺序从 1:N 修改为 N:-1:1 ---
% for i = numel(y_data_cell):-1:1
%     current_plot_styles = default_styles;
%     
%     % --- Style processing loop ---
%     % 循环体内部的所有代码都保持不变
%     for j = 1:size(style_map, 1)
%         user_prop_name = style_map{j, 1};
%         if isfield(unmatched_styles, user_prop_name)
%             style_val = unmatched_styles.(user_prop_name);
%             
%             if iscell(style_val)
%                 % 注意：这里也需要调整索引以匹配倒序
%                 raw_val = style_val{mod(i-1, numel(style_val)) + 1};
%             else
%                 raw_val = style_val;
%             end
%             
%             resolved_val = raw_val;
%             
%             % <-- Simplified logic to only handle 'colorN' -->
%             if strcmpi(user_prop_name, 'LineColor') && (ischar(raw_val) || isstring(raw_val))
%                 tokens = regexp(lower(char(raw_val)), '^color(\d+)$', 'tokens');
%                 if ~isempty(tokens)
%                     idx = str2double(tokens{1}{1});
%                     if idx > 0 && idx <= numel(default_colors)
%                         resolved_val = default_colors{idx};
%                     end
%                 end
%             end
%             current_plot_styles.(user_prop_name) = resolved_val;
%         end
%     end
%     
%     % Apply default color if not specified by user
%     if ~isfield(current_plot_styles, 'LineColor')
%         current_plot_styles.LineColor = default_colors{mod(i-1, numel(default_colors)) + 1};
%     end
% 
%     if ~isfield(current_plot_styles, 'MarkerEdgeColor')
%         current_plot_styles.MarkerEdgeColor = current_plot_styles.LineColor;
%     end
%     
%     if ~isfield(current_plot_styles, 'MarkerFaceColor')
%     current_plot_styles.MarkerFaceColor = current_plot_styles.MarkerEdgeColor;
%     end
%     % Build plot arguments
%     plot_args = {};
%     for j = 1:size(style_map, 1)
%         user_prop_name = style_map{j, 1};
%         matlab_prop_name = style_map{j, 2};
%         if isfield(current_plot_styles, user_prop_name)
%             plot_args = [plot_args, {matlab_prop_name, current_plot_styles.(user_prop_name)}];
%         end
%     end
%     
%     % MarkerInterval logic and plot command
%     if iscell(params.MarkerInterval), current_interval = params.MarkerInterval{mod(i-1, numel(params.MarkerInterval)) + 1};
%     elseif numel(params.MarkerInterval) > 1, current_interval = params.MarkerInterval(mod(i-1, numel(params.MarkerInterval)) + 1);
%     else, current_interval = params.MarkerInterval; end
%     if current_interval > 1
%         marker_indices = 1:current_interval:numel(x_data_cell{i});
%         plot_args = [plot_args, {'MarkerIndices', marker_indices}];
%     end
%     plot(ax, x_data_cell{i}, y_data_cell{i}, plot_args{:});
% end
hold(ax, 'off');

%% 6. 自定义坐标轴和标签
set(ax, 'FontName', 'Times', 'FontSize', 22, 'LineWidth', 1.2, 'Box', 'on');
set(ax, 'Layer', 'top');

if contains(lower(params.LogScale), 'x'), if any(cellfun(@(x) any(x <= 0), x_data_cell)), warning('plot_matlab:LogData', 'X-data contains non-positive values, which will be ignored on a log scale.'); end, set(ax, 'XScale', 'log'); end
if contains(lower(params.LogScale), 'y'), if any(cellfun(@(y) any(y <= 0), y_data_cell)), warning('plot_matlab:LogData', 'One or more Y-data series contain non-positive values, which will be ignored on a log scale.'); end, set(ax, 'YScale', 'log'); end

% 'LabelInterpreter' parameter
xlabel(ax, x_label,  'Interpreter', params.LabelInterpreter);
if numel(y_data_cell) == 1 && isempty(y_label) && ~isempty(y_labels_cell), y_label = y_labels_cell{1}; end
ylabel(ax, y_label, 'Interpreter', params.LabelInterpreter);

if params.ShowLegend && numel(y_data_cell) > 1 && ~isempty(y_labels_cell)
    % 'LabelInterpreter' parameter in the legend as well
    legend(ax, y_labels_cell, 'Interpreter', params.LabelInterpreter, 'Location', 'best');
end

if ~isempty(params.XLim), xlim(ax, params.XLim); end, if ~isempty(params.YLim), ylim(ax, params.YLim); end
if ~isempty(params.XTicks)
    xticks(ax, params.XTicks);
% 保留原来的 TickInterval 逻辑，但用 elseif
elseif ~strcmpi(get(ax, 'XScale'), 'log') && ~isempty(params.XTickInterval)
    current_xlim = xlim(ax); 
    xticks(ax, current_xlim(1):params.XTickInterval:current_xlim(2)); 
end
if ~strcmpi(get(ax, 'XScale'), 'log') && ~isempty(params.XTickInterval), current_xlim = xlim(ax); xticks(ax, current_xlim(1):params.XTickInterval:current_xlim(2)); end
if ~strcmpi(get(ax, 'YScale'), 'log') && ~isempty(params.YTickInterval), current_ylim = ylim(ax); yticks(ax, current_ylim(1):params.YTickInterval:current_ylim(2)); end
grid off;
%% 7. 输出图片
print(fig, params.OutputFilename, '-dpng', ['-r' num2str(params.Resolution)]);
fprintf('successfully save: %s\n', params.OutputFilename);
% close(fig);

end



% --- Helper function for evaluating expressions ---
function [result, label] = evaluate_expression(expression_str, data_cell, var_names)
    % Takes an expression string (e.g., "U*0.5 + V"), the cell array of
    % data vectors, and the cell array of variable names.
    % Returns the calculated data vector and the expression string as a label.

    % Find all words in the expression that could be variable names
    potential_vars = regexp(expression_str, '[a-zA-Z_]\w*', 'match');
    
    % Find which of these potential variables actually exist in the file
    vars_in_expr = intersect(potential_vars, var_names);
    
    if isempty(vars_in_expr)
        % Maybe it's an expression with just numbers, e.g., '3*pi/2'
        % To be safe, let's just try to evaluate it. If it fails, error out.
        try
            result = eval(expression_str);
            % It might evaluate to a scalar, so we need to make it a vector
            % matching the size of the first data column.
            if isscalar(result)
                result = result * ones(size(data_cell{1}));
            end
            label = expression_str;
            return;
        catch
             error('plot_matlab:ExpressionError', 'Expression "%s" contains no known variables from the file.', expression_str);
        end
    end
    
    % Prepare the expression for safe, element-wise evaluation
    safe_expr = expression_str;
    % Replace standard operators with their element-wise counterparts
    safe_expr = strrep(safe_expr, '*', '.*');
    safe_expr = strrep(safe_expr, '/', './');
    safe_expr = strrep(safe_expr, '^', '.^');
    
    % Build the argument list for the anonymous function
    arg_list_str = strjoin(vars_in_expr, ',');
    
    % Create the anonymous function handle from the string
    eval_func = str2func(['@(' arg_list_str ') ' safe_expr]);
    
    % Prepare the cell array of data vectors to pass to the function
    arg_data = cell(1, numel(vars_in_expr));
    for i = 1:numel(vars_in_expr)
        var_name = vars_in_expr{i};
        idx = find(strcmpi(var_names, var_name), 1);
        arg_data{i} = data_cell{idx};
    end
    
    % Execute the function with the data
    result = eval_func(arg_data{:});
    label = expression_str; % Use the original expression as the label
end

function [data_cell, var_names, file_info] = read_tecplot_plt(filename)
% read_tecplot_plt_binary Reads binary Tecplot .plt files created by 'liton_ordered_tec'.
%   This function is specifically designed to parse the binary format
%   written by the user's provided TEC_FILE class.
%
%   INPUT:
%       filename - A string containing the name of the binary .plt file.
%
%   OUTPUT:
%       data_cell - A 1xN cell array where N is the number of variables.
%                   Each cell contains a column vector of data for one variable.
%       var_names - A cell array of strings containing the variable names.
%       file_info - A struct containing metadata like Title, ZoneName, etc.

% --- 1. File Opening and Setup ---
fid = fopen(filename, 'rb'); % 'rb' for "read binary"
if fid == -1
    error('read_tecplot_plt_binary:CannotOpenFile', 'Cannot open file: %s', filename);
end
% Ensure the file is closed automatically when the function exits
cleanupObj = onCleanup(@() fclose(fid));

file_info = struct();

% --- 2. Read File Header ---
% i. Magic number, Version number
magic = fread(fid, 8, '*char')';
if ~strcmp(magic, '#!TDV112')
    error('read_tecplot_plt_binary:InvalidMagic', 'File is not a valid Tecplot binary file. Magic number is incorrect.');
end
file_info.Magic = magic;

% ii. Byte order check
fread(fid, 1, 'int32'); 

% iii. Title and variable names
file_info.FileType = fread(fid, 1, 'int32');
file_info.Title = read_null_terminated_string(fid);
num_vars = fread(fid, 1, 'int32');
var_names = cell(1, num_vars);
for i = 1:num_vars
    var_names{i} = read_null_terminated_string(fid);
end
file_info.Variables = var_names;

% --- 3. Read Zone Header ---
% This format only supports one zone, so we don't need a loop.
zone_marker = fread(fid, 1, 'float32');
if abs(zone_marker - 299.0) > 1e-6
    warning('read_tecplot_plt_binary:ZoneMarker', 'Zone marker is not 299.0 as expected.');
end

file_info.ZoneName = read_null_terminated_string(fid);
fread(fid, 1, 'int32'); % ParentZone
file_info.StrandId = fread(fid, 1, 'int32');
file_info.SolutionTime = fread(fid, 1, 'float64');
fread(fid, 1, 'int32'); % Not used
file_info.ZoneType = fread(fid, 1, 'int32');
fread(fid, 3, 'int32'); % VarLocation, RawFaceNeighbors, NumMiscFaceConnections

% Read IMax, JMax, KMax
zone_dims = fread(fid, 3, 'int32');
IMax = zone_dims(1);
JMax = zone_dims(2);
KMax = zone_dims(3);
num_points = IMax * JMax * KMax; % Total data points per variable

% Skip Auxiliary Data as it's not used in the user's case for reading
% But the writer supports it, so we must be able to skip it.
% The writer code writes "int32(1)" to indicate an aux pair, and "int32(0)" to end.
while true
    indicator = fread(fid, 1, 'int32');
    if indicator == 0 % No more aux data
        break;
    elseif indicator == 1 % Aux data pair follows
        read_null_terminated_string(fid); % Read and discard name
        fread(fid, 1, 'int32'); % Read and discard value format
        read_null_terminated_string(fid); % Read and discard value
    else
        % We've likely misread the stream, move back and assume end of aux
        fseek(fid, -4, 'cof');
        break;
    end
end


% --- 4. Read End of Header Marker ---
eoh_marker = fread(fid, 1, 'float32');
if abs(eoh_marker - 357.0) > 1e-6
    error('read_tecplot_plt_binary:EOHMarker', 'End of Header marker is not 357.0. Header parsing failed.');
end

% --- 5. Read Data Section ---
% i. Zone Data Header
data_zone_marker = fread(fid, 1, 'float32');
if abs(data_zone_marker - 299.0) > 1e-6
    warning('read_tecplot_plt_binary:DataZoneMarker', 'Data section zone marker is not 299.0.');
end

% ii. Read data types for each variable
% From the writer's gettype function: 1=single, 2=double, 3=int32, etc.
type_map = {'single', 'double', 'int32', 'int16', 'int8', 'bit1'};
data_types_int = zeros(1, num_vars);
data_types_str = cell(1, num_vars);
for i = 1:num_vars
    type_code = fread(fid, 1, 'int32');
    data_types_int(i) = type_code;
    data_types_str{i} = type_map{type_code};
end

% iii. Skip passive/shared variables info
fread(fid, 3, 'int32'); % Passive, Shared, Conn-Share

% iv. Skip Min/Max values (read 2 doubles for each variable)
fseek(fid, 2 * 8 * num_vars, 'cof');

% v. Read the actual zone data
fprintf('reading %d variables, each contains %d points \n', num_vars, num_points);
data_cell = cell(1, num_vars);
for i = 1:num_vars
    % The '*' prefix tells fread to return data in its native class
    precision = ['*', data_types_str{i}];
    data_cell{i} = fread(fid, num_points, precision);
end

end


% --- Helper function to read null-terminated strings ---
function str = read_null_terminated_string(fid)
    % Reads characters until a null character (ASCII 0) is encountered.
    str_chars = [];
    while true
        char_val = fread(fid, 1, 'int32'); % The writer uses int32 for characters
        if isempty(char_val) || char_val == 0
            break;
        end
        str_chars(end+1) = char_val;
    end
    str = char(str_chars);
end