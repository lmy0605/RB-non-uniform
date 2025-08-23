function data_tgt = conservative_remap_2d(x_src, y_src, data_src, x_tgt, y_tgt)
    nx_src = length(x_src);
    ny_src = length(y_src);
    nx_tgt = length(x_tgt);
    ny_tgt = length(y_tgt);
    num_fields = length(data_src);

    x_src = x_src(:)'; y_src = y_src(:); 
    x_tgt = x_tgt(:)'; y_tgt = y_tgt(:); 

    data_tgt = cell(size(data_src));
    for k = 1:num_fields
        data_tgt{k} = zeros(nx_tgt, ny_tgt);
    end
    % 预计算所有网格单元的边界
    x_src_bounds = calculate_cell_bounds(x_src);
    y_src_bounds = calculate_cell_bounds(y_src');
    x_tgt_bounds = calculate_cell_bounds(x_tgt);
    y_tgt_bounds = calculate_cell_bounds(y_tgt'); 

    grad_x_fields = cell(num_fields, 1);
    grad_y_fields = cell(num_fields, 1);
    for k = 1:num_fields
        [grad_x, grad_y] = compute_limited_gradients(data_src{k}, x_src, y_src);
        grad_x_fields{k} = grad_x;
        grad_y_fields{k} = grad_y;
    end

    for j = 1:ny_tgt      
        for i = 1:nx_tgt  
            % 确定可能与当前目标单元重叠的源单元的搜索范围
            i_start = max(1,      bounded_bisect(x_src, x_tgt_bounds(1, i)) - 1);
            i_end   = min(nx_src, bounded_bisect(x_src, x_tgt_bounds(2, i)) + 1);
            j_start = max(1,      bounded_bisect(y_src', y_tgt_bounds(1, j)) - 1);
            j_end   = min(ny_src, bounded_bisect(y_src', y_tgt_bounds(2, j)) + 1);

            integral_vals = zeros(1, num_fields);
            total_area = 0.0;

            for jj = j_start:j_end 
                for ii = i_start:i_end 
                   
                    [overlap_area, xc, yc] = calculate_overlap_properties(...
                        x_src_bounds(:, ii), y_src_bounds(:, jj), ...
                        x_tgt_bounds(:, i),  y_tgt_bounds(:, j));

                    if overlap_area > 1.0e-14
                        
                        for k = 1:num_fields
                            src_field = data_src{k};
                            grad_x = grad_x_fields{k};
                            grad_y = grad_y_fields{k};
                            % 二阶：使用梯度在重叠区域的中心点重建数值 
                            reconstructed_value = src_field(ii, jj) ...
                                      + grad_x(ii, jj) * (xc - x_src(ii)) ...
                                      + grad_y(ii, jj) * (yc - y_src(jj));

                       
                            integral_vals(k) = integral_vals(k) + reconstructed_value * overlap_area;
                        end
                        total_area = total_area + overlap_area;
                    end
                end
            end

            if total_area > 1.0e-12
                for k = 1:num_fields
                    data_tgt{k}(i, j) = integral_vals(k) / total_area;
                end
            else
                ii_fallback = bounded_bisect(x_src, x_tgt(i));
                jj_fallback = bounded_bisect(y_src', y_tgt(j));

                ii_fallback = min(max(ii_fallback, 1), nx_src);
                jj_fallback = min(max(jj_fallback, 1), ny_src);

                for k = 1:num_fields
                    data_tgt{k}(i, j) = data_src{k}(ii_fallback, jj_fallback);
                end
            end
        end
    end

end

function bounds = calculate_cell_bounds(arr_row)
    % 从单元中心点坐标计算一维单元的边界（左、右或上、下）
    n = length(arr_row);
    bounds = zeros(2, n); 

    bounds(1, 1) = arr_row(1) - 0.5 * (arr_row(2) - arr_row(1));
    bounds(1, 2:n) = 0.5 * (arr_row(1:n-1) + arr_row(2:n));

    bounds(2, 1:n-1) = 0.5 * (arr_row(1:n-1) + arr_row(2:n));
    bounds(2, n) = arr_row(n) + 0.5 * (arr_row(n) - arr_row(n-1));
end

function [grad_x, grad_y] = compute_limited_gradients(field, x, y)
    [nx, ny] = size(field);
    x = x(:)'; 
    y = y(:);  
    grad_x = zeros(nx, ny);
    grad_y = zeros(nx, ny);

    for i = 1:nx
        if ny > 1
            grad_y(i, 1) = (field(i, 2) - field(i, 1)) / (y(2) - y(1));
            grad_y(i, ny) = (field(i, ny) - field(i, ny-1)) / (y(ny) - y(ny-1));
            for j = 2:ny-1
                slope_B = (field(i, j) - field(i, j-1)) / (y(j) - y(j-1));
                slope_T = (field(i, j+1) - field(i, j)) / (y(j+1) - y(j));
                slope_cen = (field(i, j+1) - field(i, j-1)) / (y(j+1) - y(j-1));
                grad_y(i, j) = van_leer_limiter(slope_cen, slope_B, slope_T);
            end
        end
    end

    for j = 1:ny
        if nx > 1
            grad_x(1, j) = (field(2, j) - field(1, j)) / (x(2) - x(1));
            grad_x(nx, j) = (field(nx, j) - field(nx-1, j)) / (x(nx) - x(nx-1));
            for i = 2:nx-1
                slope_L = (field(i, j) - field(i-1, j)) / (x(i) - x(i-1));
                slope_R = (field(i+1, j) - field(i, j)) / (x(i+1) - x(i));
                slope_cen = (field(i+1, j) - field(i-1, j)) / (x(i+1) - x(i-1));
                grad_x(i, j) = van_leer_limiter(slope_cen, slope_L, slope_R);
            end
        end
    end
end

function limited_slope = van_leer_limiter(slope_cen, slope_1, slope_2)
    if slope_1 * slope_2 > 0.0
        r = slope_1 / slope_2;
        limited_slope = slope_cen * (r + abs(r)) / (1.0 + r*r);
    else
        limited_slope = 0.0;
    end
end

function [area, xc, yc] = calculate_overlap_properties(b_src_x, b_src_y, b_tgt_x, b_tgt_y)
    % 计算两个矩形单元重叠区域的面积和几何中心（质心）
    x_min = max(b_src_x(1), b_tgt_x(1));
    x_max = min(b_src_x(2), b_tgt_x(2));
    y_min = max(b_src_y(1), b_tgt_y(1));
    y_max = min(b_src_y(2), b_tgt_y(2));

    overlap_dx = x_max - x_min;
    overlap_dy = y_max - y_min;

    if overlap_dx > 0 && overlap_dy > 0
        area = overlap_dx * overlap_dy;
        xc = 0.5 * (x_min + x_max);
        yc = 0.5 * (y_min + y_max);
    else
        area = 0.0;
        xc = 0.0;
        yc = 0.0;
    end
end

function idx = bounded_bisect(arr, val)
    % 有界的二分搜索，用于快速查找一个值所在的区间索引
    n = length(arr);

    if val <= arr(1)
        idx = 1;
        return;
    elseif val >= arr(n)
        idx = n;
        return;
    end
    
    lo = 1;
    hi = n;
    while (hi - lo > 1)
        mid = floor((lo + hi) / 2);
        if val >= arr(mid)
            lo = mid;
        else
            hi = mid;
        end
    end
    idx = lo;
end