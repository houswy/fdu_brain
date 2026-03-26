function [is_place, is_object, SI_scores, MI_scores] = verify_cells_with_shuffle(behavior, imaging, object_pos)
    % 运行验证逻辑，复现识别位置细胞和对象细胞的部分
    % behavior: 行为数据
    % imaging: 钙成像数据 (temporal_weights)
    % object_pos: 对象在轨道上的位置 bin
    
    num_shuffles = 100; % 默认 100 次，根据老师要求可调整为 1000
    num_cells = size(imaging.temporal_weights, 2);
    num_frames = size(imaging.temporal_weights, 1);
    
    fprintf('正在计算原始指标...\n');
    % 计算原始 SI 和 MI
    [SI_scores, ~] = calculate_SI_internal(behavior, imaging);
    [MI_scores, ~] = calculate_MI_internal(behavior, imaging, object_pos);
    
    shuffled_SI = zeros(num_cells, num_shuffles);
    shuffled_MI = zeros(num_cells, num_shuffles);
    
    fprintf('正在进行 Shuffle 验证 (%d 次)...\n', num_shuffles);
    for s = 1:num_shuffles
        if mod(s, 10) == 0, fprintf('进度: %d%%\n', (s/num_shuffles)*100); end
        
        % 采用时间轴循环位移 (Circular Shift)
        % 这种方式保持了钙信号本身的动力学特征（如衰减常数）
        shift = randi([200, num_frames - 200]);
        shuffled_imaging = imaging;
        shuffled_imaging.temporal_weights = circshift(imaging.temporal_weights, shift, 1);
        
        shuffled_SI(:, s) = calculate_SI_internal(behavior, shuffled_imaging);
        shuffled_MI(:, s) = calculate_MI_internal(behavior, shuffled_imaging, object_pos);
    end
    
    % P-value 计算与显著性判定
    is_place = false(1, num_cells);
    is_object = false(1, num_cells);
    
    for i = 1:num_cells
        % Z-test 或直接计算 P 值
        p_SI = sum(shuffled_SI(i, :) >= SI_scores(i)) / num_shuffles;
        p_MI = sum(shuffled_MI(i, :) >= MI_scores(i)) / num_shuffles;
        
        is_place(i) = p_SI < 0.05;
        is_object(i) = p_MI < 0.05;
    end
    
    fprintf('验证完成！识别出 %d 个位置细胞, %d 个对象细胞。\n', sum(is_place), sum(is_object));
end

function [SI, p_val] = calculate_SI_internal(behavior, imaging)
    % 内部 SI 计算逻辑，参考 placecell_py.m
    num_cells = size(imaging.temporal_weights, 2);
    p_val = ones(1, num_cells); % 初始化返回值
    
    % 提取分箱活动
    total_bin = [];
    for i = 1:length(behavior)
        if behavior(i).trial_type == 2
            total_bin = [total_bin; behavior(i).location_bin];
        end
    end
    
    if isempty(total_bin)
        SI = zeros(1, num_cells);
        return;
    end
    
    use_bin = unique(total_bin);
    if length(use_bin) > 4
        use_bin = use_bin(3:end-2); % 剔除两端
    end
    
    p_bin = zeros(1, length(use_bin));
    for i = 1:length(use_bin)
        p_bin(i) = sum(total_bin == use_bin(i)) / length(total_bin);
    end
    
    SI = zeros(1, num_cells);
    for i = 1:num_cells
        cell_activity = imaging.temporal_weights(:, i);
        bin_activity = zeros(1, length(use_bin));
        
        % 确保活动数据与行为数据对齐
        data_len = min(length(cell_activity), length(total_bin));
        temp_activity = cell_activity(1:data_len);
        temp_total_bin = total_bin(1:data_len);
        
        for j = 1:length(use_bin)
            idx = temp_total_bin == use_bin(j);
            if any(idx)
                bin_activity(j) = mean(temp_activity(idx));
            else
                bin_activity(j) = 0;
            end
        end
        
        r_mean = mean(temp_activity);
        if r_mean > 0
            for j = 1:length(use_bin)
                if bin_activity(j) > 0
                    SI(i) = SI(i) + p_bin(j) * (bin_activity(j)/r_mean) * log2(bin_activity(j)/r_mean + eps);
                end
            end
        end
    end
end

function [MI, p_val] = calculate_MI_internal(behavior, imaging, object_pos)
    % 内部 MI 计算逻辑，参考 objectcell_MINI_psk.m
    num_cells = size(imaging.temporal_weights, 2);
    p_val = ones(1, num_cells); % 初始化返回值
    
    % 构造对象位置的二值序列 (0 或 1)
    total_bin = [];
    for i = 1:length(behavior)
        total_bin = [total_bin; behavior(i).location_bin];
    end
    
    if isempty(total_bin)
        MI = zeros(1, num_cells);
        return;
    end
    
    % 假设 object_pos 是物体所在的 bin
    object_vector = ismember(total_bin, [object_pos-1, object_pos, object_pos+1]);
    
    MI = zeros(1, num_cells);
    for i = 1:num_cells
        data_len = min(size(imaging.temporal_weights, 1), length(object_vector));
        cell_activity = imaging.temporal_weights(1:data_len, i);
        temp_object_vector = double(object_vector(1:data_len));
        
        if std(cell_activity) > 0 && std(temp_object_vector) > 0
            MI(i) = corr(cell_activity, temp_object_vector);
        end
    end
end
