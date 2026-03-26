function [input_strength, feedback_strength, classification] = analyze_interaction(hpc_imaging, ctx_imaging, behavior)
    % 分析皮层-海马交互，区分“视觉驱动”与“海马抽象”对象细胞
    % hpc_imaging: 海马钙成像数据
    % ctx_imaging: 新皮层钙成像数据
    % behavior: 行为数据
    
    num_ctx_cells = size(ctx_imaging.temporal_weights, 2);
    num_hpc_cells = size(hpc_imaging.temporal_weights, 2);
    
    % 1. 计算时间滞后相关性 (Cross-Correlation with Lag)
    % 目标是计算新皮层神经元与海马平均活动的相关性
    hpc_avg_activity = mean(hpc_imaging.temporal_weights, 2);
    
    input_strength = zeros(1, num_ctx_cells);
    feedback_strength = zeros(1, num_ctx_cells);
    classification = cell(1, num_ctx_cells);
    
    max_lag = 100; % 设定最大滞后帧数，对应大约 2-3 秒
    
    fprintf('正在计算时间滞后相关性与信号流向...\n');
    for i = 1:num_ctx_cells
        ctx_cell = ctx_imaging.temporal_weights(:, i);
        
        % 确保 ctx_cell 和 hpc_avg_activity 长度一致
        data_len = min(length(ctx_cell), length(hpc_avg_activity));
        temp_ctx = ctx_cell(1:data_len);
        temp_hpc = hpc_avg_activity(1:data_len);
        
        % 使用 xcov 计算互协方差，并手动进行归一化以避免长度不一导致的 'coeff' 报错
        [c, lags] = xcov(temp_ctx, temp_hpc, max_lag, 'none');
        
        % 手动归一化 (类似 'coeff' 选项)
        c = c / (norm(temp_ctx - mean(temp_ctx)) * norm(temp_hpc - mean(temp_hpc)));
        
        % 寻找峰值位置
        [max_corr, idx] = max(c);
        peak_lag = lags(idx);
        
        % 定义输入/输出强度
        % 这里的外部输入强度 (Input Strength) 简化为皮层神经元自身的独立活跃度/与物体的直接相关性
        % 实际分析中应结合外部视觉刺激信号
        input_strength(i) = var(ctx_cell); % 暂时用活跃度作为替代
        
        % 海马反馈强度 (Feedback Strength)
        % 如果 peak_lag > 0, 说明海马领先皮层 (HPC -> CTX, 反馈/抽象)
        % 如果 peak_lag < 0, 说明皮层领先海马 (CTX -> HPC, 视觉输入)
        if peak_lag > 0
            feedback_strength(i) = max_corr;
        else
            feedback_strength(i) = 0; % 仅考虑正向反馈强度
        end
        
        % 分类逻辑
        if peak_lag > 0 && max_corr > 0.3
            classification{i} = 'Hippocampal-Abstracted';
        elseif peak_lag < 0 && max_corr > 0.3
            classification{i} = 'Visual-Driven';
        else
            classification{i} = 'Independent';
        end
    end
    
    % 2. 可视化分类结果并保存
    fig2 = figure('Position', [100, 100, 800, 600], 'Visible', 'off');
    scatter(input_strength, feedback_strength, 40, 'filled');
    hold on;
    
    % 标记不同类型的点
    idx_abstract = strcmp(classification, 'Hippocampal-Abstracted');
    idx_visual = strcmp(classification, 'Visual-Driven');
    
    scatter(input_strength(idx_abstract), feedback_strength(idx_abstract), 60, 'r', 'filled');
    scatter(input_strength(idx_visual), feedback_strength(idx_visual), 60, 'b', 'filled');
    
    xlabel('外部输入强度 (Visual Input Strength)');
    ylabel('海马反馈强度 (Hippocampal Feedback Strength)');
    title('对象细胞分类: 视觉驱动 vs 海马抽象');
    legend('其他细胞', '海马抽象型 (反馈型)', '视觉驱动型 (前馈型)');
    grid on;
    
    saveas(fig2, 'object_cell_classification.png');
    close(fig2);
    
    fprintf('分析完成！已生成 2D 散点图。\n');
end
