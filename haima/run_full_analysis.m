% 主运行脚本: 皮层-海马交互与对象细胞分类分析 (haima/run_full_analysis.m)

clear; clc; close all;

% 1. 加载数据
% 根据目录结构加载海马和新皮层的数据 (从父目录查找 Data_demo_gzx)
data_path = '../Data_demo_gzx/';
fprintf('加载海马 CA1 数据...\n');
load([data_path, 'Data_hippocampal_CA1/activity/output_sort_final_2.mat']);
hpc_imaging = output_sort_final_2;

fprintf('加载新皮层数据...\n');
load([data_path, 'Data_neocortex/activity/output_sort1_final_2.mat']);
ctx_imaging = output_sort1_final_2;

fprintf('加载行为数据...\n');
load([data_path, 'Data_neocortex/behavior/behavior_20251106_mode2.mat']);
behavior = behavior_20251106_mode2;

% 2. 自动优化参数 (确定最佳 Trial 筛选标准和步长)
fprintf('\n=== 正在进行参数自动优化 ===\n');
[best_speed_th, best_n_step] = optimize_parameters(behavior, ctx_imaging);

% 根据优化结果筛选 Trial
valid_trials_idx = [];
for t = 1:length(behavior)
    if isfield(behavior(t), 'speed')
        avg_speed = mean(behavior(t).speed);
    else
        avg_speed = std(behavior(t).location_bin); 
    end
    if avg_speed > best_speed_th
        valid_trials_idx = [valid_trials_idx, t];
    end
end
behavior_optimized = behavior(valid_trials_idx);
fprintf('优化完成：筛选出 %d 个高质量 Trial，建议 RNN 步长为 %d-step。\n', length(behavior_optimized), best_n_step);

% 3. 执行位置细胞和对象细胞的 Shuffle 验证 (使用优化后的 Trial)
% 假设对象位置在轨道 bin 15
object_pos = 15; 
fprintf('\n=== 开始识别位置细胞和对象细胞 (使用优化参数) ===\n');
[is_place_hpc, is_object_hpc, SI_hpc, MI_hpc] = verify_cells_with_shuffle(behavior_optimized, hpc_imaging, object_pos);
[is_place_ctx, is_object_ctx, SI_ctx, MI_ctx] = verify_cells_with_shuffle(behavior_optimized, ctx_imaging, object_pos);

% 4. 分析皮层-海马交互与信号流向 (结合最佳步长逻辑)
fprintf('\n=== 分析皮层-海马交互流向与分类 ===\n');
[input_strength, feedback_strength, classification] = analyze_interaction(hpc_imaging, ctx_imaging, behavior_optimized);

% 5. 结果汇总与可视化
fprintf('\n=== 汇总结果 ===\n');
fprintf('优化参数: 速度阈值 > %.1f, RNN 步长 = %d-step\n', best_speed_th, best_n_step);
fprintf('海马 CA1: 位置细胞 %d 个, 对象细胞 %d 个\n', sum(is_place_hpc), sum(is_object_hpc));
fprintf('新皮层: 位置细胞 %d 个, 对象细胞 %d 个\n', sum(is_place_ctx), sum(is_object_ctx));

% 分析新皮层对象细胞的分类分布
num_abstract = sum(strcmp(classification, 'Hippocampal-Abstracted'));
num_visual = sum(strcmp(classification, 'Visual-Driven'));
fprintf('新皮层对象细胞中: %d 个为海马抽象型, %d 个为视觉驱动型\n', num_abstract, num_visual);

% 保存分析结果
save('analysis_results.mat', 'is_place_hpc', 'is_object_hpc', 'is_place_ctx', 'is_object_ctx', ...
    'input_strength', 'feedback_strength', 'classification', 'best_speed_th', 'best_n_step');

fprintf('\n所有分析已完成！结果已保存到 haima/analysis_results.mat\n');
