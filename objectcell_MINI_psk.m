function object_cell_idx = objectcell_MINI_psk(behavior, imaging, object)

cellActivity = double(imaging.temporal_weights');
%function object_assembly = Objectassembly_psk(celllist, behavior, imaging, plot_flag)
a=randperm(size(imaging.temporal_weights,2));
celllist=sort(a);
trial_num = find([behavior.trial_type]==2 | [behavior.trial_type]==1);
%trial_num = find([behavior.trial_type]==2);
trial_num_trash = find([behavior.trial_type]==0);
%trial_num_trash = find([behavior.trial_type]==0 | [behavior.trial_type]==1);

K=object;
% K=K(trial_num,:);

behavior_use=behavior;
behavior_use([trial_num_trash]) = [];

object_frame = [];
for i = 1:length(trial_num)
    idx = find(ismember(behavior_use(i).location_bin, [K(i)-1, K(i), K(i)+1]));
    bins = behavior_use(i).MINI_sto-behavior_use(i).MINI_sta+1;
    H = [];
    H = zeros(bins,1)';
    trial_sta = behavior_use(i).MINI_frames(idx(1)) - behavior_use(i).MINI_sta + 1;
    trial_sto = behavior_use(i).MINI_frames(idx(end)) - behavior_use(i).MINI_sta + 1;
    H(trial_sta:trial_sto) = 1;
    object_frame = [object_frame, H];
end
x = object_frame;

cell_activity = [];
cell_activity_trial = struct();
cell_binary = struct();
y = struct();
% for i = 1:length(trial_num)
%     fieldName = sprintf('Trial_%d', i);  % 生成字段名，如"field1", "field2"
%     assembly_activity_trial.(fieldName) = [];         % 初始化字段值为空数组（可替换为其他默认值）
% end
for j = 1:size(cellActivity,1)
    row_data = cellActivity(j, :)
    current_row = [];
    current_row_norm = [];
    for i = 1:length(trial_num)
        G = [];
        G = row_data(behavior_use(i).MINI_sta:behavior_use(i).MINI_sto)
        J = [];
        J = rescale(G);
        J(J>=0.8) = 1;
        J(J<0.8) = 0;
        J(isnan(J)) = 0;
        fieldName = sprintf('Trial_%d', i);  % 生成字段名，如"field1", "field2"
        cell_activity_trial(j).(fieldName) = G;
        cell_binary(j).(fieldName) = J;
        current_row = [current_row, G];
        current_row_norm = [current_row_norm, J];
    end
    cell_activity(j,:) = current_row;
    y(j).activity = current_row_norm;
end

%calculate MI
for i=1:length(y)
    MI(i) = mutualinfo(x, y(i).activity);
end

%shuffle
nPermutations = 1000; % 置换次数
%MI_shuffle = zeros(nPermutations, 1);
for i = 1:nPermutations
    s = rand;
    s_shift = round(s*length(x)*0.3)+1;
    x_permuted = [x(s_shift:length(x)), x(1:s_shift-1)]; % 随机置换x
    for j=1:length(MI)
        MIs(j,i) = mutualinfo(x_permuted, y(j).activity);
    end
end
MIs=MIs'

for i=1:length(MI)
    sigma = std(MIs(:,i));
    mu = median(MIs(:,i));
    [~, test(i)] = ztest(MI(i), mu, sigma, "Tail","right");
end

object_cell_idx = find(test<0.05);

end