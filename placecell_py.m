function place_cell_idx = placecell_py(behavior, imaging)
%imaging.temporal_weights = [];
%imaging.temporal_weights = assemblyActivity';
%place_assembly_idx = placecell_py(behavior, imaging)

    num_frame = size(imaging.temporal_weights, 1);
    total_trial = length(behavior);
    if abs(num_frame - behavior(total_trial).WF_sto)<abs(num_frame - behavior(total_trial).MINI_sto)
        frame_sta = [behavior.WF_sta];
        frame_sto = [behavior.WF_sto];
        choice = 1;
    else
        frame_sta = [behavior.MINI_sta];
        frame_sto = [behavior.MINI_sto];
        choice = 2;
    end
    
    % add threshold for speed
    total_bin=[];
    for i = 1:total_trial
        if behavior(i).trial_type == 2
            temp_bin=behavior(i).location_bin;
            total_bin=[total_bin;temp_bin];
        end
    end
    num_bin = unique(total_bin);
    use_bin = num_bin(3:end-2);

    for i = 1:total_trial
        if behavior(i).trial_type == 2
            temp_bin=behavior(i).location_bin;
            remove_idx = [find(temp_bin==0); find(temp_bin==1); find(temp_bin==max(num_bin)); find(temp_bin==max(num_bin)-1)];
            temp_bin_crop = temp_bin;
            temp_bin_crop(remove_idx)=[];
            temp_bin_crop_diff = diff(temp_bin_crop);
            temp_bin_crop_idx_diff = diff(find(temp_bin_crop_diff==1));
            if sum(temp_bin_crop_idx_diff>300)>0
                behavior(i).trial_type = 3;
            end
        end
    end

    total_bin=[];
    for i = 1:total_trial
        if behavior(i).trial_type == 2
            temp_bin=behavior(i).location_bin;
            total_bin=[total_bin;temp_bin];
        end
    end
    num_bin = unique(total_bin);
    use_bin = num_bin(3:end-2);

    remove_idx = [find(total_bin==0); find(total_bin==1); find(total_bin==max(num_bin)); find(total_bin==max(num_bin)-1)];
    total_bin_crop = total_bin;
    total_bin_crop(remove_idx)=[];

    for i = 1:length(use_bin)
        p(i) = sum(total_bin_crop==use_bin(i))/length(total_bin_crop);
    end
    

    
    
    disp('calculating spatial infomation...')
    correct_trial = find([behavior.trial_type]==2);
    calcium_bin_raw = zeros([length(use_bin),size(imaging.temporal_weights, 2), length(correct_trial)]);
    for i = 1:length(correct_trial)
        if choice==1
            frames = behavior(correct_trial(i)).WF_frames;
        else
            frames = behavior(correct_trial(i)).MINI_frames;
        end
        temp_bin = behavior(correct_trial(i)).location_bin;
        temp_activity = imaging.temporal_weights(frames, :);
        for j = 1:length(use_bin)
            temp_idx = temp_bin==use_bin(j);
            if isempty(temp_idx)
                bin_activity = nan(1,size(imaging.temporal_weights, 2));
            else    
                sub_activity = temp_activity(temp_idx,:);
                bin_activity = mean(sub_activity, 1);
            end
            calcium_bin_raw(j,:,i) = bin_activity;
        end
    end
    for i = 1:size(calcium_bin_raw,2)
        for j = 1:size(calcium_bin_raw,3)
            temp = calcium_bin_raw(:,i,j);
            norm_temp = rescale(temp);
            calcium_bin_norm(:,i,j)=norm_temp;
        end
    end
    calcium_bin=nanmean(calcium_bin_norm,3);

%     for i = 1:size(calcium_bin_raw,2)
%         temp = squeeze(calcium_bin_raw(:,i,:));
%         r_mean(i) = nanmean(temp,'all');
%     end

    for i = 1:size(calcium_bin_norm,2)
        temp = squeeze(calcium_bin_norm(:,i,:));
        r_mean(i) = nanmean(temp,'all');
    end

   % r_mean_2=mean(calcium_bin,1);
    
    for i = 1:size(imaging.temporal_weights, 2)
        temp_SI = [];
        for j = 1: length(use_bin)
            temp_SI(j) = p(j)*calcium_bin(j, i)/r_mean(i)*log2(calcium_bin(j, i)/r_mean(i));
        end
        SI(i)=nansum(temp_SI);
    end
       
    % add random shuffle
    disp('shuffle test...')
    for s = 1:100
        m = rng;
        shuffle = round(rand(1)*num_frame/2);
        shuffle_activity = [imaging.temporal_weights(shuffle:num_frame,:);imaging.temporal_weights(1:shuffle-1,:)];
        %shuffle_activity = [imaging.temporal_weights(shuffle:num_frame,:);imaging.temporal_weights(1:shuffle,:)];
        calcium_bin_raw_shuffle = zeros([length(use_bin),size(imaging.temporal_weights, 2), length(correct_trial)]);   
        for i = 1:length(correct_trial)
            if choice==1
                frames = behavior(correct_trial(i)).WF_frames;
            else
                frames = behavior(correct_trial(i)).MINI_frames;
            end
            temp_bin = behavior(correct_trial(i)).location_bin;
            temp_activity = shuffle_activity(frames, :);
            for j = 1:length(use_bin)
                temp_idx = temp_bin==use_bin(j);
                if isempty(temp_idx)
                    bin_activity = nan(1,size(imaging.temporal_weights, 2));
                else    
                    sub_activity = temp_activity(temp_idx,:);
                    bin_activity = mean(sub_activity, 1);
                end
                calcium_bin_raw_shuffle(j,:,i) = bin_activity;
            end
        end
        for i = 1:size(calcium_bin_raw_shuffle,2)
            for j = 1:size(calcium_bin_raw_shuffle,3)
                temp = calcium_bin_raw_shuffle(:,i,j);
                norm_temp = rescale(temp);
                calcium_bin_norm_shuffle(:,i,j)=norm_temp;
            end
        end
        calcium_bin_shuffle=nanmean(calcium_bin_norm_shuffle,3);
    
%         for i = 1:size(calcium_bin_raw_shuffle,2)
%             temp = squeeze(calcium_bin_raw_shuffle(:,i,:));
%             r_mean_shuffle(i) = nanmean(temp,'all');
%         end

        for i = 1:size(calcium_bin_norm_shuffle,2)
             temp = squeeze(calcium_bin_norm_shuffle(:,i,:));
             r_mean_shuffle(i) = nanmean(temp,'all');
         end
        
        for i = 1:size(imaging.temporal_weights, 2)
            temp_SI = [];
            for j = 1: length(use_bin)
                temp_SI(j) = p(j)*calcium_bin_shuffle(j, i)/r_mean_shuffle(i)*log2(calcium_bin_shuffle(j, i)/r_mean_shuffle(i));
            end
            SI_shuffle(i,s)=nansum(temp_SI);
        end
    end

    for i = 1:size(imaging.temporal_weights, 2)
        sigma = std(SI_shuffle(i,:));
        mu = median(SI_shuffle(i,:));
        [~, test(i)] = ztest(SI(i), mu, sigma, "Tail","right");
    end
    
    place_cell_idx = find(test<0.05);
   % place_cell_idx = find(test<0.05/size(imaging.temporal_weights, 2));
end


    %% place cell description
%     place_cell_output = struct();
%     for i = 1:length(place_cell_idx)
%         temp_calcium_bin = squeeze(calcium_bin_raw(:,place_cell_idx(i),:));
%         temp_calcium_bin_ave = nanmean(temp_calcium_bin,2);
%         plot(temp_calcium_bin_ave);
%         hold on
%     end