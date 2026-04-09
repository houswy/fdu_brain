function output = prepdata(output_sort_final,  input_DLC, alignment)

    changed_a = alignment;
     
    % get binning
    input_DLC = input_DLC(changed_a(:,1),:);
    min_loc = min(input_DLC(:,2));
    max_loc = max(input_DLC(:,2));
    bin_num = 30;
    bin_width = (max_loc-min_loc)/bin_num;
    DLC_bin = floor((input_DLC(:,2)-min_loc)./bin_width)+1;
    start_bin = mode(DLC_bin(DLC_bin<5));
    start_loc = mean(input_DLC(DLC_bin==start_bin,2));
    end_bin = mode(DLC_bin(DLC_bin>25));
    end_loc = mean(input_DLC(DLC_bin==end_bin, 2));
    crop_DLC = input_DLC(:,2);
    crop_DLC(input_DLC(:,2)<start_loc) = start_loc;
    crop_DLC(input_DLC(:,2)>end_loc) = end_loc;
    start_loc = start_loc-1;
    end_loc = end_loc+1;
    bin_width_crop = (end_loc-start_loc)/bin_num;
    DLC_bin_crop = ceil((crop_DLC-start_loc)./bin_width_crop);
    
    % get forward running trials
    
    i = 1;
    forward_sta_temp=[];
    while i <=length(DLC_bin_crop)
        temp_idx= find(DLC_bin_crop(i:end)==3,1);
        temp_f_idx = find(DLC_bin_crop(i+temp_idx:end)==4,1);
        temp_r_idx = find(DLC_bin_crop(i+temp_idx:end)==2,1);
        if temp_r_idx<temp_f_idx
            i = i+temp_r_idx;
            continue
        end
        forward_sta_temp=[forward_sta_temp, i+temp_idx];
        i = i +temp_idx+temp_f_idx;
    end
    forward_lap_idx=[];
    for i = 1:length(forward_sta_temp)
        temp_idx = find(DLC_bin_crop(forward_sta_temp(i):end)==28,1);
        if isempty(temp_idx)
            continue
        end
        if (i<length(forward_sta_temp)) && (temp_idx>(forward_sta_temp(i+1)-forward_sta_temp(i)))
            continue
        end
        temp_DLC_bin = DLC_bin_crop(forward_sta_temp(i): forward_sta_temp(i)+temp_idx);
        rev_dir = diff(temp_DLC_bin)==-1;
        if sum(rev_dir)>0
            continue
        end
     forward_lap_idx=[forward_lap_idx;[forward_sta_temp(i), forward_sta_temp(i)+temp_idx]];
    end
    for i = 1:length(forward_lap_idx)
        forward_lap_frame(i,1)= changed_a(forward_lap_idx(i,1),2);
        forward_lap_frame(i,2)= changed_a(forward_lap_idx(i,2),2);
    end
    
    forward_lap=[];
    for i = 1:length(forward_lap_idx)
        forward_lap(i).DLC_raw = crop_DLC(forward_lap_idx(i,1):forward_lap_idx(i,2));
        forward_lap(i).DLC=DLC_bin_crop(forward_lap_idx(i,1):forward_lap_idx(i,2));
        forward_lap(i).MINI=changed_a(forward_lap_idx(i,1):forward_lap_idx(i,2),2);
        forward_lap(i).Calcium=output_sort_final.temporal_weights(forward_lap_frame(i,1): forward_lap_frame(i,2), :);
        forward_lap(i).Speed = mean(diff([forward_lap(i).DLC_raw]));
        for j = 3:28
            temp = [forward_lap(i).DLC]==j;
            temp_frame = forward_lap(i).MINI(temp);
            temp_calcium = output_sort_final.temporal_weights(temp_frame, :);
            forward_lap(i).Calcium_bin(:, j-2) = mean(temp_calcium,1);
        end
    end
    
    % get reverse running trials
    
    i = 1;
    reverse_sta_temp=[];
    while i <=length(DLC_bin_crop)
        temp_idx= find(DLC_bin_crop(i:end)==28,1);
        temp_f_idx = find(DLC_bin_crop(i+temp_idx:end)==27,1);
        temp_r_idx = find(DLC_bin_crop(i+temp_idx:end)==29,1);
        if temp_r_idx<temp_f_idx
            i = i+temp_r_idx;
            continue
        end
        reverse_sta_temp=[reverse_sta_temp, i+temp_idx];
        i = i +temp_idx+temp_f_idx;
    end
    reverse_lap_idx=[];
    for i = 1:length(reverse_sta_temp)
        temp_idx = find(DLC_bin_crop(reverse_sta_temp(i):end)==3,1);
        if isempty(temp_idx)
            continue
        end
        if (i<length(reverse_sta_temp)) && (temp_idx>(reverse_sta_temp(i+1)-reverse_sta_temp(i)))
            continue
        end
        temp_DLC_bin = DLC_bin_crop(reverse_sta_temp(i): reverse_sta_temp(i)+temp_idx);
        rev_dir = diff(temp_DLC_bin)==1;
        if sum(rev_dir)>0
            continue
        end
     reverse_lap_idx=[reverse_lap_idx;[reverse_sta_temp(i), reverse_sta_temp(i)+temp_idx]];
    end
    for i = 1:length(reverse_lap_idx)
        reverse_lap_frame(i,1)= changed_a(reverse_lap_idx(i,1),2);
        reverse_lap_frame(i,2)= changed_a(reverse_lap_idx(i,2),2);
    end
    
    reverse_lap=[];
    for i = 1:length(reverse_lap_idx)
        reverse_lap(i).DLC_raw = crop_DLC(reverse_lap_idx(i,1):reverse_lap_idx(i,2));
        reverse_lap(i).DLC=DLC_bin_crop(reverse_lap_idx(i,1):reverse_lap_idx(i,2));
        reverse_lap(i).MINI=changed_a(reverse_lap_idx(i,1):reverse_lap_idx(i,2),2);
        reverse_lap(i).Calcium=output_sort_final.temporal_weights(reverse_lap_frame(i,1): reverse_lap_frame(i,2), :);
        reverse_lap(i).Speed = mean(diff([reverse_lap(i).DLC_raw]));
        for j = 3:28
            temp = [reverse_lap(i).DLC]==j;
            temp_frame = reverse_lap(i).MINI(temp);
            temp_calcium = output_sort_final.temporal_weights(temp_frame, :);
            reverse_lap(i).Calcium_bin(:, j-2) = mean(temp_calcium,1);
        end
    end

    output.forward = forward_lap;
    output.reverse = reverse_lap;
    output.bin = DLC_bin_crop;
    output.align = changed_a;
    output.raw_calcium = output_sort_final.temporal_weights;
    output.totaltime = size(input_DLC,1);
    output.totallap = length(output.forward) + length(output.reverse);

end 

%%sort peaks    
%     for i = 1:length(forward_lap)
%         calcium_forward(:,:,i)=forward_lap(i).Calcium_bin;
%     end
%     calcium_forward_mean = mean(calcium_forward,3);
%     [~, peak_idx] = max(calcium_forward_mean, [], 2);
%     [~, sort_peak] = sort(peak_idx, "ascend");
%     
%     for i = 1:length(reverse_lap)
%         calcium_reverse(:,:,i)=reverse_lap(i).Calcium_bin;
%     end
%     calcium_reverse_mean = nanmean(calcium_reverse,3);
%     
%     % reward zone
%     
%     reward_end_idx=[];
%     for i = 1: length(forward_lap_idx)
%         temp_idx = find(reverse_lap_idx(:,1)>forward_lap_idx(i,2),1);
%         if isempty(temp_idx)
%             continue
%         end
%         temp_DLC = DLC_bin_crop(forward_lap_idx(i,2): reverse_lap_idx(temp_idx,1));
%         if sum(temp_DLC<20)>0
%             continue
%         end
%         reward_end_idx=[reward_end_idx; [forward_lap_idx(i,2),reverse_lap_idx(temp_idx,1)]];
%     end
%     
%     for i = 1:length(reward_end_idx)
%         reward_end_frame(i,1)= changed_a(reward_end_idx(i,1),2);
%         reward_end_frame(i,2)= changed_a(reward_end_idx(i,2),2);
%     end
%     
%     
%     reward_sta_idx=[];
%     for i = 1: length(reverse_lap_idx)
%         temp_idx = find(forward_lap_idx(:,1)>reverse_lap_idx(i,2),1);
%         if isempty(temp_idx)
%             continue
%         end
%         temp_DLC = DLC_bin_crop(reverse_lap_idx(i,2): forward_lap_idx(temp_idx,1));
%         if sum(temp_DLC>10)>0
%             continue
%         end
%         reward_sta_idx=[reward_sta_idx; [reverse_lap_idx(i,2),forward_lap_idx(temp_idx,1)]];
%     end
%     
%     for i = 1:length(reward_sta_idx)
%         reward_sta_frame(i,1)= changed_a(reward_sta_idx(i,1),2);
%         reward_sta_frame(i,2)= changed_a(reward_sta_idx(i,2),2);
%     end
% 
% %






