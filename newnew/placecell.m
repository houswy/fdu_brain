function output = placecell(prepdata)
    %treat forward and reverse as different positions
    total_position = length(prepdata.bin);
    for i = 1:length(prepdata.forward)
        temp_position = prepdata.forward(i).DLC;
        for j = 3:28
            loc_count = sum(temp_position==j);
            position_occupy(j-2, i) = loc_count;
        end
    end
    p_forward = sum(position_occupy, 2)./total_position;
    position_occupy = [];
    for i = 1:length(prepdata.reverse)
        temp_position = prepdata.reverse(i).DLC;
        for j = 3:28
            loc_count = sum(temp_position==j);
            position_occupy(j-2, i) = loc_count;
        end
    end
    p_reverse = sum(position_occupy, 2)./total_position;
    r_mean = mean(prepdata.raw_calcium(prepdata.align(:,2), :), 1);
   frame_bin_F={};
        for j = 3:28
            bin_frame=[];
             for i = 1:length(prepdata.forward)
                 temp = prepdata.forward(i).DLC==j;
                 temp_frame = unique(prepdata.forward(i).MINI(temp));
                 bin_frame = [bin_frame; temp_frame];
             end
            frame_bin_F{j-2}=bin_frame;
        end

    for i = 1:length(frame_bin_F)
        temp_frame = frame_bin_F{i};
        temp_calcium = mean(prepdata.raw_calcium(temp_frame, :), 1);
        calcium_bin_F(i, :) = temp_calcium;
    end

       frame_bin_R={};
        for j = 3:28
            bin_frame=[];
             for i = 1:length(prepdata.reverse)
                 temp = prepdata.reverse(i).DLC==j;
                 temp_frame = unique(prepdata.reverse(i).MINI(temp));
                 bin_frame = [bin_frame; temp_frame];
             end
            frame_bin_R{j-2}=bin_frame;
        end
    for i = 1:length(frame_bin_R)
        temp_frame = frame_bin_R{i};
        temp_calcium = mean(prepdata.raw_calcium(temp_frame, :), 1);
        calcium_bin_R(i, :) = temp_calcium;
    end

    calcium_bin = [calcium_bin_F; calcium_bin_R];
    p = [p_forward; p_reverse];

    cell_num = length(r_mean);
    for i = 1: cell_num
        temp_SI = [];
        for j = 1:length(p)
            temp_SI (j) = p(j)*calcium_bin(j, i)/r_mean(i)*log2(calcium_bin(j, i)/r_mean(i));
        end
        SI(i) = nansum(temp_SI);
        temp_SI = [];
        for j = 1:length(p_forward)
            temp_SI (j) = p_forward(j)*calcium_bin_F(j, i)/r_mean(i)*log2(calcium_bin_F(j, i)/r_mean(i));
        end
        SI_F(i) = nansum(temp_SI);
        temp_SI = [];
        for j = 1:length(p_reverse)
            temp_SI (j) = p_reverse(j)*calcium_bin_R(j, i)/r_mean(i)*log2(calcium_bin_R(j, i)/r_mean(i));
        end
        SI_R(i) = nansum(temp_SI);        
    end
    
    total_frame = size(prepdata.raw_calcium,1);
    for m = 1:500
        shift = round(rand*size(prepdata.raw_calcium,1)/2);
        for i = 1:length(frame_bin_F)
            temp_frame = frame_bin_F{i}+shift;
            temp_frame(temp_frame>total_frame)= temp_frame(temp_frame>total_frame) - total_frame;
            temp_calcium = mean(prepdata.raw_calcium(temp_frame, :), 1);
            calcium_bin_F_shuffle(i, :) = temp_calcium;
        end
        for i = 1:length(frame_bin_R)
            temp_frame = frame_bin_R{i}+shift;
            temp_frame(temp_frame>total_frame)= temp_frame(temp_frame>total_frame) - total_frame;
            temp_calcium = mean(prepdata.raw_calcium(temp_frame, :), 1);
            calcium_bin_R_shuffle(i, :) = temp_calcium;
        end
    
        calcium_bin_shuffle = [calcium_bin_F_shuffle; calcium_bin_R_shuffle];

        for i = 1: cell_num
        temp_SI = [];
            for j = 1:length(p)
                temp_SI (j) = p(j)*calcium_bin_shuffle(j, i)/r_mean(i)*log2(calcium_bin_shuffle(j, i)/r_mean(i));
            end
        SI_shuffle(m,i) = nansum(temp_SI);

         temp_SI = [];
            for j = 1:length(p_forward)
                temp_SI (j) = p_forward(j)*calcium_bin_F_shuffle(j, i)/r_mean(i)*log2(calcium_bin_F_shuffle(j, i)/r_mean(i));
            end
        SI_F_shuffle(m,i) = nansum(temp_SI);       

           temp_SI = [];
            for j = 1:length(p_reverse)
                temp_SI (j) = p_reverse(j)*calcium_bin_R_shuffle(j, i)/r_mean(i)*log2(calcium_bin_R_shuffle(j, i)/r_mean(i));
            end
        SI_R_shuffle(m,i) = nansum(temp_SI);       
        end      
    end

        sigma = std(SI_shuffle(:,i));
        mu = mean(SI_shuffle(:, i));
        sigma_F=std(SI_F_shuffle(:,i));
        mu_F = mean(SI_F_shuffle(:, i));
        sigma_R=std(SI_R_shuffle(:,i));
        mu_R = mean(SI_R_shuffle(:, i));

    for i = 1:cell_num
        [~, test(i)] = ztest(SI(i), mu, sigma, "Tail","right");
        [~, test_F(i)] = ztest(SI_F(i), mu_F, sigma_F, "Tail","right");
        [~, test_R(i)] = ztest(SI_R(i), mu_R, sigma_R, "Tail","right");
    end

    place_cell_Com = find(test<0.01);
    place_cell_F = find(test_F<0.01);
    place_cell_R = find(test_R<0.01);
    list = unique([place_cell_Com, place_cell_F, place_cell_R]);
    
    for  i = 1:length(prepdata.forward)
        forward_calcium_bin(:,:,i) = prepdata.forward(i).Calcium_bin(list, :);
    end
    forward_calcium_bin_mean = mean(forward_calcium_bin, 3);
    [~, peak_loc] = max(forward_calcium_bin_mean, [], 2);
    [~, peak_sort_F] = sort(peak_loc, "ascend");
    for  i = 1:length(prepdata.reverse)
        reverse_calcium_bin(:,:,i) = prepdata.reverse(i).Calcium_bin(list, :);
    end
    reverse_calcium_bin_mean = mean(reverse_calcium_bin, 3);
    [~, peak_loc] = max(reverse_calcium_bin_mean, [], 2);
    [~, peak_sort_R] = sort(peak_loc, "ascend");

    % calculate norm response
    for i = 1:length(prepdata.forward)
        baseline = min(prepdata.forward(i).Calcium(:, list), [],1)';
        peak = max(prepdata.forward(i).Calcium(:, list), [],1)';
        forward_calcium_bin_norm(:,:,i) = ((forward_calcium_bin(:,:,i)-baseline)./(peak-baseline));
    end

    forward_calcium_bin_template = nanmean(forward_calcium_bin_norm, 3);
    for i = 1:length(prepdata.reverse)
        baseline = min(prepdata.reverse(i).Calcium(:, list), [],1)';
        peak = max(prepdata.reverse(i).Calcium(:, list), [],1)';
        reverse_calcium_bin_norm(:,:,i) = ((reverse_calcium_bin(:,:,i)-baseline)./(peak-baseline));
    end

    reverse_calcium_bin_template = nanmean(reverse_calcium_bin_norm, 3);
    
    % correlation to the template response
    for i = 1: length(list)
        lap_corr= [];
        for j = 1:length(prepdata.forward)
            temp_template = forward_calcium_bin_template(i,:);
            temp_target = prepdata.forward(j).Calcium_bin(list(i),:);
            temp_corr = corrcoef(temp_template', temp_target',"Rows","complete");
            lap_corr=[lap_corr, temp_corr(1,2)];
        end
        cell_corr{i,1}=lap_corr;
        lap_corr= [];
        for j = 1:length(prepdata.reverse)
            temp_template = reverse_calcium_bin_template(i,:);
            temp_target = prepdata.reverse(j).Calcium_bin(list(i),:);
            temp_corr = corrcoef(temp_template', temp_target',"Rows","complete");
            lap_corr=[lap_corr, temp_corr(1,2)];
        end
        cell_corr{i,2}=lap_corr;
    end

    for i = 1:length(list)
        [~,p_corr(i)]=ttest2(cell_corr{i,1}, cell_corr{i,2});
    end

    % response size
     for i = 1: length(list)
        lap_response= [];
        for j = 1:length(prepdata.forward)
            temp_target = prepdata.forward(j).Calcium(:,list(i));
            temp_response = max(temp_target)-min(temp_target);
            lap_response=[lap_response, temp_response];
        end
        cell_response{i,1}=lap_response;
        lap_response= [];
        for j = 1:length(prepdata.reverse)
            temp_target = prepdata.reverse(j).Calcium(:,list(i));
            temp_response = max(temp_target)-min(temp_target);
            lap_response=[lap_response, temp_response];
        end
        cell_response{i,2}=lap_response;
     end

     for i = 1:length(list)
        [~,p_rsp(i)]=ttest2(cell_response{i,1}, cell_response{i,2});
     end

    place_cell_N = setdiff(place_cell_Com,place_cell_R);
    place_cell_NN = setdiff(place_cell_N,place_cell_F);
    SI_RR = SI_R(place_cell_R);
    SI_FF = SI_F(place_cell_F);
    SI_N = SI(place_cell_NN);
    SI_U = union(SI_FF,SI_RR);
    SI_UN = union(SI_U,SI_N)';

     % output
     for i=1:length(list)
         output(i).Cell_Number = list(i);
         if p_rsp(i)<0.05
             if mean(cell_response{i,1})>mean(cell_response{i,2})
                 output(i).DIR = 'F';
             else
                 output(i).DIR = 'R';
             end
         elseif p_corr(i) <0.05
             if mean(cell_corr{i,1})>mean(cell_corr{i,2})
                 output(i).DIR = 'F';
             else
                 output(i).DIR = 'R';
             end
         else
             output(i).DIR = 'N';
         end
         if output(i).DIR == 'F'
             temp_norm = squeeze(forward_calcium_bin_norm(i, :,:));
             temp_ave = forward_calcium_bin_template(i,:);
             temp_corr = nanmean(cell_corr{i,1});
         elseif output(i).DIR == 'R'
             temp_norm = squeeze(reverse_calcium_bin_norm(i, :,:));
             temp_ave = reverse_calcium_bin_template(i,:);
             temp_corr = nanmean(cell_corr{i,2});
         elseif output(i).DIR =='N'
             temp_norm = [squeeze(forward_calcium_bin_norm(i, :,:)), squeeze(reverse_calcium_bin_norm(i, :,:))];
             temp_ave = (forward_calcium_bin_template(i,:)+reverse_calcium_bin_template(i,:))./2;
             temp_corr = (nanmean(cell_corr{i,1})+nanmean(cell_corr{i,2}))/2;
         end
         output(i).Lap_norm = temp_norm';
         output(i).Ave_norm = temp_ave;
         output(i).Correlation = temp_corr;
     end

     for i = 1:length(list)
         [peak, ~] = max(output(i).Ave_norm);
         baseline = min(output(i).Ave_norm);
         [~, loc_r, w, ~] = findpeaks([baseline, baseline, output(i).Ave_norm,baseline,baseline], "MinPeakHeight", peak*0.95, "MinPeakDistance", 10;
         output(i).Peak = loc_r;
         output(i).Width = w;
     end

%     for i = 1:length(list)
%         if  output(i).DIR == 'F'
%             temp_cell_number = find(place_cell_F,i)
%             temp_SI_F = SI_F(:,temp_cell_number);
%             output(i).SI = temp_SI_F;
%         elseif  output(i).DIR == 'R'
%             temp_SI_R = SI_R(:,i);
%             output(i).SI = temp_SI_R;
%         elseif  output(i).DIR == 'N'
%             temp_SI_N = SI(:,i);
%             output(i).SI = temp_SI_N;
%         end
%     end


%     place_cell_SI = SI(list);
%     place_cell_SI = place_cell_SI'
%     place_cell_SI = mean(place_cell_SI)
%     P_value = 2*normcdf(place_cell_SI');
%     place_cells_P_value_SI = mean(P_value);

end    