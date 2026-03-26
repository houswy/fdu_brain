list=dir('F:\906\WF\906-20260325-mode1+mode2');

use_list=[];
for i = 1:length(list)
    temp_name = list(i).name;
    if contains(temp_name, '.tif')
        use_list=[use_list, i];
    end
end

for i=1:length(use_list)
    temp_name=list(use_list(i)).name;
    sta_idx = strfind(temp_name,'_t');
    sto_idx = strfind(temp_name,'.tif');
    temp_num = temp_name(sta_idx+2:sto_idx-1);
    list_img_num(i) = str2num(temp_num);
end

[~, sort_img_list]=sort(list_img_num,'ascend');

temp_img = imread(list(use_list(sort_img_list(1))).name);
num_img = length(use_list);

% temp_img_spatial_ds = imresize(temp_img,0.5);
 xdim = size(temp_img,1);
 ydim = size(temp_img,2);

h5create('combined.h5','/image',[xdim, ydim, num_img],'Datatype','uint16');
for  i = 1:num_img
    disp(i)
    
    temp_img = imread(list(use_list(sort_img_list(i))).name);
    %temp_img_spatial_ds = imresize(temp_img,0.5);
    h5write('combined.h5','/image', temp_img,[1,1,i],[xdim, ydim,1]);
end
    