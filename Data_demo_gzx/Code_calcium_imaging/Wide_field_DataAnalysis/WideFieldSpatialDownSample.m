input = 'combined_mc.h5';
%input = 'combined.h5';
info = h5info(input);
xdim = info.Datasets.Dataspace.Size(1);
ydim = info.Datasets.Dataspace.Size(2);
setname = ['/',info.Datasets.Name];

frame = info.Datasets.Dataspace.Size(3);
M=[];
for i = 1:frame
    disp(i)
    temp = h5read(input, setname,[1,1,i],[xdim,ydim,1]);
    %temp_ds = downsample_space(temp, 4); %修改这个数字可以调整缩放抽样比例
    temp_ds = downsample_space(temp, 3); %修改这个数字可以调整缩放抽样比例
    M(:,:,i) = temp_ds;
end
M = single(M);

xdim_ds = size(M,1);
ydim_ds = size(M,2);
idx = strfind(input,'.h5');
output = [input(1:idx-1),'_ds.h5'];
h5create(output,setname,[xdim_ds,ydim_ds,frame],"Datatype","single");
h5write(output,setname,M);





data = h5read('combined_mc.h5', '/datasetname');
numFrames = size(data, 3); % 假设数据的第一维代表帧
data(1:2:numFrames, :) = NaN; % 将奇数帧设置为NaN
evenFrames = data(2:2:end, :); % 读取偶数帧的数据
h5write('combined_mc_fixed.h5', '/newdatasetname', evenFrames);

input = 'combined_mc.h5';
info = h5info(input);
setname = ['/',info.Datasets.Name];
frame = info.Datasets.Dataspace.Size(3);
info(1:2:frame, :) = NaN; % 将奇数帧设置为NaN
evenFrames = info(2:2:end, :); % 读取偶数帧的数据
h5write('combined_mc_fixed.h5', setname, evenFrames);




input = 'combined_mc.h5';
info = h5info(input);
xdim = info.Datasets.Dataspace.Size(1);
ydim = info.Datasets.Dataspace.Size(2);
setname = ['/',info.Datasets.Name];

frame = info.Datasets.Dataspace.Size(3);
M=[];
for i = 2:2:frame
    disp(i)
    temp = h5read(input, setname,[1,1,i],[xdim,ydim,1]);
    temp_ds = downsample_space(temp, 2); %修改这个数字可以调整缩放抽样比例
    M(:,:,i/2) = temp_ds;
end
M = single(M);

xdim_ds = size(M,1);
ydim_ds = size(M,2);
idx = strfind(input,'.h5');
output = [input(1:idx-1),'_ds.h5'];
h5create(output,setname,[xdim_ds,ydim_ds,frame/2],"Datatype","single");
h5write(output,setname,M);
