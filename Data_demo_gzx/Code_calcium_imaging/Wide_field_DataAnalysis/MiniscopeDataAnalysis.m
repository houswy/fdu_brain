%% Write avi to h5
V = VideoReader("0.avi");
%numframes = 2000;
  numFrames = 0;
     while hasFrame(V)
         readFrame(V);
         numFrames = numFrames + 1;
     end
h5create('0.h5','/images', [600,600,inf],'Datatype','uint8','Chunksize',[600,600,1000]);
for i = 1: numFrames 
    temp_img = read(V,i);
    h5write('0.h5','/images',temp_img(:,:,1),[1,1,i],[600,600,1]);
end  
%% Motion correction
input = 'combined.h5';
dataset_name='/image';
name_idx = strfind(input,'.h5');
output_name = strcat(input(1:name_idx-1),'_mc.h5');

info = h5info(input);
xdim = info.Datasets.Dataspace.Size(1);
ydim = info.Datasets.Dataspace.Size(2);
num_frame = info.Datasets.Dataspace.Size(3);
chunk = 1000;

temp = h5read(input,dataset_name,[1,1,1],[xdim,ydim,100]);
template = mean(temp,3);
chunk_num = ceil(num_frame/chunk);
h5create(output_name,dataset_name,[xdim,ydim,num_frame],"Chunksize",[xdim,ydim,1000],"Datatype","uint8");

for i = 1:chunk_num
    status = strcat('Processing frame:',num2str(i,2),'000');
    disp(status)
    frame_start = (i-1)*chunk+1;
    frame_end = min(i*chunk,num_frame);
    Yf = h5read(input,dataset_name,[1,1,frame_start],[xdim,ydim,frame_end-frame_start+1]);
    [d1,d2,T] = size(Yf);

    % high pass
    if (0)    
        hLarge = fspecial('average', 40);
        hSmall = fspecial('average', 2); 
        for t = 1:T
            Y(:,:,t) = filter2(hSmall,Yf(:,:,t)) - filter2(hLarge, Yf(:,:,t));
        end
        %Ypc = Yf - Y;
        bound = size(hLarge,1);
    else
        gSig = 7; 
        gSiz = 3*gSig; 
        psf = fspecial('gaussian', round(2*gSiz), gSig);
        ind_nonzero = (psf(:)>=max(psf(:,1)));
        psf = psf-mean(psf(ind_nonzero));
        psf(~ind_nonzero) = 0;   % only use pixels within the center disk
        %Y = imfilter(Yf,psf,'same');
        %bound = 2*ceil(gSiz/2);
        Y = imfilter(Yf,psf,'symmetric');
        bound = 0;
    end

options_r = NoRMCorreSetParms('d1',d1-bound,'d2',d2-bound,'bin_width',100,'max_shift',20,'iter',1,'correct_bidir',false,'upd_template',false);
tic; [~,shifts1,~] = normcorre_batch(Y(bound/2+1:end-bound/2,bound/2+1:end-bound/2,:),options_r); toc % register filtered data
tic; Mr = apply_shifts(Yf,shifts1,options_r,bound/2,bound/2); toc % apply shifts to full dataset

    h5write(output_name,dataset_name,Mr,[1,1,frame_start],[xdim,ydim,T]);
end
%% Down-sample with low RAM
input = 'combined_mc.h5';
dataset_name='/image';
name_idx = strfind(input,'.h5');
output_name = strcat(input(1:name_idx-1),'_ds.h5');

info = h5info(input);
xdim = info.Datasets.Dataspace.Size(1);
ydim = info.Datasets.Dataspace.Size(2);
num_frame = info.Datasets.Dataspace.Size(3);
chunk = 1000;
%downsample_factor = 3;
downsample_factor = 2;
temp = h5read(input,dataset_name,[1,1,1],[xdim,ydim,100]);
temp = single(temp);
temp_ds = downsample_space(temp,downsample_factor);
chunk_num = ceil(num_frame/chunk);
h5create(output_name,dataset_name,[size(temp_ds,1),size(temp_ds,2),num_frame],"Chunksize",[size(temp_ds,1),size(temp_ds,2),1000],"Datatype","single" );

for i = 1:chunk_num
    status = strcat('Processing frame:',num2str(i,2),'000');
    disp(status)
    frame_start = (i-1)*chunk+1;
    frame_end = min(i*chunk,num_frame);
    M = h5read(input,dataset_name,[1,1,frame_start],[xdim,ydim,frame_end-frame_start+1]);
    M = single(M);
    Md = downsample_space(M,downsample_factor);
    h5write(output_name,dataset_name,Md,[1,1,frame_start],[size(temp_ds,1),size(temp_ds,2),frame_end-frame_start+1]);
end
%%
% Extract
M = h5read('combined_mc_ds.h5','/image');
% M = h5read('0_mc_ds.h5','/image');
config = [];
config.use_gpu=0;
config = get_defaults(config);
config.avg_cell_radius=3;%改变extract识别的cell的半径大小
config.num_partitions_x=1;
config.num_partitions_y=1;
config.cellfind_min_snr = 2;
config.thresholds.T_min_snr=1;
config.spatial_highpass_cutoff =5;
config.max_iter=10;
%config.remove_duplicate_cells = 0; 
config.visualize_cellfinding=1;
output=extractor(M,config);
%% 
% Check cell, remove the bad cells
cell_check(output, M);
%% 
% Output_sort_final, final robust
output_sort1 = output;
[traces,filters]=final_robust_run(M,output_sort1,[output_sort1.user_labels]==1);
%[traces,filters]=final_robust_run(M,output_sort1,[output_sort1.user_labels]==-2);
%[traces,filters]=final_robust_run(M,output_sort1,([output_sort1.user_labels]==1)|([output_sort1.user_labels]==-2));
output_sort1_final=output_sort1;
output_sort1_final.spatial_weights=filters;
output_sort1_final.temporal_weights=traces;
%%
%plot figure, show the quality
figure, imshow(max(M,[],3),[0 350]);
plot_cells_overlay(output_sort1_final.spatial_weights,[1,0,0],[]);
plot_cells_overlay(output_sort1_final.spatial_weights(:,:,18),[1,0,0],[]);