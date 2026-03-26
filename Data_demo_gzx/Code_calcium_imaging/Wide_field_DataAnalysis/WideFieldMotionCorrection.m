input = 'combined.h5';
dataset_name='/image';
name_idx = strfind(input,'.h5');
output_name = strcat(input(1:name_idx-1),'_mc.h5');

info = h5info(input);
xdim = info.Datasets.Dataspace.Size(1);
ydim = info.Datasets.Dataspace.Size(2);
num_frame = info.Datasets.Dataspace.Size(3);
chunk = 400;

temp = h5read(input,dataset_name,[1,1,1],[xdim,ydim,100]);
template = mean(temp,3);
chunk_num = ceil(num_frame/chunk);
h5create(output_name,dataset_name,[xdim,ydim,num_frame],"Chunksize",[xdim,ydim,400],"Datatype","single");

for i = 1:chunk_num
    status = strcat('Processing frame:',num2str(i,2),'000');
    disp(status)
    frame_start = (i-1)*chunk+1;
    frame_end = min(i*chunk,num_frame);
    Yf = h5read(input,dataset_name,[1,1,frame_start],[xdim,ydim,frame_end-frame_start+1]);
    [d1,d2,T] = size(Yf);
    Yf = single(Yf);
    % high pass
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

options_r = NoRMCorreSetParms('d1',d1-bound,'d2',d2-bound,'bin_width',100,'max_shift',20,'iter',1,'correct_bidir',false,'upd_template',false);
tic; [~,shifts1,~] = normcorre_batch(Y(bound/2+1:end-bound/2,bound/2+1:end-bound/2,:),options_r,template); toc % register filtered data
tic; Mr = apply_shifts(Yf,shifts1,options_r,bound/2,bound/2); toc % apply shifts to full dataset

    h5write(output_name,dataset_name,Mr,[1,1,frame_start],[xdim,ydim,T]);
end