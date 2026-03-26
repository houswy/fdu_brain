# %%
import numpy as np
from rnn_class import RNN
import matplotlib.pyplot as plt
from numpy.linalg import inv
from utils_admm import solve_corrn_admm_gpu
from utils_admm import solve_corrn_admm
from sklearn.decomposition import PCA

from utils_admm import solve_corrn_admm_lr_gpu

# %%
from scipy.io import loadmat
df=loadmat('data/C597.mat')

# %%
firing_rate = df['firing_rate']
opts = {};
opts['g'] = 1.8;
opts['n_rec'] = firing_rate.shape[1]
opts['sigma_input'] = 0
opts['sigma_conversion'] = 0
opts['alpha'] = 0.1
opts['num_cores'] = 8
r_in=firing_rate[:-1,:]
r_tar=firing_rate[1:,:]

# %%
u_in=df['input'][:-1,:]
opts['bin_max'] = u_in.max()
u_in_times = np.zeros((u_in.shape[1],2))
u_in_times[[0,1,u_in.shape[1]-2,u_in.shape[1]-1],0] = opts['bin_max']
u_in_times[2:u_in.shape[1]-3,1] = opts['bin_max']
u_in = u_in @ u_in_times
opts['n_in'] = u_in.shape[1]
temp_input = df['input'][:-1,:]

# %%
w = solve_corrn_admm_lr_gpu(r_in,r_tar,u_in = u_in, alph = opts['alpha'], K_rank = 100,
                l2 = 1e-4, threshold = 1, rho = 10000,
                verbose = 2,mask = None,
                num_iters = 20)

# %%
theta = w.T;
w_rec = theta[:,:-1*u_in.shape[1]];
w_in = theta[:,-1*u_in.shape[1]:];
m1 = RNN(opts)
m1.rnn['w_rec'] = w_rec
m1.rnn['w_in']  = w_in

# %%
plt.imshow(theta,cmap = 'jet',interpolation = 'none'); 
plt.colorbar()
plt.show()

# %%
# one step predict
n=1
r_pred=np.zeros(r_in.shape)
for i in range(r_in.shape[0]):
    temp_r=m1.get_time_evolution(r_in[i,:],T=n,u=u_in[i:i+n,:])
    r_pred[i,:]=temp_r[-1,:]

# %%
start_idx=df['start']
end_idx=df['stop']

# %%
# forward/reverse input to get binned firing rate
r_pred_bin=np.zeros([temp_input.shape[1],r_in.shape[1],start_idx.shape[0]])
r_GT_bin=np.zeros([temp_input.shape[1],r_in.shape[1],start_idx.shape[0]])
for i in range(start_idx.shape[0]):
    
    temp_pred_bin=r_pred[start_idx[i][0]:end_idx[i][0],:]
    temp_GT_bin=r_in[start_idx[i][0]:end_idx[i][0],:]

    bin_times=temp_input[start_idx[i][0]:end_idx[i][0],:].astype(int).T
    denom=np.sum(bin_times,axis=1)[:,np.newaxis]
    denom[denom==0]=1

    temp_bin=bin_times@temp_pred_bin
    temp_bin=(temp_bin/denom)
    r_pred_bin[:,:,i]=temp_bin

    temp_bin=bin_times@temp_GT_bin
    temp_GT_bin=(temp_bin/denom)
    r_GT_bin[:,:,i]=temp_GT_bin

# %%
r_pred_mean=np.nanmean(r_pred_bin,axis=2)
r_GT_mean=np.nanmean(r_GT_bin,axis=2)
gt_pred_rmse=np.sqrt(np.mean((r_GT_mean-r_pred_mean)**2,axis=1))
gt_pred_rmse=np.mean(gt_pred_rmse)
print('GT pred RMSE:',gt_pred_rmse)

# %%
import matplotlib.pyplot as plt
fig = plt.figure()
ax = fig.add_subplot(121)
sorted_bin=r_GT_mean[:,np.argsort(np.argmax(r_GT_mean,axis=0))]
im1=ax.imshow(sorted_bin.T,cmap = 'jet',interpolation = 'none',aspect='auto');
ax.set_title('GT')
ax.set_xlabel('bin')
ax.set_ylabel('Neuron index')
plt.colorbar(im1,ax=ax,label='FR')

ax = fig.add_subplot(122)
sorted_bin=r_pred_mean[:,np.argsort(np.argmax(r_GT_mean,axis=0))]
im2=ax.imshow(sorted_bin.T,cmap = 'jet',interpolation = 'none',aspect='auto');
if n!='multi':
    ax.set_title('train %d input CoRNN %d step RMSE %f' % (opts['n_in'],n,gt_pred_rmse))
else:
    ax.set_title('train %d input CoRNN multi step RMSE %f' % (opts['n_in'],gt_pred_rmse))
ax.set_xlabel('bin')
ax.set_ylabel('Neuron index')
plt.colorbar(im2,ax=ax,label='FR')

plt.tight_layout()
plt.show()
