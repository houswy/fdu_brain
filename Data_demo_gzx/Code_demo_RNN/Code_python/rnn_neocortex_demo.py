# %% 导入必要的库
import numpy as np
from rnn_class import RNN
import matplotlib.pyplot as plt
from numpy.linalg import inv
from utils_admm import solve_corrn_admm_lr_gpu
from sklearn.decomposition import PCA
from scipy.io import loadmat

# %% 加载新皮层数据
df = loadmat('data/C597.mat')

# %% 数据预处理
firing_rate = df['firing_rate']
opts = {}
opts['g'] = 1.8
opts['n_rec'] = firing_rate.shape[1]
opts['sigma_input'] = 0
opts['sigma_conversion'] = 0
opts['alpha'] = 0.1
opts['num_cores'] = 8
r_in = firing_rate[:-1,:]
r_tar = firing_rate[1:,:]

# %% 处理位置输入
u_in = df['input'][:-1,:]
opts['bin_max'] = u_in.max()
u_in_times = np.zeros((u_in.shape[1],2))
u_in_times[[0,1,u_in.shape[1]-2,u_in.shape[1]-1],0] = opts['bin_max']
u_in_times[2:u_in.shape[1]-3,1] = opts['bin_max']
u_in = u_in @ u_in_times
opts['n_in'] = u_in.shape[1]
temp_input = df['input'][:-1,:]

# %% 使用ADMM算法求解RNN权重
print('训练新皮层RNN模型...')
w = solve_corrn_admm_lr_gpu(r_in, r_tar, u_in=u_in, alph=opts['alpha'], K_rank=100,
                l2=1e-4, threshold=1, rho=10000,
                verbose=2, mask=None,
                num_iters=20)

# %% 创建RNN模型并设置权重
theta = w.T
w_rec = theta[:,:-1*u_in.shape[1]]
w_in = theta[:,-1*u_in.shape[1]:]
m1 = RNN(opts)
m1.rnn['w_rec'] = w_rec
m1.rnn['w_in'] = w_in

# %% 可视化权重矩阵
plt.figure(figsize=(10, 8))
plt.imshow(theta, cmap='jet', interpolation='none')
plt.colorbar(label='权重值')
plt.title('新皮层RNN权重矩阵')
plt.xlabel('输入维度')
plt.ylabel('神经元索引')
plt.tight_layout()
plt.savefig('neocortex_weights.png')
plt.show()

# %% 进行一步预测
n = 1
r_pred = np.zeros(r_in.shape)
for i in range(r_in.shape[0]):
    temp_r = m1.get_time_evolution(r_in[i,:], T=n, u=u_in[i:i+n,:])
    r_pred[i,:] = temp_r[-1,:]

# %% 计算分箱的firing rate
start_idx = df['start']
end_idx = df['stop']

r_pred_bin = np.zeros([temp_input.shape[1], r_in.shape[1], start_idx.shape[0]])
r_GT_bin = np.zeros([temp_input.shape[1], r_in.shape[1], start_idx.shape[0]])
for i in range(start_idx.shape[0]):
    temp_pred_bin = r_pred[start_idx[i][0]:end_idx[i][0],:]
    temp_GT_bin = r_in[start_idx[i][0]:end_idx[i][0],:]

    bin_times = temp_input[start_idx[i][0]:end_idx[i][0],:].astype(int).T
    denom = np.sum(bin_times, axis=1)[:,np.newaxis]
    denom[denom==0] = 1

    temp_bin = bin_times @ temp_pred_bin
    temp_bin = (temp_bin/denom)
    r_pred_bin[:,:,i] = temp_bin

    temp_bin = bin_times @ temp_GT_bin
    temp_GT_bin = (temp_bin/denom)
    r_GT_bin[:,:,i] = temp_GT_bin

# %% 计算预测误差
r_pred_mean = np.nanmean(r_pred_bin, axis=2)
r_GT_mean = np.nanmean(r_GT_bin, axis=2)
gt_pred_rmse = np.sqrt(np.mean((r_GT_mean - r_pred_mean)**2, axis=1))
gt_pred_rmse = np.mean(gt_pred_rmse)
print('新皮层RNN预测RMSE:', gt_pred_rmse)

# %% 可视化结果
fig = plt.figure(figsize=(12, 6))
ax = fig.add_subplot(121)
sorted_bin = r_GT_mean[:,np.argsort(np.argmax(r_GT_mean, axis=0))]
im1 = ax.imshow(sorted_bin.T, cmap='jet', interpolation='none', aspect='auto')
ax.set_title('新皮层 - 真实Firing Rate')
ax.set_xlabel('位置分箱')
ax.set_ylabel('神经元索引')
plt.colorbar(im1, ax=ax, label='Firing Rate')

ax = fig.add_subplot(122)
sorted_bin = r_pred_mean[:,np.argsort(np.argmax(r_GT_mean, axis=0))]
im2 = ax.imshow(sorted_bin.T, cmap='jet', interpolation='none', aspect='auto')
ax.set_title('新皮层 - 预测Firing Rate (RMSE: %.6f)' % gt_pred_rmse)
ax.set_xlabel('位置分箱')
ax.set_ylabel('神经元索引')
plt.colorbar(im2, ax=ax, label='Firing Rate')

plt.tight_layout()
plt.savefig('neocortex_prediction.png')
plt.show()

# %% 保存结果
np.savez('neocortex_rnn_results.npz', 
         r_pred=r_pred, 
         r_GT=r_in, 
         r_pred_mean=r_pred_mean, 
         r_GT_mean=r_GT_mean, 
         rmse=gt_pred_rmse, 
         w_rec=w_rec, 
         w_in=w_in)

print('新皮层RNN分析完成！结果已保存到 neocortex_rnn_results.npz')
