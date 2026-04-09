预处理主要步骤，1. bin一下数据，bin_size 2. 用高斯卷积一下，sigma，3. 把数据标准化到-1，1
后面RNN拟合的时候，需要固定alpha = bin_size/sigma . 为了让RNN的时间参数和数据匹配上
这几个是matlab的函数哈