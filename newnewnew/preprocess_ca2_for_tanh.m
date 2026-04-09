function out = preprocess_ca2_for_tanh(X, Fs, varargin)
% PREPROCESS_CA2_FOR_TANH  Bin -> Gaussian smooth -> per-cell normalize to ~[-1,1].
%
% out = preprocess_ca2_for_tanh(X, Fs, 'BinSec',0.1,'SigmaSec',1.0,'TruncSigma',3)
%
% Inputs
%   X   : (T x N) Ca2+ traces
%   Fs  : sampling rate of X (Hz), i.e., frames per second
%
% Name-Value
%   'BinSec'       (0.1)   Bin duration (s)
%   'SigmaSec'     (1.0)   Gaussian kernel sigma (s) applied AFTER binning
%   'TruncSigma'   (3)     Kernel half-width in sigmas (3 => ±3σ)
%   'Eps'          (1e-8)  Small constant to avoid divide-by-zero
%
% Outputs (struct)
%   out.X_bin      : (Tb x N) binned traces
%   out.X_smooth   : (Tb x N) smoothed traces
%   out.X_norm     : (Tb x N) normalized traces (~[-1,1])
%   out.t_bin      : (Tb x 1) time vector for binned data (s)
%   out.kernel     : gaussian kernel (in bins) used for convolution
%   out.params     : parameters used
%
% Notes
% - Binning uses mean within each bin.
% - Smoothing uses normalized Gaussian kernel with sum(kernel)=1 (preserves scale).
% - Normalization: X_centered = X - mean(X); then divide by max(abs()) per cell.

% -------- parse inputs --------
p = inputParser;
p.addParameter('BinSec', 0.1);
p.addParameter('SigmaSec', 1.0);
p.addParameter('TruncSigma', 3);
p.addParameter('Eps', 1e-8);
p.parse(varargin{:});
opts = p.Results;

[T, N] = size(X);

% -------- 1) binning --------
binSec   = opts.BinSec;
binSamp  = max(1, round(binSec * Fs));          % samples/frames per bin
Tb       = floor(T / binSamp);                  % number of full bins

X_use = X(1:Tb*binSamp, :);
X_reshaped = reshape(X_use, binSamp, Tb, N);    % (binSamp x Tb x N)
X_bin = squeeze(mean(X_reshaped, 1));           % (Tb x N)

% time vector for binned signal (bin centers)
t_bin = ((0:Tb-1)' + 0.5) * (binSamp / Fs);

% -------- 2) Gaussian convolution (in bins) --------
dt_bin = binSamp / Fs;                          % actual bin duration (s)
sigma_bins = opts.SigmaSec / dt_bin;            % sigma in bins
halfWidth = max(1, ceil(opts.TruncSigma * sigma_bins));

xk = (-halfWidth:halfWidth);
kernel = exp(-0.5 * (xk / sigma_bins).^2);
kernel = kernel / sum(kernel);                  % normalize to sum=1

% Convolve along time for each cell (same length output)
X_smooth = zeros(size(X_bin));
for i = 1:N
    X_smooth(:, i) = conv(X_bin(:, i), kernel, 'same');
end

% -------- 3) per-cell normalize to ~[-1,1] --------
mu = mean(X_smooth, 1);                         % (1 x N)
X_center = X_smooth - mu;

scale = max(abs(X_center), [], 1);              % max absolute per cell
scale = max(scale, opts.Eps);

X_norm = X_center ./ scale;                     % roughly in [-1,1]

% -------- package outputs --------
out = struct();
out.X_bin    = X_bin;
out.X_smooth = X_smooth;
out.X_norm   = X_norm;
out.t_bin    = t_bin;
out.kernel   = kernel(:);
out.params   = struct('Fs', Fs, ...
                      'BinSec', binSec, ...
                      'BinSamp', binSamp, ...
                      'SigmaSec', opts.SigmaSec, ...
                      'SigmaBins', sigma_bins, ...
                      'TruncSigma', opts.TruncSigma, ...
                      'Eps', opts.Eps);
end
