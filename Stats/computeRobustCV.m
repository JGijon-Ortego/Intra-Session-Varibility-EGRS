function robustCV = computeRobustCV(data, dim)
%COMPUTEROBUSTCV Computes the robust coefficient of variation (robust CV).
%
%   robustCV = computeRobustCV(data, dim) computes, along dimension dim,
%   the robust CV of data as:
%
%       robust CV = 1.4826 * MAD / median
%
%   where MAD = median(|x - median(x)|) is the median absolute deviation.
%   The constant 1.4826 (= 1/norminv(0.75)) makes the MAD a consistent
%   estimator of the standard deviation for normally distributed data, so
%   the robust CV is comparable with the classical CV (SD/mean) while being
%   insensitive to outliers. It is typically used to quantify the
%   variability of a parameter across epochs.
%
% Inputs:
%   data: numeric array of any size (e.g., parameters x subjects x ROIs x epochs).
%   dim:  (optional) dimension along which the robust CV is computed
%         (e.g., the epochs dimension). Default: first non-singleton
%         dimension.
%
% Outputs:
%   robustCV: array with the same size as data, except that dimension dim
%             has length 1.
%
% Example:
%   % param: [parameters x subjects x ROIs x epochs]
%   cv = computeRobustCV(param, 4);   % robust CV across epochs
%   cvROIs = mean(cv, 3);             % averaged across ROIs
%
% Author: Jorge Gijón-Ortego

if nargin < 2
    dim = find(size(data) ~= 1, 1);
    if isempty(dim), dim = 1; end
end

medianValue = median(data, dim);
mad = median(abs(medianValue - data), dim);
robustCV = 1.4826 * (mad ./ medianValue);
end
