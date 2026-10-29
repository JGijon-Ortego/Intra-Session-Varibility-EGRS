%% Example: statistical comparison of the robust CV between conditions
%% (bootstrap + Benjamini-Hochberg FDR correction).
%
% For each of the 40 parameters (see example_pipeline.m), this script
% compares the robust coefficient of variation (RCV) between the classical
% resting state (RS) and the externally-guided resting state (EGRS):
%   1. Globally (RCV averaged across ROIs), and
%   2. For each ROI.
% The statistic of each comparison is the mean paired difference of the RCV
% (dif_RCV) applied to the [subjects x 2] matrix [RS EGRS]. Its distribution
% is estimated with n_boots bootstrap resamples. The two-sided p-value is
% 2*min(p, 1-p), with p = proportion of bootstrap values <= 0, and the 95%
% confidence interval is given by the 2.5 and 97.5 percentiles. P-values are
% then corrected for multiple comparisons with the Benjamini-Hochberg FDR
% procedure (fdr_BH).
%
% INPUT DATA (not included in this repository)
% The robust CV (RCV) results must be downloaded or computed beforehand
% (see example_pipeline.m, which computes them with computeRobustCV), and
% loaded into the workspace with these names and dimensions:
%
%   cv_rob_rs_all:   [40 x subjects]         RCV of RS averaged across ROIs.
%   cv_rob_egrs_all: [40 x subjects]         RCV of EGRS averaged across ROIs.
%   cv_rob_rs:       [40 x subjects x ROIs]  RCV of RS for each ROI.
%   cv_rob_egrs:     [40 x subjects x ROIs]  RCV of EGRS for each ROI.
%
% where 40 is the number of parameters and ROIs = 68.
%
% Outputs:
%   ic95:             [40 x 2]         95% CI of the global comparison.
%   p_values:         [1 x 40]         p-values of the global comparison.
%   p_values_FDR:     FDR-adjusted p-values of the global comparison.
%   ic95_ROI:         [40 x ROIs x 2]  95% CI for each ROI.
%   p_values_ROI:     [40 x ROIs]      p-values for each ROI.
%   p_values_ROI_FDR: [40 x ROIs]      FDR-adjusted p-values for each ROI
%                                      (correction over all 40*ROIs tests).
%
% Requirements: this repository on the MATLAB path (dif_RCV, fdr_BH) and
% the Statistics and Parallel Computing toolboxes.

close all; clc
addpath(genpath('path/to/this/repository'));

%% LOAD THE RCV RESULTS HERE
% Download the RCV results and load them, e.g.:
% load('path/to/RCV_results.mat');   % must contain the four variables above

%% PARAMETERS
n_boots = 50000;     % number of bootstrap resamples
seed = 42;           % random seed (reproducibility)
rng(seed)

n_params = size(cv_rob_rs_all, 1);   % 40 parameters
n_ROIs = size(cv_rob_rs, 3);         % 68 ROIs

%% BOOTSTRAP COMPARISONS
ic95 = zeros(n_params, 2);
p_values = zeros(1, n_params);
ic95_ROI = zeros(n_params, n_ROIs, 2);
p_values_ROI = zeros(n_params, n_ROIs);

for param = 1:n_params
    tic
    fprintf('Parameter: %d\n', param)

    % Global comparison (RCV averaged across ROIs)
    data = [cv_rob_rs_all(param,:)' cv_rob_egrs_all(param,:)'];
    distribution_boot = bootstrp(n_boots, @(x) dif_RCV(x), data, Options=statset(UseParallel=true));
    p_val = mean(distribution_boot <= 0);
    ic95(param,:) = prctile(distribution_boot, [2.5 97.5]);
    p_values(param) = 2 * min(p_val, 1 - p_val);

    % Comparison for each ROI
    for ROI = 1:n_ROIs
        data_ROI = [squeeze(cv_rob_rs(param,:,ROI))' squeeze(cv_rob_egrs(param,:,ROI))'];
        distribution_boot_ROI = bootstrp(n_boots, @(x) dif_RCV(x), data_ROI, Options=statset(UseParallel=true));
        p_val_ROI = mean(distribution_boot_ROI <= 0);
        ic95_ROI(param,ROI,:) = prctile(distribution_boot_ROI, [2.5 97.5]);
        p_values_ROI(param,ROI) = 2 * min(p_val_ROI, 1 - p_val_ROI);
    end
    toc
end

%% MULTIPLE COMPARISONS CORRECTION (Benjamini-Hochberg FDR)
% Global comparison: 40 tests
p_values_FDR = fdr_BH(p_values, 0.05);

% ROI comparisons: 40 x 68 tests corrected together
p_values_ROI_FDR = reshape(p_values_ROI, [n_params*n_ROIs 1]);
p_values_ROI_FDR = fdr_BH(p_values_ROI_FDR, 0.05);
p_values_ROI_FDR = reshape(p_values_ROI_FDR, [n_params n_ROIs]);
