%% Example: local activation, connectivity and graph parameters per epoch,
%% and their robust CV across epochs.
%
% Example pipeline of the analysis. For each subject, ROI and epoch it
% computes 40 parameters in two conditions: externally-guided resting state
% (EGRS) and classical resting state (RS). The epoch-to-epoch variability of
% each parameter is then quantified with the robust coefficient of
% variation (computeRobustCV).
%
% INPUT DATA (not included in this repository)
% The signals must be loaded BEFORE running this script, as two arrays:
%
%   EGRS: [subjects x ROIs x samples]  EGRS (externally-guided) signals.
%   RS:   [subjects x ROIs x samples]  classical resting-state signals.
%
% Both arrays must have the same size. For example, 47 subjects x 68 ROIs x
% 90000 samples (90000 samples = 180 s at fs = 500 Hz). The number of
% samples must be a multiple of fs*epochLength.
% Each subject's signals are ROI time series, e.g.:
%
%   load('subject01.mat');          % contains 'rawdata' [ROIs x samples]
%   EGRS(1,:,:) = rawdata;
%
% Parameters (rows of param_egrs / param_rs; 6 bands each: delta, theta,
% alpha, beta1, beta2, gamma):
%   1        sample entropy
%   2        Lempel-Ziv complexity
%   3        Shannon entropy
%   4        median frequency
%   5-10     normalized clustering coefficient (AEC)
%   11-16    normalized path length (AEC)
%   17-22    node strength (AEC)
%   23-28    normalized clustering coefficient (PLI)
%   29-34    normalized path length (PLI)
%   35-40    node strength (PLI)
%
% Clustering coefficient and path length are normalized by the mean of
% nSurrogates surrogate graphs.
%
% Outputs:
%   param_egrs, param_rs:           [40 x subjects x ROIs x epochs].
%   cv_rob_egrs, cv_rob_rs:         robust CV across epochs [40 x subjects x ROIs].
%   cv_rob_egrs_all, cv_rob_rs_all: robust CV averaged across ROIs [40 x subjects].
%
% Requirements: this repository and the Brain Connectivity Toolbox (BCT) on
% the MATLAB path (used by computePathLength), and the Signal Processing,
% Statistics and Parallel Computing toolboxes.

% Note: do not use 'clear all' here, since it would delete EGRS and RS.
close all; clc

addpath(genpath('path/to/this/repository'));
addpath(genpath('path/to/BCT'));

%% LOAD SIGNALS HERE
% Load the arrays EGRS and RS described above (e.g., from the .mat files of
% each subject). The data are not available in this repository.
% Example with random data to test the script:
% EGRS = randn(2, 68, 5000);  RS = randn(2, 68, 5000);

%% PARAMETERS
fs = 500;
epochLength = 5;
m = 1;
r = 0.1;
band = [1 70];
nSurrogates = 50;
Den = 1;

% Band-pass FIR filters (one per frequency band)
liveFilterDelta=designfilt('bandpassfir','DesignMethod','window','Window','hamming','FilterOrder',1000,'CutoffFrequency1',1,'CutoffFrequency2',4,'SampleRate',fs);
liveFilterTheta=designfilt('bandpassfir','DesignMethod','window','Window','hamming','FilterOrder',1000,'CutoffFrequency1',4,'CutoffFrequency2',8,'SampleRate',fs);
liveFilterAlpha=designfilt('bandpassfir','DesignMethod','window','Window','hamming','FilterOrder',1000,'CutoffFrequency1',8,'CutoffFrequency2',13,'SampleRate',fs);
liveFilterBeta1=designfilt('bandpassfir','DesignMethod','window','Window','hamming','FilterOrder',1000,'CutoffFrequency1',13,'CutoffFrequency2',19,'SampleRate',fs);
liveFilterBeta2=designfilt('bandpassfir','DesignMethod','window','Window','hamming','FilterOrder',1000,'CutoffFrequency1',19,'CutoffFrequency2',30,'SampleRate',fs);
liveFilterGamma=designfilt('bandpassfir','DesignMethod','window','Window','hamming','FilterOrder',1000,'CutoffFrequency1',30,'CutoffFrequency2',70,'SampleRate',fs);

num_subjects = size(EGRS, 1);
num_ROIs = size(EGRS, 2);

    % Number of epochs (computed outside the parfor)
    n_epochs = size(EGRS,3)/(fs*epochLength);
    
    % Preallocation (required for parfor and for speed)
    param_egrs = zeros(40, num_subjects, num_ROIs, n_epochs);
    param_rs = zeros(40, num_subjects, num_ROIs, n_epochs);
    
    parfor subj = 1:num_subjects
        tic
        fprintf('Subject: %d\n', subj)
        
        % Local temporary variables to avoid indexing issues inside parfor
        temp_param_egrs = zeros(40, num_ROIs, n_epochs);
        temp_param_rs = zeros(40, num_ROIs, n_epochs);
        
        current_egrs = squeeze(EGRS(subj,:,:));
        current_rs = squeeze(RS(subj,:,:));
        
        % Preallocate to avoid dynamic growth
        eeg_egrs = zeros(num_ROIs, fs*epochLength, n_epochs);
        eeg_rs = zeros(num_ROIs, fs*epochLength, n_epochs);
        
        for epochs = 1:n_epochs
            idx1 = (epochs-1)*fs*epochLength+1;
            idx2 = epochs*fs*epochLength;
            eeg_egrs(:,:,epochs) = current_egrs(:, idx1:idx2);
            eeg_rs(:,:,epochs) = current_rs(:, idx1:idx2);
        end
        
        % Preallocate PSD to avoid size errors
        [PSD_test_egrs, f_test_egrs] = pwelch(squeeze(eeg_egrs(:,:,1))',[],[],[],fs);
        PSD_egrs = zeros(size(PSD_test_egrs,1), num_ROIs, n_epochs);
        f_egrs = zeros(size(f_test_egrs,1), 1, n_epochs);
        PSD_rs = zeros(size(PSD_test_egrs,1), num_ROIs, n_epochs);
        f_rs = zeros(size(f_test_egrs,1), 1, n_epochs);

        for epochs2 = 1:n_epochs
            [PSD_egrs(:,:,epochs2), f_egrs(:,:,epochs2)] = pwelch(squeeze(eeg_egrs(:,:,epochs2))',[],[],[],fs);
            [PSD_rs(:,:,epochs2), f_rs(:,:,epochs2)] = pwelch(squeeze(eeg_rs(:,:,epochs2))',[],[],[],fs);
        end

        for epochs3 = 1:n_epochs
            for ROI = 1:num_ROIs
                temp_param_egrs(1,ROI,epochs3) = computeSampleEntropy(squeeze(eeg_egrs(ROI,:,epochs3)),m,r);
                temp_param_rs(1,ROI,epochs3) = computeSampleEntropy(squeeze(eeg_rs(ROI,:,epochs3)),m,r);
                temp_param_egrs(2,ROI,epochs3) = computeLZC(squeeze(eeg_egrs(ROI,:,epochs3)));
                temp_param_rs(2,ROI,epochs3) = computeLZC(squeeze(eeg_rs(ROI,:,epochs3)));
                temp_param_egrs(3,ROI,epochs3) = computeShannonEntropy(squeeze(PSD_egrs(:,ROI,epochs3)),squeeze(f_egrs(:,:,epochs3)),band);
                temp_param_rs(3,ROI,epochs3) = computeShannonEntropy(squeeze(PSD_rs(:,ROI,epochs3)),squeeze(f_rs(:,:,epochs3)),band);
                temp_param_egrs(4,ROI,epochs3) = computeMedianFrequency(squeeze(PSD_egrs(:,ROI,epochs3)),squeeze(f_egrs(:,:,epochs3)),band);
                temp_param_rs(4,ROI,epochs3) = computeMedianFrequency(squeeze(PSD_rs(:,ROI,epochs3)),squeeze(f_rs(:,:,epochs3)),band);
            end
        end
        
        %% EGRS: re-referencing, filtering and functional connectivity
        EEG = double(current_egrs)';
        signal_average = nanmean(EEG,2);
        EEG_rereferenced = EEG - repmat(signal_average,1,size(EEG,2));
        EEG_zscored = zscore(EEG_rereferenced);

        cleansignalsFiltered = zeros(6, size(EEG_zscored,1), size(EEG_zscored,2));
        cleansignalsFiltered(1,:,:)=filtfilt(liveFilterDelta.Coefficients,Den,EEG_zscored);
        cleansignalsFiltered(2,:,:)=filtfilt(liveFilterTheta.Coefficients,Den,EEG_zscored);
        cleansignalsFiltered(3,:,:)=filtfilt(liveFilterAlpha.Coefficients,Den,EEG_zscored);
        cleansignalsFiltered(4,:,:)=filtfilt(liveFilterBeta1.Coefficients,Den,EEG_zscored);
        cleansignalsFiltered(5,:,:)=filtfilt(liveFilterBeta2.Coefficients,Den,EEG_zscored);
        cleansignalsFiltered(6,:,:)=filtfilt(liveFilterGamma.Coefficients,Den,EEG_zscored);

        % Preallocate as [Time x Channels x Epochs x Bands]
        vectorEEG_Filtered_egrs = zeros(fs*epochLength, num_ROIs, n_epochs, 6);
        for b=1:6
            for iEpoch = 1:n_epochs
                idx1 = (iEpoch-1)*fs*epochLength + 1;
                idx2 = iEpoch*fs*epochLength;
                % squeeze converts 1x500x68 into 500x68 (Time x Channels)
                vectorEEG_Filtered_egrs(:,:,iEpoch,b) = squeeze(cleansignalsFiltered(b,idx1:idx2, :));
            end
        end

        AEC_EGRS = zeros(num_ROIs, num_ROIs, n_epochs, 6);
        PLI_EGRS = zeros(num_ROIs, num_ROIs, n_epochs, 6);
        for b=1:6
            for trial=1:n_epochs
                data_trial = vectorEEG_Filtered_egrs(:,:,trial,b);
                signal_ort = orthogonalizeChannels(data_trial);
                conn = zeros(num_ROIs, num_ROIs);
                for nROI = 1:num_ROIs
                    AEC1 = computeAEC(squeeze(signal_ort(:,:,nROI)));
                    conn(:,nROI) = AEC1(:,nROI);
                end
                AECglobal = (triu(conn) + tril(conn).') / 2;
                AEC_EGRS(:,:,trial,b) = abs(triu(AECglobal,1)+AECglobal.');
                PLI_EGRS(:,:,trial,b) = computePLI(data_trial, num_ROIs);
            end
        end

        %% RS: re-referencing, filtering and functional connectivity
        EEG = double(current_rs)';
        signal_average = nanmean(EEG,2);
        EEG_rereferenced = EEG - repmat(signal_average,1,size(EEG,2));
        EEG_zscored = zscore(EEG_rereferenced);

        cleansignalsFiltered(1,:,:)=filtfilt(liveFilterDelta.Coefficients,Den,EEG_zscored);
        cleansignalsFiltered(2,:,:)=filtfilt(liveFilterTheta.Coefficients,Den,EEG_zscored);
        cleansignalsFiltered(3,:,:)=filtfilt(liveFilterAlpha.Coefficients,Den,EEG_zscored);
        cleansignalsFiltered(4,:,:)=filtfilt(liveFilterBeta1.Coefficients,Den,EEG_zscored);
        cleansignalsFiltered(5,:,:)=filtfilt(liveFilterBeta2.Coefficients,Den,EEG_zscored);
        cleansignalsFiltered(6,:,:)=filtfilt(liveFilterGamma.Coefficients,Den,EEG_zscored);

        % Preallocate as [Time x Channels x Epochs x Bands]
        vectorEEG_Filtered_rs = zeros(fs*epochLength, num_ROIs, n_epochs, 6);
        for b=1:6
            for iEpoch = 1:n_epochs
                idx1 = (iEpoch-1)*fs*epochLength + 1;
                idx2 = iEpoch*fs*epochLength;
                % squeeze converts 1x500x68 into 500x68 (Time x Channels)
                vectorEEG_Filtered_rs(:,:,iEpoch,b) = squeeze(cleansignalsFiltered(b,idx1:idx2, :));
            end
        end

        AEC_RS = zeros(num_ROIs, num_ROIs, n_epochs, 6);
        PLI_RS = zeros(num_ROIs, num_ROIs, n_epochs, 6);
        for b=1:6
            for trial=1:n_epochs
                data_trial = vectorEEG_Filtered_rs(:,:,trial,b);
                signal_ort = orthogonalizeChannels(data_trial);
                conn = zeros(num_ROIs, num_ROIs);
                for nROI = 1:num_ROIs
                    AEC1 = computeAEC(squeeze(signal_ort(:,:,nROI)));
                    conn(:,nROI) = AEC1(:,nROI);
                end
                AECglobal = (triu(conn) + tril(conn).') / 2;
                AEC_RS(:,:,trial,b) = abs(triu(AECglobal,1)+AECglobal.');
                PLI_RS(:,:,trial,b) = computePLI(data_trial, num_ROIs);
            end
        end

        %% GRAPH METRICS AND SURROGATES (no clear statements inside parfor)
        for b=1:6
            for trial=1:n_epochs
                % AEC, EGRS
                [~,ClC_nodal] = computeClusteringCoefficient(squeeze(AEC_EGRS(:,:,trial,b)));
                [~,~,~,~,~,D] = computePathLength((squeeze(AEC_EGRS(:,:,trial,b))));
                PL_nodal = nanmean(D);

                PL_surrogate = zeros(num_ROIs, nSurrogates);
                ClC_surrogate = zeros(num_ROIs, nSurrogates);
                for iSurrogate = 1:nSurrogates
                    surrogateMatrix = surrogateGraph(squeeze(AEC_EGRS(:,:,trial,b)));
                    [~,~,~,~,~,D_surrogate] = computePathLength(surrogateMatrix);
                    PL_surrogate(:,iSurrogate) = nanmean(D_surrogate);

                    surrogateMatrix = surrogateGraph(squeeze(AEC_EGRS(:,:,trial,b)));
                    [~,ClC_nodal_surrogate] = computeClusteringCoefficient(surrogateMatrix); 
                    ClC_surrogate(:,iSurrogate) = ClC_nodal_surrogate;
                end
                temp_param_egrs(4+b,:,trial) = ClC_nodal./nanmean(ClC_surrogate,2);
                temp_param_egrs(10+b,:,trial) = PL_nodal'./nanmean(PL_surrogate,2);
                temp_param_egrs(16+b,:,trial) = computeNodeStrength(squeeze(AEC_EGRS(:,:,trial,b)));

                % AEC, RS
                [~,ClC_nodal] = computeClusteringCoefficient(squeeze(AEC_RS(:,:,trial,b)));
                [~,~,~,~,~,D] = computePathLength((squeeze(AEC_RS(:,:,trial,b))));
                PL_nodal = nanmean(D);
                for iSurrogate = 1:nSurrogates
                    surrogateMatrix = surrogateGraph(squeeze(AEC_RS(:,:,trial,b)));
                    [~,~,~,~,~,D_surrogate] = computePathLength(surrogateMatrix);
                    PL_surrogate(:,iSurrogate) = nanmean(D_surrogate);

                    surrogateMatrix = surrogateGraph(squeeze(AEC_RS(:,:,trial,b)));
                    [~,ClC_nodal_surrogate] = computeClusteringCoefficient(surrogateMatrix); 
                    ClC_surrogate(:,iSurrogate) = ClC_nodal_surrogate; 
                end
                temp_param_rs(4+b,:,trial) = ClC_nodal./nanmean(ClC_surrogate,2);
                temp_param_rs(10+b,:,trial) = PL_nodal'./nanmean(PL_surrogate,2);
                temp_param_rs(16+b,:,trial) = computeNodeStrength(squeeze(AEC_RS(:,:,trial,b)));

                % PLI, EGRS
                [~,ClC_nodal] = computeClusteringCoefficient(squeeze(PLI_EGRS(:,:,trial,b)));
                [~,~,~,~,~,D] = computePathLength((squeeze(PLI_EGRS(:,:,trial,b))));
                PL_nodal = nanmean(D); 
                for iSurrogate = 1:nSurrogates
                    surrogateMatrix = surrogateGraph(squeeze(PLI_EGRS(:,:,trial,b)));
                    [~,~,~,~,~,D_surrogate] = computePathLength(surrogateMatrix); 
                    PL_surrogate(:,iSurrogate) = nanmean(D_surrogate);

                    surrogateMatrix = surrogateGraph(squeeze(PLI_EGRS(:,:,trial,b)));
                    [~,ClC_nodal_surrogate] = computeClusteringCoefficient(surrogateMatrix); 
                    ClC_surrogate(:,iSurrogate) = ClC_nodal_surrogate; 
                end
                temp_param_egrs(22+b,:,trial) = ClC_nodal./nanmean(ClC_surrogate,2); 
                temp_param_egrs(28+b,:,trial) = PL_nodal'./nanmean(PL_surrogate,2); 
                temp_param_egrs(34+b,:,trial) = computeNodeStrength(squeeze(PLI_EGRS(:,:,trial,b)));

                % PLI, RS
                [~,ClC_nodal] = computeClusteringCoefficient(squeeze(PLI_RS(:,:,trial,b)));
                [~,~,~,~,~,D] = computePathLength((squeeze(PLI_RS(:,:,trial,b))));
                PL_nodal = nanmean(D); 
                for iSurrogate = 1:nSurrogates
                    surrogateMatrix = surrogateGraph(squeeze(PLI_RS(:,:,trial,b)));
                    [~,~,~,~,~,D_surrogate] = computePathLength(surrogateMatrix);
                    PL_surrogate(:,iSurrogate) = nanmean(D_surrogate);

                    surrogateMatrix = surrogateGraph(squeeze(PLI_RS(:,:,trial,b)));
                    [~,ClC_nodal_surrogate] = computeClusteringCoefficient(surrogateMatrix); 
                    ClC_surrogate(:,iSurrogate) = ClC_nodal_surrogate; 
                end
                temp_param_rs(22+b,:,trial) = ClC_nodal./nanmean(ClC_surrogate,2); 
                temp_param_rs(28+b,:,trial) = PL_nodal'./nanmean(PL_surrogate,2); 
                temp_param_rs(34+b,:,trial) = computeNodeStrength(squeeze(PLI_RS(:,:,trial,b)));
            end
        end

        % Assign the temporary variables to the global arrays at the end of the worker
        param_egrs(:, subj, :, :) = temp_param_egrs;
        param_rs(:, subj, :, :) = temp_param_rs;
        toc
    end % END OF PARFOR

%% ROBUST CV ACROSS EPOCHS (4th dimension of param)
cv_rob_egrs = computeRobustCV(param_egrs, 4);
cv_rob_egrs_all = mean(cv_rob_egrs, 3);   % averaged across ROIs

cv_rob_rs = computeRobustCV(param_rs, 4);
cv_rob_rs_all = mean(cv_rob_rs, 3);
