function AEC = computeAEC(signal)
%COMPUTEAEC Computes the amplitude envelope correlation (AEC) between channels.
%
% Inputs:
%   signal: matrix (samples x channels) on which the metric is computed.
%           Filter the signal beforehand (e.g., in the band of interest).
%
% Outputs:
%   AEC: channels x channels matrix of amplitude envelope correlations.

% Analytic signal via the Hilbert transform
analyticSignal = hilbert(signal);
envelope = abs(analyticSignal);

% Log of the squared envelope (power envelope)
logPowerEnvelope = log(envelope.^2);

% Linear correlation between the envelopes of all channel pairs
AEC = corr(logPowerEnvelope, logPowerEnvelope);
end
