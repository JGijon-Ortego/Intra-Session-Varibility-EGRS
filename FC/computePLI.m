function PLI = computePLI(signal, numChannels)
%COMPUTEPLI Computes the Phase Lag Index (PLI) between all channel pairs.
%
% Inputs:
%   signal:      matrix (time x channels). Filter the signal beforehand.
%   numChannels: number of channels (columns of signal).
%
% Outputs:
%   PLI: channels x channels matrix of phase lag index values.

phaseSignal = angle(hilbert(signal));

% Build all channel pairs in a vectorized way
angles1 = reshape(repmat(phaseSignal, numChannels, 1), size(phaseSignal, 1), numChannels * numChannels);
angles2 = repmat(phaseSignal, 1, numChannels);

PLIvector = abs(mean(sign(sin(angles1 - angles2))));
PLI = reshape(PLIvector, [numChannels, numChannels]);

end
