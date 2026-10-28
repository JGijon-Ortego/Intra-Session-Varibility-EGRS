function orthSignal = orthogonalizeChannels(signal)
%ORTHOGONALIZECHANNELS Pairwise orthogonalization of each channel with
% respect to all the others (to remove zero-lag leakage).
% Based on O'Neill et al. (2015).
%
% Inputs:
%   signal: matrix (samples x channels).
%
% Outputs:
%   orthSignal: array (samples x channels x channels). orthSignal(:,j,i)
%               is channel j orthogonalized with respect to channel i.
%               orthSignal(:,i,i) is the original channel i.

numChannels = size(signal, 2);
orthSignal = NaN(size(signal, 1), numChannels, numChannels);

for refChannel = 1:numChannels
    orthSignal(:, refChannel, refChannel) = signal(:, refChannel);
    for targetChannel = 1:numChannels
        if refChannel ~= targetChannel
            ref = signal(:, refChannel);
            target = signal(:, targetChannel);
            beta = target.' * pinv(ref.');
            orthSignal(:, targetChannel, refChannel) = target - beta * ref;
        end
    end
end
end
