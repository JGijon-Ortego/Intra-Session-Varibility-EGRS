function shannonEntropy = computeShannonEntropy(PSD, f, band)
%COMPUTESHANNONENTROPY Computes the normalized Shannon entropy of the
% given distribution.
%
%   shannonEntropy = computeShannonEntropy(PSD, f, band) computes the
%   Shannon entropy of the distribution PSD, indexed by the frequency
%   vector f, between the frequencies given in band ([fmin fmax]).

% Indices of the frequencies within the band of interest
bandIdx = find(f >= band(1) & f <= band(2));

% Check whether the PSD is made of NaNs (artifact)
if isnan(PSD(bandIdx(3)))
    shannonEntropy = NaN;
else
    shannonEntropy = [];

    % Total power in the band
    totalPower = sum(PSD(bandIdx));
    % Probability density function in the band
    pdf = PSD(bandIdx) / totalPower;

    % Normalized Shannon entropy
    shannonEntropy = -nansum(pdf(:) .* log(pdf(:))) / log(length(pdf(:)));
end
end
