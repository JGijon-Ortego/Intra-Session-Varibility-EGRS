function medianFrequency = computeMedianFrequency(PSD, f, band)
%COMPUTEMEDIANFREQUENCY Computes the median frequency of the distribution
% contained in PSD.
%
%   medianFrequency = computeMedianFrequency(PSD, f, band) computes the
%   median frequency of the distribution PSD, indexed by the frequency
%   vector f, between the frequencies given in band ([fmin fmax]).

medianFrequency = [];

% Indices of the frequencies within the band of interest
bandIdx = find((f >= band(1)) & (f <= band(2)));
% Total power in the band
totalPower = sum(PSD(bandIdx));
% Cumulative power
cumulativePower = cumsum(PSD(bandIdx));

% Index at which half of the total power is reached
halfIdx = max(find(cumulativePower <= (totalPower / 2)));
medianIdx = bandIdx(halfIdx);
% If no index was selected, more than 50% of the power lies in the first value
if isempty(medianIdx)
    medianIdx = bandIdx(1);
end

% Frequency at the median power (median frequency)
medianFrequency = f(medianIdx);
end
