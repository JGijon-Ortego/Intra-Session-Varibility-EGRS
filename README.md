# EEG/MEG analysis functions (MATLAB)

MATLAB functions for functional connectivity, graph-theoretical analysis, local activation/complexity and statistics of EEG/MEG signals.

## Structure and functions

|---|---|---|
| 'FC' | ‘computeAEC’ | Amplitude envelope correlation between channels |
| ‘FC’ | ‘computePLI’ | Phase lag index between all channel pairs |
| ‘FC’ | ‘orthogonalizeChannels’ | Pairwise orthogonalization to remove zero-lag leakage (O'Neill et al., 2015) |
| ‘Graph’ | ‘computeNodeStrength’ | Node strength (weighted degree), undirected graphs |
| ‘Graph’ | ‘computeClusteringCoefficient’ | Nodal and global clustering coefficient (weighted, directed) |
| ‘Graph’ | ‘computePathLength’ | Characteristic path length, efficiency, eccentricity, radius, diameter |
| ‘Graph’ | ‘surrogateGraph’ | Surrogate graph by shuffling connections |
| ‘Local_Activ’ | ‘computeMedianFrequency’ | Median frequency of a PSD within a band |
| ‘Local_Activ’ | ‘computeShannonEntropy’ | Normalized Shannon entropy of a PSD within a band |
| ‘Local_Activ’ | ‘computeLZC’ | Lempel-Ziv complexity |
| ‘Local_Activ’ | ‘computeSampleEntropy’ | Sample entropy (SampEn) |
| ‘Stats’ | ‘fdr_BH’ | Benjamini-Hochberg FDR correction for multiple comparisons |
| ‘Stats’ | ‘computeRobustCV’ | Robust coefficient of variation, 1.4826 x MAD / median, along a chosen dimension (e.g., across epochs) |

Signals must be filtered in the band of interest before computing AEC, PLI or LZC.

## Installation

‘‘‘matlab
addpath(genpath('path/to/repo'));
‘‘‘

## Dependencies

- MATLAB with Signal Processing Toolbox (‘hilbert’) and Statistics and Machine Learning Toolbox (‘corr’, ‘pdist’, ‘nansum’ in older versions).
- [Brain Connectivity Toolbox](https://sites.google.com/site/bctnet) (‘distance_wei’, ‘charpath’), required by ‘computePathLength’. Install it separately and add it to the MATLAB path.
- ‘fdr_BH’ uses ‘cummin’ with the ‘'reverse'‘ option (MATLAB R2023a or later).
- ‘computeLZC’ uses ‘findstr’, which is deprecated in recent MATLAB versions but still available.

## Usage

‘‘‘matlab
aec = computeAEC(signal);                       % signal: samples x channels
pli = computePLI(signal, size(signal,2));
[lambda, eff] = computePathLength(pli);
[clcMean, clcNodal] = computeClusteringCoefficient(pli);
lzc = computeLZC(signal(:,1));
mf  = computeMedianFrequency(psd, f, [8 13]);
[pAdj, ~, h] = fdr_BH(pvalues, 0.05);
cv  = computeRobustCV(param, 4);                % param: [parameters x subjects x ROIs x epochs]
‘‘‘

## Citation

See ‘CITATION.cff’. A DOI is issued through Zenodo for each release.

## License

MIT (see ‘LICENSE’) for the original code in this repository. Portions derived from the Brain Connectivity Toolbox are under the BSD 3-Clause License (see ‘THIRD_PARTY_LICENSES.md’).

## Credits

Repository author: Jorge Gijón-Ortego (‘computeRobustCV’ and the repository organization).

Several functions are adaptations of code written by other authors:

| Function | Original author(s) |
|---|---|
| ‘computeNodeStrength’, ‘surrogateGraph’ | Javier Gomez-Pilar |
| ‘computeClusteringCoefficient’, ‘computePathLength’ | Javier Gomez-Pilar; based on the Brain Connectivity Toolbox (Rubinov and Sporns) |
| ‘computeMedianFrequency’, ‘computeShannonEntropy’ | Jesus Poza Crespo |
| ‘computeLZC’ | Daniel Abasolo (original code); Carlos Gomez Pena (modifications) |
| `computeSampleEntropy`, `fdr_BH` | V. Martinez-Cagigal |
