# Intra-session variability of EEG signal features in externally-guided and classical resting state paradigms

MATLAB functions used in the paper *Intra-session variability of EEG signal features in externally-guided and classical resting state paradigms* (Gijón-Ortego J, et al.; journal, year and DOI to be added).

The code computes the EEG signal features analyzed in the study and the statistics used to quantify their variability across epochs.

## Contents

**'FC/'** Functional connectivity
- 'computeAEC': amplitude envelope correlation.
- 'computePLI': phase lag index.
- 'orthogonalizeChannels': pairwise orthogonalization to remove zero-lag leakage (O'Neill et al., 2015).

**'Graph/'** Graph metrics
- 'computeNodeStrength': node strength (sum of the weights of the links of each node).
- 'computeClusteringCoefficient': nodal and global clustering coefficient.
- 'computePathLength': characteristic path length, efficiency, eccentricity, radius and diameter.
- 'surrogateGraph': surrogate graph obtained by shuffling connections.

**'Local_Activ/'** Local activation and complexity
- 'computeMedianFrequency': median frequency of the power spectrum.
- 'computeShannonEntropy': normalized Shannon entropy of the power spectrum.
- 'computeLZC': Lempel-Ziv complexity.
- 'computeSampleEntropy': sample entropy.

**'Stats/'** Statistics
- 'computeRobustCV': robust coefficient of variation, RCV (1.4826 x MAD / median).
- 'dif_RCV': mean paired difference of the RCV between two conditions (statistic used in the bootstrap comparison).
- 'fdr_BH': Benjamini-Hochberg false discovery rate correction.

**'Examples/'** Example scripts (the signals and the results are not included; each script describes the input variables it expects)
- 'example_pipeline': how to compute the parameters (local activation, functional connectivity and graph metrics) for each epoch, and the RCV across epochs.
- 'example_statistics': statistical analysis used to compare the RCV between conditions (bootstrap of the mean paired difference, 95% confidence intervals and Benjamini-Hochberg FDR correction). It starts from the RCV results, which must be downloaded or computed beforehand with 'example_pipeline'.

## Requirements

- MATLAB with the Signal Processing and Statistics and Machine Learning toolboxes (the example also uses the Parallel Computing toolbox).
- [Brain Connectivity Toolbox](https://sites.google.com/site/bctnet), needed by 'computePathLength'.
- MATLAB R2023a or later for 'fdr_BH' (uses 'cummin' with the 'reverse' option).

## Usage

    addpath(genpath('path/to/repo'));

    aec = computeAEC(signal);                    % signal: samples x channels
    pli = computePLI(signal, size(signal, 2));
    str = computeNodeStrength(pli);
    [lambda, efficiency] = computePathLength(pli);
    cv  = computeRobustCV(param, 4);             % robust CV across the 4th dimension (epochs)

Signals must be band-pass filtered before computing AEC, PLI or Lempel-Ziv complexity.

## Citation

If you use this code, please cite the paper above and this repository (see 'CITATION.cff').

## License

MIT (see 'LICENSE'). Parts derived from the Brain Connectivity Toolbox are under the BSD 3-Clause License (see 'THIRD_PARTY_LICENSES.md').

## Credits

Repository author: Jorge Gijón-Ortego.

Some functions are adaptations of code by other authors:
- Javier Gomez-Pilar: 'surrogateGraph', 'computeClusteringCoefficient', 'computePathLength'.
- Jesús Poza Crespo: 'computeMedianFrequency', 'computeShannonEntropy'.
- Daniel Abasolo and Carlos Gómez Peña: 'computeLZC'.
- V. Martínez-Cagigal: 'computeSampleEntropy', 'fdr_BH'.
- Rubinov and Sporns (Brain Connectivity Toolbox): basis of the clustering coefficient and path length functions; 'computeNodeStrength' is equivalent to 'strengths_und'.
