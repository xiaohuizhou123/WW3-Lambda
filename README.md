

# WW3-Lambda
Wave breaking model within the framework of Phillips' Lambda implemented in WaveWatchIII.

Romero, Leonel. 2019. “Distribution of Surface Wave Breaking Fronts.” Geophysical Research Letters, September, 2019GL083408. https://doi.org/10.1029/2019GL083408.

# Release v1.1.0 (April 19, 2023)
The latest release v1.1.0 fixes a major issue related to the direction of the modulation transfer function of the dissipation which originally used a scalar mean direction weighted by the energy spectrum which led to unrealistic growth of the spectral tail in conditions with misaligned winds and dominant waves. This latest release uses a scale-dependent mean wave direction weighed by the cumulative mean squared slope.
This fix ony updated the following files: w3srcxmd.ftn w3srcemd.ftn ww3_ounp.ftn  ww3_outp.ftn  gx_outp.ftn
