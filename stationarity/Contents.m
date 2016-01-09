% RRP - Stationary signal processing on graphs
%
%   Abstract
%   --------
%   
%   Graphs are a central tool in machine learning and information
%   processing as they allow to conveniently capture the structure of
%   complex datasets. In this context, it is of high importance to develop
%   flexible models of signals defined over graphs or networks. In this
%   paper, we generalize the traditional concept of wide sense stationarity
%   to signals defined over the vertices of arbitrary weighted undirected
%   graphs. We show that stationarity is intimately linked to statistical
%   invariance under a localization operator reminiscent of translation. We
%   prove that stationary signals are characterized by a well-defined Power
%   Spectral Density that can be efficiently estimated even for large
%   graphs. We leverage this new concept to derive Wiener-type estimation
%   procedures of noisy and partially observed signals and illustrate the
%   performance of this new model for denoising and regression.
%
%   Authors: Nathanael Perraudin and Pierre Vandergheynst
%
%   Date: January 2016
%
%
%   Contents
%   --------
%
%   Figures of the paper
%      ESTIMATION_PSD - Estimation of the power spectrum density
%      ESTIMATION_PSD_SCALABILITY - Study the scalability of the PSD estimation method
%      DEMO_GENERATE_SAMPLES - Generation of USPS digit
%      SYNTHETIC_WIENER_DECONVOLUTION - Wiener deconvolution experiment on a synthetic dataset
%      SYNTHETIC_WIENER_INPAINTING - Wiener inpainting experiment on a synthetic dataset
%      EXPERIMENT_MOLENE_HUMIDITY - In-painting on the Molene humidity dataset
%      EXPERIMENT_MOLENE_TEMPERATURE - In-painting on the Molene temperature dataset
%      EXPERIMENT_USPS_INPAINTING - In-painting on the USPS dataset
%
%
%   Functions 
%      GSP_EXPERIMENTAL_PSD - Experimental power density function
%      GSP_DESIGN_TRANSLATEs - Create a filterbank by uniformly translating a window
%      GSP_STATIONARITY_COV - Covariance matrix from graph stationary data
%      GSP_STATIONARITY_RATIO - Assert the stationarity level of some data
%      GSP_PSD_ESTIMATION - Estimation of the Power spectrum density
%      GSP_WINER_OPTIMIZATION -Solve wiener optimization problem
%      GSP_WINER_L2 - Solve wiener optimization problem with l2 fidelity term
%
%  For help, bug reports, suggestions etc. please send email to
%  nathanael (dot) perraudin (at) epfl (dot) ch

