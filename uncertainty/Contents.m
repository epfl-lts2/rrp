% Global and Local Uncertainty Principles for Signals on Graphs - Code
%
%   This package contains the code to reproduce all the figures of the
%   paper: 
%   
%   Global and Local Uncertainty Principles for Signals on Graphs
%
%   Authors: Nathanael Perraudin, Benjamin Ricaud, David I Shuman, Pierre
%   Vandergheynst 
%
%   ArXiv: http://arxiv.org/abs/1603.03030
%
%   Availlable experiments:
%    GLOBAL_ILLUSTRATION           -  Laplacian eigenvector localization
%    COMET_COHERENCE               -  Evolution of the coherence for the comet graph
%    MODIFIED_PATH_EIGENVECTORS    -  Display some eigenvectors for a modified path
%    MODIFIED_PATH_COHERENCE       -  Illustration on the modified path graph 
%    NORM_OF_LOCALIZATION_OPERATOR -  Illustration of the localization operator
%    SPREAD_GABOR                  -  Display the spread of different Gabor transforms
%    LOCAL_UNCERTAINTY             -  Local uncertainty illustration
%    ATOM_LOCALIZATION             -  Example of the localization of some atoms
%    MODIFIED_PATH_GABOR           -  Show global and local bounds on the modified path
%    ADAPTATED_SAMPLING            -  Adaptation of the sampling according to local uncertainy
%
%
%   Abstract of the paper
%   ---------------------
%
%   Uncertainty principles such as Heisenberg's provide limits on the
%   time-frequency concentration of a signal, and constitute an important
%   theoretical tool for designing and evaluating linear signal transforms.
%   Generalizations of such principles to the graph setting can inform
%   dictionary design for graph signals, lead to algorithms for
%   reconstructing missing information from graph signals via sparse
%   representations, and yield new graph analysis tools. While previous
%   work has focused on generalizing notions of spreads of a graph signal
%   in the vertex and graph spectral domains, our approach is to generalize
%   the methods of Lieb in order to develop uncertainty principles that
%   provide limits on the concentration of the analysis coefficients of any
%   graph signal under a dictionary transform whose atoms are jointly
%   localized in the vertex and graph spectral domains. One challenge we
%   highlight is that due to the inhomogeneity of the underlying graph data
%   domain, the local structure in a single small region of the graph can
%   drastically affect the uncertainty bounds for signals concentrated in
%   different regions of the graph, limiting the information provided by
%   global uncertainty principles. Accordingly, we suggest a new way to
%   incorporate a notion of locality, and develop local uncertainty
%   principles that bound the concentration of the analysis coefficients of
%   each atom of a localized graph spectral filter frame in terms of
%   quantities that depend on the local structure of the graph around the
%   center vertex of the given atom. Finally, we demonstrate how our
%   proposed local uncertainty measures can improve the random sampling of
%   graph signals.
%
%   References: perraudin2016global

% Undocumented code
%    UNCERTAINTY_BOUNDS            -  Compute the global uncertainty bound of some graphs
