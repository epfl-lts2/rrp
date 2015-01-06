%   A fast Griffin lim algorithm
%
%   AN EXTENDED GRIFFIN LIM ALGORTITHM
%
%   Paper: Perraudin Nathanael, Balazs Peter
%   
%   Demonstration matlab file:  Perraudin Nathanael
%
%   In this paper, we present a new algorithm to estimate a signal from its
%   short-time Fourier transform modulus (STFTM). This algorithm is
%   computationally simple and is obtained by an acceleration of the
%   well-known Griffin-Lim algorithm (GLA). Before deriving the algorithm,
%   we will give a new interpretation of the GLA and formulate the phase
%   recovery problem in an optimization form. We then present some
%   experimental results where the new algorithm is tested on various
%   signals. It shows not only significant improvement in speed of
%   convergence but it does as well recover the signals with a smaller
%   error than the traditional GLA.
%
%   Conference paper:
%   https://lts2research.epfl.ch/unlocbox/notes/unlocbox-note-007.pdf
%
%   Cite this paper
%   https://lts2research.epfl.ch/unlocbox/notes/unlocbox-note-007.bib
%
%  Availlable experiments
%    RR_ALPHA - Test the influence of the parameter alpha
%    RR_SIGNALS - Phase reconstruction experiments on different signals
%    RR_SPECTRO_MULT - Spectrogram multiplication experiment on different signals
%
%  For help, bug reports, suggestions etc. please send email to
%  unlocbox-help@lists.sourceforge.net
%