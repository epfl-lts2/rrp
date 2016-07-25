% RRP - Audio in-painting
%
%   This package requires the LTFAT and the GSPBox toolbox to be executed.
%
%   Abstract
%   --------
%   
%   In this contribution, we present a method to compensate for long
%   duration data gaps in audio signals, in particular music. To achieve
%   this task, a similarity graph is constructed, based on a short-time
%   Fourier analysis of reliable signal segments, e.g. the uncorrupted
%   remainder of the music piece, and the temporal regions adjacent to the
%   unreliable section of the signal. A suitable candidate segment is then
%   selected through an optimization scheme and smoothly inserted into the
%   gap.
%
%   * Paper   : Audio inpainting with similarity graphs
%   * Authors : Nathanael Perraudin, Nicki Holighaus, Piotr Majdak, Peter Balazs 
%   * Date    : June 2016 
%   * ArXiv   : http://arxiv.org/abs/1607.06667
%   * Online demo : https://lts2.epfl.ch/web-audio-inpainting/
%
%   In order to run this script you need to install the LTFAT and the
%   GSPBox toolboxes. You can download them at:
%
%   * http://ltfat.github.io/
%   * https://lts2.epfl.ch/gsp/
%
%
%   Contents
%   --------
%
%   A demo is availlable throught the *file ai_demo_inpaint*. The function
%   *ai_conf* contains the default parameters. 
%
%   Audio inpainting functions
%      AI_AUDIO_INPAINT     - Inpaint an audio signal
%      AI_TIME_AUDIO_GRAPH  - Create a graph of similarities from an audio signal
%      AI_AUDIO_INPAINT_MULTIPLE_HOLES - Inpaint an audio signal with multiple holes
%      AI_PLOT_GRAPH        - Plot the transitions graph
%      AI_START             - Start the audio inpainting box
%      AI_COMPUTE_FEATURES  - Compute a time-frequency features matrix
%      AI_FIND_TRANSITIONS  - Find two optimal transitions
%
%  For help, bug reports, suggestions etc. please send email to
%  nathanael (dot) perraudin (at) epfl (dot) ch
%


