%ADAPTATED_SAMPLING Adaptation of the sampling according to local uncertainy
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
%   This experiment
%   ---------------
%
%   In order to motivate our uncertainty Theorem from a practical signal
%   processing point of view, we use it to optimize the sampling of a
%   signal over a graph. To asses the quality of the sampling, we solve a
%   small inpainting problem where only a part of a signal is measured and
%   the goal is to reconstruct the entire signal. Assuming that the signal
%   varies smoothly in the vertex domain, we can formulate the inverse
%   problem as:
%
%   .. argmin_x  x^T L x      s. t.    y = Mx,
%
%   .. math:: \mathop{\rm argmin}_x  x^T L x \hspace{0.25cm} \text{ s. t. } \hspace{0.25cm} y = Mx,
%
%   where $y$ is the observed signal, $M$ the inpainting masking operator
%   and $x^T L x$ the graph Tikhonov regularizer ($L$ being the
%   Laplacian). In order to generate the original signal, we filter
%   Gaussian noise on the graph with a low pass kernel $\hat{h}$. The
%   frequency content of the resulting signal will be close to the shape of
%   the filter $\hat{h}$. For this example, we use the low pass kernel
%   $\hat{h}(x) = \frac{1}{1+\frac{100}{\lmax} x}$ to generate the smooth
%   signal.
%
%   For a given number of measurements, the traditional idea is to randomly
%   sample the graph. Under that strategy, the measurements are distributed
%   across the network. Alternatively, we can use our local uncertainty
%   principles to create an adapted mask. The intuitive idea that nodes
%   with less uncertainty (higher local sparsity values) should be sampled
%   with higher probability because their value can be inferred less easily
%   from other nodes. Another way to picture this fact is the following.
%   Imagine that we want to infer a quantity over a random sensor network.
%   In the more densely populated parts of the network, the measurements
%   are more correlated and redundant. As result, a lower sampling rate is
%   necessary. On the contrary, in the parts where there are fewer sensors,
%   the information has less redundancy and a higher sampling rate is
%   necessary. The heat kernel $\hat{g}(x)=e^{-\tau x}$ is a convenient
%   choice to probe the local uncertainty of a graph, because
%   $\hat{g}^{2}(x)=e^{-2\tau~x}$ is also a heat kernel, resulting in a
%   sparsity level depending only on $\|T_{j}g^{2}\|_{2}$. Indeed we have
%   $||T_{j}g^{2}||_{1}=\sqrt{N}$. The local uncertainty bound of our
%   Theorem becomes:
%
%   .. O_{1}(j)= || T_j g^2 ||_1 / || T_j g^2 ||_2 = 
%              = \sqrt{N} / || T_j g^2 ||_2.
%
%   .. math:: O_{1}(j)=\frac{\|T_{j}g^{2}\|_{1}}{\|T_{j}g^{2}\|_{2}}=\frac{\sqrt{N}}{\|T_{j}g^{2}\|_{2}}.
%   
%   Based on this measure, we design a second random sampled mask with a
%   probability proportional to $||T_ig^2||_2$; that is, the higher the
%   overlap level at vertex $j$, the smaller the probability that vertex
%   $j$ is chosen as a sampling point, and vice-versa. For each sampling
%   ratio, we performed $100$ experiments and averaged the results. For
%   each experiment, we also randomly generated new graphs. The experiment
%   was carried out using open-source code: the UNLocBoX and the GSPBox.
%
%   For the sensor graph, we observe that our local measure of uncertainty
%   varies smoothly on the graph and is higher in the more dense part.
%   Thus, the likelihood of sampling poorly connected vertices is higher
%   than the likelihood of sampling well connected vertices. For the
%   community graph, we observe that the uncertainty is highly related to
%   the size of the community. The larger the community, the larger the
%   uncertainty (or, equivalently, the smaller the local sparsity value).
%   In both cases, the adapted, non-uniform random sampling performs better
%   than random uniform sampling.
%
%
%   .. figure::
%
%      Reconstructions resulting
%
%      Non uniform sampling performs better
%
%   .. figure::
%
%      Local uncertainty
%
%      The local uncertainty is inversly proportional to $||T_i g||_2$.
%  
%   .. figure::
%
%      Norm of the atoms
%
%      The random non-uniform sampling distribution is proportional to
%      $||T_i g||_2$.
%
%   
%   References: perraudin2016global
%



%% Initialisation
clear;
close all;

global SAVE;
global type;

% plotting paramter
paramplot.save = SAVE;
paramplot.savefig = 0;

paramplot.position = [100 100 300 220];
gsp_reset_seed(0)

%% Parameter


N = 300;
Ntest = 100;
Ns = 20:10:200; 
Nsmax = max(Ns);
h = @(G) @(x) 1./(1+100/G.lmax*x);
if ~numel(type)
    type = 'community'; % or sensor
end
%%
error1 = zeros(Ntest,length(Ns));
error2 = zeros(Ntest,length(Ns));


%%
for ii = 1:Ntest

    switch type
        case 'community'    
            G = gsp_community(N);
        case 'sensor'
            G = gsp_sensor(N);
        otherwise
            error('Unknown type of graph')
    end

   G = gsp_estimate_lmax(G);


    g = gsp_design_heat(G,8);
      ntig = gsp_norm_tig(G,@(x) g(x).^2);

    % Perform the sampling
    sample_non_uniform = randsampleWRW((1:N),Nsmax,ntig);    
    sample_uniform = randsample(N,Nsmax,0);

    s = gsp_filter_analysis(G,h(G),randn(G.N,1));
    for jj = 1:length(Ns);
        % Create the masks
        mask_non_uniform = zeros(N,1);
        mask_non_uniform(sample_non_uniform(1:Ns(jj))) = 1;
        mask_uniform = zeros(N,1);
        mask_uniform(sample_uniform(1:Ns(jj))) = 1;
        
        % Input signals
        y1 = mask_non_uniform .* s;
        y2 = mask_uniform .* s;

        param = struct;
        param.verbose = 0;
        param.exact = 1;
        
        % Solve the problem
        sol_non_uniform = gsp_regression_tik(G ,mask_non_uniform, y1 , 0, param );
        sol_uniform = gsp_regression_tik(G ,mask_uniform, y2 , 0, param );

%         sol1 = y1;
%         sol1(logical(1-mask1)) = sol1_m;
%         sol2 = y2;
%         sol2(logical(1-mask2)) = sol2_m;

        error1(ii,jj) = norm(sol_non_uniform-s)/norm(s);
        error2(ii,jj) = norm(sol_uniform-s)/norm(s);

    end

end

%%

switch G.type
    case 'sensor'
        figure(1)
        plot(Ns,mean(error1),'rx-',Ns,mean(error2),'bx-')
        legend('Adapted sampling','Random sampling')
        xlabel('Number of measurements')
        ylabel('Relative recovery error')
        if ~SAVE
            title('Sensor graph - Overlap level')
        end
        
        gsp_reset_seed(3)
        Nsf = 50;
        clear G
        G = gsp_sensor(N);
        G = gsp_estimate_lmax(G);
        g = gsp_design_heat(G,5);
        ntig = gsp_norm_tig(G,@(x) (g(x)).^2);
        
        % Perform an aditional experiment
        % Perform the sampling
        sample_non_uniform = randsampleWRW((1:N),Nsmax,ntig);    
        sample_uniform = randsample(N,Nsmax,0);

        s = gsp_filter_analysis(G,h(G),randn(G.N,1));

        % Create the masks
        mask_non_uniform = zeros(N,1);
        mask_non_uniform(sample_non_uniform(1:Nsf)) = 1;
        mask_uniform = zeros(N,1);
        mask_uniform(sample_uniform(1:Nsf)) = 1;

        % Input signals
        y1 = mask_non_uniform .* s;
        y2 = mask_uniform .* s;

        param = struct;
        param.verbose = 0;
        param.exact = 1;

        % Solve the problem
        sol_non_uniform = gsp_regression_tik(G ,mask_non_uniform, y1 , 0, param );
        sol_uniform = gsp_regression_tik(G ,mask_uniform, y2 , 0, param );
        
        gsp_plotfig('adapted_sampling_rel_error_sensor',paramplot)
        figure(2)
        gsp_plot_signal(G, 1./ntig)
        if ~SAVE
            title('Sensor graph - Uncertainty/ O_1')
        end
        gsp_plotfig('adapted_sampling_uncertainty_sensor',paramplot)
        
        figure(3)
        gsp_plot_signal(G, ntig)
        if ~SAVE
            title('Sensor graph - ||T_ig||')
        end
        gsp_plotfig('adapted_sampling_sparsity_sensor',paramplot)

        figure(4)
        gsp_plot_signal(G,s)
        caxis([min(s), max(s)]);
        if ~SAVE
            title('Original signal');
        end
        gsp_plotfig('single_experiment_adapted_ori',paramplot)

        figure(5)
        % y2(logical(1-mask_non_uniform)) = min(y2);
        gsp_plot_signal(G,mask_non_uniform)
        % y1(logical(1-mask_uniform)) = min(y1);
        % gsp_plot_signal(G,y1)
        if ~SAVE
            title('Mask adaptated sampling');
        end
        gsp_plotfig('single_experiment_adapted_mask',paramplot)


        figure(6)
        gsp_plot_signal(G,mask_uniform)
        if ~SAVE
            title('Mask adaptated sampling');
        end
        gsp_plotfig('single_experiment_uniform_mask',paramplot)

        figure(7)
        gsp_plot_signal(G,sol_uniform)
        caxis([min(s), max(s)]);
        if ~SAVE
            title('Uniform sampling reconstruction');
        end
        gsp_plotfig('single_experiment_adapted_rec_un',paramplot)

        figure(8)
        gsp_plot_signal(G,sol_non_uniform)
        caxis([min(s), max(s)]);
        if ~SAVE
            title('Adaptated sampling reconstruction');
        end
        gsp_plotfig('single_experiment_adapted_rec_opt',paramplot)

        
    case 'community'

        figure(1)
        plot(Ns,mean(error1),'rx-',Ns,mean(error2),'bx-')
        legend('Adapted sampling','Random sampling','Location','Best')
        xlabel('Number of measurements')
        ylabel('Relative recovery error')
        if ~SAVE
            title('Community graph - Error level')
        end
        gsp_plotfig('adapted_sampling_rel_error_cummunity',paramplot)

        figure(2)
        gsp_plot_signal(G,1./ntig)
        if ~SAVE
            title('Community graph - Uncertainty/O_1')
        end
        gsp_plotfig('adapted_sampling_uncertainty_cummunity',paramplot)
        
        figure(3)
        gsp_plot_signal(G,ntig)
        if ~SAVE
            title('Community graph - ||T_ig||')
        end
        gsp_plotfig('adapted_sampling_sparsity_cummunity',paramplot)


    otherwise
        
        figure(1)
        plot(Ns,mean(error1),'rx-',Ns,mean(error2),'bx-')
        legend('Adaptated sampling','Random sampling','Location','Best')
        xlabel('Number of measurements')
        ylabel('Relative recovery error')


        figure(2)
        gsp_plot_signal(G,1./ntig)
        title('Local uncertainty estimation')
end


%%
fprintf('Relative error adaptated sampling: %0.3f\n', mean(error1(:)));
fprintf('Relative error random sampling: %0.3f\n', mean(error2(:)));







