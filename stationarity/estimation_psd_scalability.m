%ESTIMATION_PSD_SCALABILITY Study the scalability of the PSD estimation method
%
%   Authors: Nathanael Perraudin and Pierre Vandergheynst
%
%   Date: January 2016
%
%   Paper: Stationary signal processing on graphs
%
%   Abstract of the paper
%   ---------------------
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
%   This experiment
%   ---------------
%
%   Figure 1 shows the time of PSD-estimation algorithm with respect of the
%   number of node. We use a $10$-nearest neighbors graph (sensor type) and
%   only $1$ signal with $M=30$ filters. For this experiment $K_2$ is set
%   to $10$ and the Chebysheff polynomial order is $30$. We average the
%   result over $10$ experiments.
%
%   .. figure::
%
%      Results
%
%      Computation time versus size of the graph. We use $M=30$ filters.
%      The algorithm scales linearly with the number of edges
%
%   Complexity analysis
%   -------------------
%
%   This method shows itself to be very powerful when the number of nodes
%   $N$ is greater than a few thousands. Diagonalizing the Laplacian
%   requires $\mathcal{O}(N^3)$ operations, while for a fixed error, the
%   approximation scales with the number of edges of the graph:
%   $\mathcal{O}(|E|)$, (which is proportional to $N$ in many graphs). In
%   fact, this estimation necessitates $ (K + K_2) M $ filtering operations
%   (with $M$ the number of points where the PSD). The final computation
%   cost of the method is thus $\mathcal{O}\left((K + K_2) M |E|\right)$.
%
%   References: perraudin2016stationary

clear
close all;

Ndata = 1;
Nfilt = 30;
Ntry = 10;
paramsensor.use_flann = 1;

N = 10000:10000:60000;

time = zeros(numel(N),Ntry);

for jj = 1:Ntry
    for ii = 1:numel(N) 

        % Create a graph
        x = rand(N(ii),3);
        G = gsp_nn_graph(x,paramsensor);
        G = gsp_estimate_lmax(G);


        % Take a band limited signal
        s0 = @(x) sin(x*4*pi/G.lmax).*(x<G.lmax/4);
        s = @(x) (s0(x)).^2;
        w = randn(N(ii),Ndata); % White noise
        x = gsp_filter(G,s0,w); %Original signal

        % PSD estimation
        fprintf('PSD estimation for N = %i  : ',N(ii))
        t = tic;
        param.Nfilt = Nfilt;
        psd = gsp_psd_estimation(G,x,param);
        time(ii,jj) = toc(t);
        fprintf(' %f s\n',time(ii));

    end
end


%% Plotting the result

figure(1)
h = errorbar(N,mean(time,2),var(time,0,2));
h(1).LineWidth = 1;

xlabel('Number of nodes')
ylabel('Time (s)')
title('Scalabilty')
axis([10000 60000 0 5])
paramplot.position = [100 100 300 200];
gsp_plotfig('psd_time',paramplot)
% figure(2)
% plot(N,time)

% %%
%     for ii = 1:numel(N) 
%         G = gsp_nd_sensor( N(ii),3,paramsensor);
%         Ne(ii) = G.Ne;
%     end
%     %%
% figure(3)
% plot(N,Ne)
