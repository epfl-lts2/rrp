% Graph comparizon

clear;
close all

% gsp_reset_seed(0)


alpha = 0.2;
k = 20;

%% Load the data

[x, y,xx, yy] = load_usps_full();
X = [x, xx];
Y = [y;yy];

Nso = size(X,2);
ny = 16;
nx = 16;
%%

data0 = zeros(nx*ny,nx*ny*Nso);
kk = 0;
X2 = reshape(X,nx,ny,[]);
for ii = 1:nx
    for jj =1:ny
        
        data0(:,(1:Nso)+kk*Nso) = reshape(circshift(circshift(X2,ii,1),jj,2),nx*ny,[]);
        kk = kk+1;
    end
    
end
data1 = X;
data3 = X(:,Y == 3 );
data7 = X(:,Y == 7 );
data9 = X(:,Y == 9 );



%% Covariance matrices for the full data

covM0 = gsp_stationarity_cov(data0);
covM1 = gsp_stationarity_cov(data1);
covM3 = gsp_stationarity_cov(data3);
covM7 = gsp_stationarity_cov(data7);
covM9 = gsp_stationarity_cov(data9);

% %% Covariance matrices for the sampled data
% sampled_stationarity_cov = @(X,Ns) gsp_stationarity_cov(X(:,randi(size(X,2),Ns,1)));
% covM0s = sampled_stationarity_cov(data0,Ns);
% covM1s = sampled_stationarity_cov(data1,Ns);
% covM2s = sampled_stationarity_cov(data2,Ns);


%% Graphs creation



Gg = gsp_2dgrid(nx,ny);
Gg2 = gsp_torus(nx,ny);
Gg3 = gsp_ring(nx*ny);
Gg = gsp_compute_fourier_basis(Gg);
Gg2 = gsp_compute_fourier_basis(Gg2);
Gg3 = gsp_compute_fourier_basis(Gg3);


param.k = k;
param.use_full = 1;
param.sigma = alpha*size(data0,2);
G0n = gsp_nn_graph(data0,param);

param.sigma = alpha*size(data1,2);
G1n = gsp_nn_graph(data1,param);

param.sigma = alpha*size(data3,2);
G3n = gsp_nn_graph(data3,param);

param.sigma = alpha*size(data7,2);
G7n = gsp_nn_graph(data7,param);

param.sigma = alpha*size(data9,2);
G9n = gsp_nn_graph(data9,param);


G0n = gsp_compute_fourier_basis(G0n);
G1n = gsp_compute_fourier_basis(G1n);

G3n = gsp_compute_fourier_basis(G3n);
G7n = gsp_compute_fourier_basis(G7n);
G9n = gsp_compute_fourier_basis(G9n);



%%
% 
% rg02 = gsp_stationarity_ratio(Gg2,covM0);
% rg03 = gsp_stationarity_ratio(Gg3,covM0);
% rg12 = gsp_stationarity_ratio(Gg2,covM1);
% rg22 = gsp_stationarity_ratio(Gg2,covM3);


rg0 = gsp_stationarity_ratio(Gg,covM0);
rg1 = gsp_stationarity_ratio(Gg,covM1);
rg3 = gsp_stationarity_ratio(Gg,covM3);
rg7 = gsp_stationarity_ratio(Gg,covM7);
rg9 = gsp_stationarity_ratio(Gg,covM9);


% r0 = gsp_stationarity_ratio(G0,covM0);
% r1 = gsp_stationarity_ratio(G1,covM1);
% r2 = gsp_stationarity_ratio(G2,covM2);

r0n = gsp_stationarity_ratio(G0n,covM0);
r1n = gsp_stationarity_ratio(G1n,covM1);
r3n = gsp_stationarity_ratio(G3n,covM3);
r7n = gsp_stationarity_ratio(G7n,covM7);
r9n = gsp_stationarity_ratio(G9n,covM9);




