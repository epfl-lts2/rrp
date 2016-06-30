function [x1,x2] = pdanalysis_dgt(y,param1,param2,gam,varargin)

% An implementation of Algorithm-2 for the analysis prior,
% from 'Algorithms for Audio Decomposition Using Mixed Norms'
%  
% Ilker Bayram
% ibayram@itu.edu.tr
% Istanbul Teknik Universitesi, 2012


win1 = param1.win;
Hop1 = param1.Hop;
N1 = length(win1);
lam1 = param1.lam;
W1 = param1.W;
L1 = param1.L;
SW1 = param1.SW;
SL1 = param1.SL;
CHA1 = N1;


win2 = param2.win;
Hop2 = param2.Hop;
N2 = length(win2);
lam2 = param2.lam;
W2 = param2.W;
L2 = param2.L;
SW2 = param2.SW;
SL2 = param2.SL;
CHA2 = N2;

MAX_ITER = 100;
if ~isempty(varargin),
    MAX_ITER = varargin{1};
end

% initialize

[w1,Ls1] = dgt(0*y,win1,Hop1,CHA1);
for r = 1:SW1:W1,
    for c = 1:SL1:L1,
        z{r,c} = w1;
    end
end

[w2,Ls2] = dgt(0*y,win2,Hop2,CHA2);
for r = 1:SW2:W2,
    for c = 1:SL2:L2,
        t{r,c} = 0*w2;
    end
end

sd1 = size(w1);

sd2 = size(w2);

x1 = zeros(size(y));
x2 = zeros(size(y));
xx1 = x1;
xx2 = x2;

wb = waitbar(0);
%figure;
for iter = 1:MAX_ITER,
    waitbar(iter/MAX_ITER,wb,strcat(num2str(round(100*iter/MAX_ITER)),'%'));    

    % primal step
    
    w1 = 0*z{1,1};
    for r = 1:SW1:W1,
        for c = 1:SL1:L1,
            w1 = w1 + z{r,c};
        end
    end
    v1 = idgt(w1,win1,Hop1,Ls1);    
    %v1 = v1(1:length(y));
    v1 = real(v1(1:length(y)));
    
    w2 = 0*t{1,1};
    for r = 1:SW2:W2,
        for c = 1:SL2:L2,
            w2 = w2 + t{r,c};
        end
    end
    v2 = idgt(w2,win2,Hop2,Ls2);
    %v2 = v2(1:length(y));
    v2 = real(v2(1:length(y)));
    
    u1 = x1 + gam*(y - lam1*v1);
    u2 = x2 + gam*(y - lam2*v2);
          
    x1 = ((gam + 1)*u1 - gam*u2)/(2*gam + 1);    
    x2 = ((gam + 1)*u2 - gam*u1)/(2*gam + 1);
    
    % dual step        
    w = dgt((2*x1 - xx1),win1,Hop1,CHA1);
    xx1 = x1;
    %subplot(1,2,1);imagesc(-(abs(w(end/2:-1:1,:))).^(1/4));colormap(gray);title('Transient Component');drawnow;
    W = W1;
    L = L1;
    SW = SW1;
    SL = SL1;
    lam = lam1;
    sd = sd1;
    for r = 1:SW:W,
        for c = 1:SL:L,
            
            u = (z{r,c} + lam*gam*w);
            % u is to be projected...
            
            % compute the energies of the blocks
            M1 = W * floor( (sd(1)- r + 1)/W);
            N1 = L * floor( (sd(2)- c + 1)/L);
            
            v = (abs(u(r:r+M1-1, c:c+N1-1))).^2;
            e = upfirdn(v,[0 ones(1,W)],1,W);
            e = upfirdn(e',[0 ones(1,L)],1,L);
            e = e(2:end,2:end)';
            e = sqrt(e);
            e = min(1./e,1);
            e = kron(e,ones(W,L));
            z{r,c}(r:r+M1-1,c:c+N1-1) = u(r:r+M1-1,c:c+N1-1).*e;                       
        end
    end
    
    w = dgt((2*x2 - xx2),win2,Hop2,CHA2);
    xx2 = x2;
    %subplot(1,2,2);imagesc(-(abs(w(end/2:-1:1,:))).^(1/4));colormap(gray);title('Tonal Component');drawnow;    
    W = W2;
    L = L2;
    SW = SW2;
    SL = SL2;
    lam = lam2;
    sd = sd2;
    for r = 1:SW:W,
        for c = 1:SL:L,
            
            u = (t{r,c} + lam*gam*w);
            % u is to be projected...
            
            % compute the energies of the blocks
            M1 = W * floor( (sd(1) - r + 1)/W);
            N1 = L * floor( (sd(2) - c + 1)/L);
            
            v = (abs(u(r:r+M1-1, c:c+N1-1))).^2;
            e = upfirdn(v,[0 ones(1,W)],1,W);
            e = upfirdn(e',[0 ones(1,L)],1,L);
            e = e(2:end,2:end)';
            e = sqrt(e);
            e = min(1./e,1);
            e = kron(e,ones(W,L));
            t{r,c}(r:r+M1-1,c:c+N1-1) = (u(r:r+M1-1,c:c+N1-1).*e);
        end
    end    
end
close(wb);