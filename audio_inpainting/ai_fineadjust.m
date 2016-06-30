function [shift,shiftr] = ai_fineadjust(shole,times,ratio,param)

switch param.finetune
    case {'none'}
        shift = 0;
        shiftr = shift;
    case {'wave'}
        anew = ratio*param.a;
        win_lengthnew = ratio*param.win_length;
        if param.keeplength
            exthole = shole( times(1,1)-win_lengthnew:times(2,2)+win_lengthnew);
            extcand = shole( times(1,2)-2*anew-win_lengthnew:times(2,1)+2*anew+win_lengthnew);
            temp = conv(extcand,flipud(exthole),'same');
            [~,shift] = max(temp(ceil(end/2)+(-2*anew+1:2*anew-1)));
            shift = shift-2*anew;
            shiftr = shift;
        else
            Ls = length(shole);
            ll = times(1,1)-win_lengthnew;
            lr = times(1,1)+win_lengthnew;            
            extholel = [zeros(1-ll,1);shole( max(1,ll):min(Ls,lr));zeros(lr-Ls,1)];
            
            ll = times(2,2)-win_lengthnew;
            lr = times(2,2)+win_lengthnew;
            extholer = [zeros(1-ll,1);shole( max(1,ll):min(Ls,lr));zeros(lr-Ls,1)];
            
            ll = times(1,2)-win_lengthnew-anew;
            lr = times(1,2)+win_lengthnew+anew;            
            extcandl = [zeros(1-ll,1);shole( max(1,ll):min(Ls,lr));zeros(lr-Ls,1)];
            
            ll = times(2,1)-win_lengthnew-anew;
            lr = times(2,1)+win_lengthnew+anew;   
            extcandr = [zeros(1-ll,1);shole( max(1,ll):min(Ls,lr));zeros(lr-Ls,1)];
            
            temp = conv(extcandl,flipud(extholel),'same');
            [~,shift] = max(temp(ceil(end/2)+(-anew+1:anew-1)));
            shift = shift-anew;
            
            temp = conv(extcandr,flipud(extholer),'same');
            [~,shiftr] = max(temp(ceil(end/2)+(-anew+1:anew-1)));
            shiftr = shiftr-anew;
        end
    case {'beat'}
        fs = param.fs;
        exthole = shole( [times(1,1)+(-fs:fs),times(2,2)+(-fs:fs)]);
        extcand = shole( [times(1,2)+(-fs:fs),times(2,1)+(-fs:fs)]);
        shift = rhythm_correction(extcand,exthole);
        shiftr = shift;
    otherwise
        error('Unknown method');
end

end

function shift = rhythm_correction(s1,s2)

% parameters for the transient component
%%% STFT parameters
N = 512;
Hop = 32;
win = bartlett(N);
win = NormalizeW(win,Hop);
param1.win = win;
param1.Hop = Hop;
%%% mixed norm parameters
param1.W = 16;
param1.L = 2;
param1.SW = 4;
param1.SL = 2;
param1.lam = 0.006;

% parameters for the tonal component
%%% STFT parameters
N = 1024;
Hop = 256;
win = hamming(N);
win = NormalizeW(win,Hop); % the actual window used (a Parseval frame)
param2.win = win;
param2.Hop = Hop;
%%% mixed norm parameters
param2.W = 2;
param2.L = 16;
param2.SW = 2;
param2.SL = 4;
param2.lam = 0.005;

% specify gamma
g1 = param1.lam*ceil(param1.W*param1.L/(param1.SW*param1.SL));
g2 = param2.lam*ceil(param2.W*param2.L/(param2.SW*param2.SL));
gam = min(1/g1,1/g2);

MAX_ITER = 30;

% compute transient harmonic decomposition
tic
[x11,x12] = PDAnalysis_dgt(s1,param1,param2,gam,MAX_ITER);
[x21,x22] = PDAnalysis_dgt(s2,param1,param2,gam,MAX_ITER);
toc

% Parameters for the rhythm detection
LKern = 2048;
kern = fftshift(pgauss(LKern,LKern/16));
kern = kern/norm(kern,1);
hop = param1.Hop;
smhop = 4;

ANA = dgt(x11,win,smhop,length(win)/16);

X = max(20*log10(abs(ANA))+80,0);
Y = X-circshift(X,[0,1]);
Z = sum(abs(Y),1);

smoothpeaksc = conv(Z,kern,'same');

p = peakpick(smoothpeaksc,1.5,100,3);

DIF = diff(p);
DIF = DIF < median(diff(p))/2;
FIRST = (p(1) == 1);
DIF = [FIRST,DIF];
P = p(DIF< 1);

deltac = zeros(size(ANA,2),1);
deltac(P) = 1;

%---------------------------------------------------

ANA = dgt(x21,win,smhop,length(win)/16);

X = max(20*log10(abs(ANA))+80,0);
Y = X-circshift(X,[0,1]);
Z = sum(abs(Y),1);

smoothpeaksh = conv(Z,kern,'same');

p = peakpick(smoothpeaksh,1.5,200,3);

DIF = diff(p);
DIF = DIF < median(diff(p))/2;
FIRST = (p(1) == 1);
DIF = [FIRST,DIF];
P = p(DIF< 1);

deltah = zeros(size(ANA,2),1);
deltah(P) = 1;

%----------------------------------------------------
% Compare the patterns

kern2 = fftshift(pgauss(LKern,LKern/256));

% Smoothen one of the delta trains
smdelta = conv(deltah,kern2,'same');

XY = conv(deltac,flipud(smdelta),'same');
XZ = XY(ceil(end/2)+(-2*hop/smhop+1:2*hop/smhop-1));
[j,i] = max(XZ);
shift = (i-2*hop/smhop)*smhop;

%figure(51);
%plot([deltac,circshift(deltah,[shift,0])]);
end