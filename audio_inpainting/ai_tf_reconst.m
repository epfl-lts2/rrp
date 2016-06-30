function srec = ai_tf_reconst(sinput,times,ratio,param)

anew = param.a*ratio;
win_lengthnew = param.win_length*ratio;
Mnew = win_lengthnew;
[A,B] = gabframebounds(firwin(param.win,win_lengthnew),anew,Mnew);

t1 = times(2,2) -times(1,1);
t2 = times(2,1) -times(1,2);
TX = max(t1,t2);

if param.keeplength
    TF1 = dgtreal(sinput(times(1,1)-2*win_lengthnew:times(1,1)+TX+2*win_lengthnew-1),...
        firwin(param.win,win_lengthnew),anew,Mnew,'timeinv');
    TF2 = dgtreal(sinput( times(1,2)-2*win_lengthnew:times(1,2)+TX+2*win_lengthnew-1),...
        firwin(param.win,win_lengthnew),anew,Mnew,'timeinv');

    bordersize = win_lengthnew/anew;
    TF1(:,2*bordersize+(1:size(TF2,2)-4*bordersize))=TF2(:,2*bordersize+1:end-2*bordersize);

    sout = idgtreal(TF1,firwin(param.win,win_lengthnew),anew,Mnew,TX+4*win_lengthnew,'timeinv')/B;
    sout = sout(win_lengthnew+1:end-win_lengthnew);

    srec = [sinput(1:times(1,1)-win_lengthnew-1);sout;...
            sinput(times(1,1)+TX+win_lengthnew:end)];

    %fs = 44100; sgram(srec( times(1,1)-2*fs : times(1,1)+(2*fs+TX)),'dynrange',100)
    %fs = 44100; sound(srec( times(1,1)-2*fs : times(1,1)+(2*fs+TX)),44100)
else
    Ls = length(sinput);
    ll = times(1,1)-2*win_lengthnew;
    lr = times(1,1)+2*win_lengthnew-1;            
    TF11 = dgtreal([zeros(1-ll,1);sinput(max(1,ll):min(Ls,lr));zeros(lr-Ls,1)],...
        firwin(param.win,win_lengthnew),anew,Mnew,'timeinv');
    
    Ls = length(sinput);
    ll = times(2,2)-2*win_lengthnew;
    lr = times(2,2)+2*win_lengthnew-1;       
    TF12 = dgtreal([zeros(1-ll,1);sinput(max(1,ll):min(Ls,lr));zeros(lr-Ls,1)],...
        firwin(param.win,win_lengthnew),anew,Mnew,'timeinv');
    
    Ls = length(sinput);
    ll = times(1,2)-2*win_lengthnew;
    lr = times(1,2)+2*win_lengthnew-1;       
    TF21= dgtreal([zeros(1-ll,1);sinput(max(1,ll):min(Ls,lr));zeros(lr-Ls,1)],...
        firwin(param.win,win_lengthnew),anew,Mnew,'timeinv');
    
    Ls = length(sinput);
    ll = times(2,1)-2*win_lengthnew;
    lr = times(2,1)+2*win_lengthnew-1;       
    TF22= dgtreal([zeros(1-ll,1);sinput(max(1,ll):min(Ls,lr));zeros(lr-Ls,1)],...
        firwin(param.win,win_lengthnew),anew,Mnew,'timeinv');

    bordersize = win_lengthnew/anew;
    TF31=[TF11(:,1:2*bordersize),TF21(:,2*bordersize+1:end)];
    TF32=[TF22(:,1:2*bordersize),TF12(:,2*bordersize+1:end)];  
    
    sout1 = idgtreal(TF31,firwin(param.win,win_lengthnew),anew,Mnew,4*win_lengthnew,'timeinv')/B;
    sout1 = sout1(win_lengthnew+1:end-win_lengthnew);
    
    sout2 = idgtreal(TF32,firwin(param.win,win_lengthnew),anew,Mnew,4*win_lengthnew,'timeinv')/B;
    sout2 = sout2(win_lengthnew+1:end-win_lengthnew);
    
    %stest1 = [sinput(times(1,1)+(-win_lengthnew:-1));sinput(times(1,2)+(0:win_lengthnew-1))];
    %stest2 = [sinput(times(2,1)+(-win_lengthnew:-1));sinput(times(2,2)+(0:win_lengthnew-1))];
    %
    %figure(1); plot(sout1-stest1);
    %figure(2); plot(sout2-stest2);
    
    srec = [sinput(1:times(1,1)-win_lengthnew-1);sout1;...
            sinput(times(1,2)+win_lengthnew:times(2,1)-win_lengthnew-1); 
            sout2;sinput(times(2,2)+win_lengthnew:end)];
    %stest2 = [sinput(1:times(1,1)-1);sinput(times(1,2):times(2,1)-1);...
    %    sinput(times(2,2):end)];
    %plot(srec-stest2);
        
    %fs = 44100; sgram(srec( times(1,1)-2*fs : times(1,1)+(2*fs+TX)),'dynrange',100);
    %fs = 44100; sound(srec( times(1,1)-2*fs : times(1,1)+(2*fs+TX)),44100);
end
if param.verbose
    disp(['Original length: ',num2str(length(sinput)),', new length: ',num2str(length(srec))]);
end

end