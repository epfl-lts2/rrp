function a=plot_time(g,varagin)

Lg = length(g);


if mod(Lg,2)
    a =-(Lg-1)/2:(Lg-1)/2;
else
    a = -Lg/2:(Lg/2-1);
    
end

if nargin>1
    plot(a,fftshift(g),varagin);
else
    plot(a,fftshift(g));
end

xlim([min(a),max(a)]);
xlabel('Time (in samples)')


setplottime(g)

end

