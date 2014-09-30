function [  ] = ylatex( s )
%YLATEX This function put the string s as latex label on the y label


    latexint = get(gca,'DefaultTextInterpreter');
    set(gca,'DefaultTextInterpreter','tex');
    ylabel(s);
    set(gca,'DefaultTextInterpreter',latexint);

end

