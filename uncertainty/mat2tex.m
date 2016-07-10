function mat2tex(a,incenv,format)
% MAT2TEX(A,INCENV,FORMAT) - Nicki Holighaus 05.11.10
% returns the latex code for 
% the matrix in 'tabular' or 'array' format.
% A is the MATLAB/octave matrix to be converted.
% 'FORMAT' is an optional string parameter representing the 
% fprint format
% %d : Decimal notation (e.g. integer)
% %e : Exponential notation (% in 3.1415e+00)
% %E : Exponential notation (using an uppercase E as in 3.1415E+00)
% %f : Fixed-point notation
% %g : The more compact of %e or %f. Insignificant zeros do not print, sorry.
% %G : Same as %g, but using an uppercase E
% etc. See C reference for more information.
% Numbers (two number to be exact, with a decimal between them) can be injected between
% the '%' and the letter, as per a C reference of your choice.
% Default format is '%8.4f';
%
% INCENV is an additional parameter specifying the output format:
% 0       : array without equation environment
% 1       : array with equation environment (default)
% 'table' : table format 
%
% Examples: 
% 
% mat2tex(A,'table','%5.2f')
% mat2tex(A,0)
% mat2tex(A)


% This file is based on MATRIX2TEX.M, Copyright Joseph C. Slater, 2005

texstr='';
m=size(a,1);
n=size(a,2);
if nargin<=2
	format='%8.4f';
    ar = round(a);
    noint = ones(m,n)-(ar == a);
    if ( sum(noint(:)) == 0 )
      format='%d';
    end
    clear ar noint;
    if nargin==1
        incenv = 1;
    end
end

TAB = 0;

if (isstr(incenv) == 1 && isempty(strfind(incenv,'table')) == 0)
   TAB = 1;
end
    

eol='\\\\';
formataug=[format '&'];

if TAB == 1
    
    texstr = [texstr '\\begin{table}[thb] \n'];
    texstr = [texstr '\\begin{center}\\begin{tabular}{ |c|'];
    
    N=n;
    
    while N >= 1
        texstr = [texstr 'c'];
        N = N-1;
    end
    
    texstr = [texstr '|} \n \\hline \n'];
    
    N=n;
    
    while N >= 1
        texstr = [texstr ' & '];
        N = N-1;
    end
    
    texstr = [texstr eol ' \n \\hline \n'];
    
    for i=1:m
        texstr=[texstr  ' & $ '];
        for j=1:n
            
            if real(a(i,j))~=0|imag(a(i,j))==0
                texstr=[texstr sprintf(format,real(a(i,j)))];
            end
            if sign(imag(a(i,j)))==1&real(a(i,j))~=0&imag(a(i,j))~=0
                texstr=[texstr '+'];
            elseif sign(imag(a(i,j)))==-1&real(a(i,j))~=0&imag(a(i,j))~=0
                texstr=[texstr '-'];
            end
            if imag(a(i,j))~=0
                texstr=[texstr sprintf(format,abs(imag(a(i,j)))) 'i '];
            end
            if j==n
                texstr=[texstr  ' $ ' eol ' \n'];
            else
                texstr=[texstr ' $ & $ '];
            end
        end
    end
    
    texstr = [texstr ' \\hline \n \\end{tabular}\\end{center} \n \\caption{ } \n'];
    texstr = [texstr '\\end{table} \n'];
    
else
    
    if incenv == 1
        texstr = [texstr '\\begin{equation} \n'];
    end
    
    texstr = [texstr '\\left(\\begin{array}{'];

    N=n;
    while N >= 1
        texstr = [texstr 'c'];
        N = N-1;
    end
    
    texstr = [texstr '} \n'];
    
    for i=1:m
        for j=1:n
            if real(a(i,j))~=0|imag(a(i,j))==0
                texstr=[texstr sprintf(format,real(a(i,j)))];
            end
            if sign(imag(a(i,j)))==1&real(a(i,j))~=0&imag(a(i,j))~=0
                texstr=[texstr '+'];
            elseif sign(imag(a(i,j)))==-1&real(a(i,j))~=0&imag(a(i,j))~=0
                texstr=[texstr '-'];
            end
            if imag(a(i,j))~=0
                texstr=[texstr sprintf(format,abs(imag(a(i,j)))) 'i '];
            end
            if j==n&i==m
                texstr=[texstr  ' \n'];
            elseif j==n
                texstr=[texstr  eol ' \n'];
            else
                texstr=[texstr '&'];
            end
        end
    end
    
    texstr = [texstr '\\end{array}\\right) \n'];
    
    if incenv == 1
        texstr = [texstr '\\end{equation} \n'];
    end
    
end

%Octave/matlab differentiator
if exist('ver')==2
	sprintf(texstr)
else
	fprintf(texstr)
end