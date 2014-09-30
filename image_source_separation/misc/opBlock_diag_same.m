function op = opBlock_diag_same(J,op)
% Create a block Diagonal operator by the J times the repetition of the
% same operator op.
%
% This function is not optimal! 
%
% Autor: Nathana?l Perraudin 
% nathanael.perraudin@epfl.ch
%


info = op([],0);

m = info{1};            % Total number of rows
n = info{2};            % Total number of columns
 

op = @(x,mode) opBlockDiag_intrnl(m,n,op,x,mode,J);
end


function y = opBlockDiag_intrnl(m,n,op, x, mode,J)

    
   if mode==1 % Direct operator
   y  = zeros(m*J,1);
   kx = 0;
   ky = 0;
   for i=1:J

      
      
      y(ky+1:ky+m) = op(x(kx+1:kx+n),mode);
      kx              = kx + n;
      ky              = ky + m;
   end
   
   elseif mode==2 % inverse operator
   y  = zeros(n*J,1);
   kx = 0;
   ky = 0;
   for i=1:J

      
      
      y(ky+1:ky+n) = op(x(kx+1:kx+m),mode);
      kx              = kx + m;
      ky              = ky + n;
   end
   else
       error('Error here, unknow mode!\n 1 = operator, 2 = ajoint operator')
   end
   

end
