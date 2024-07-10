% sld   Symmetric Logarithmic Derivative. Usage: sld(rho,A), where
%       rho is a density matrix and A is an operator.
%       The form sld(rho,A,threshold) makes it possible
%       to define the threshold below which an eigenvalue is
%       consdiered zero. The default value is 1e-20.

function f=sld(rho,J,varargin)

% Make rho surely Hermitian
% A little bit of non-Hermitianicity can lead to
% non-orthogonal eigenvalues for eig
rho=(rho+rho')/2;

%Threshold=1e-10;
Threshold=1e-20;

if length(varargin)==1
   Threshold=varargin{1};
end %if

   [v,d]=eig(rho);
   [sy,sx]=size(rho);
   f=0;
   for n=1:sx
       for m=1:sx
           lambdan=d(n,n);
           lambdam=d(m,m);
           if abs(lambdam+lambdan)>Threshold;
               f=f+i*(lambdam-lambdan)/(lambdam+lambdan)*(ket(v(:,m))*bra(v(:,n)))*braket(v(:,m),J,v(:,n));
           end %if      
       end %for
   end %for

   f=f*2; % Because of missing factor of 2
   
end %if
    
