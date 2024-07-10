% pop   Gives the p operator as matrix. Usage: p(d),
%       where d is the size of the matrix.

function m=pop(d);

m=zeros(d,d);
for n=1:d-1
    m(n,n+1)=sqrt(n);    
end %for
m=(m-m')/sqrt(2)/i;