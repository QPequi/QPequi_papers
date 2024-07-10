% xop   Gives the x operator as matrix. Usage: x(d),
%       where d is the size of the matrix.

function m=xop(d);

m=zeros(d,d);
for n=1:d-1
    m(n,n+1)=sqrt(n);    
end %for
m=(m+m')/sqrt(2);