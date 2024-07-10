% sud   SU(d) generators for a qudit
%   sud(d) gives the array of d^2-1 local orthogonal traceless observables
%   for a d-state system. obs=sud(d) gives back a d x d x (d^2-1)
%   dimensional matrix. The kth observable can be accessed as obs(:,:,k).
%   These observables fulfill trace(A_k*A_l)=0 for k>l and
%   trace(A^2)=2. 
%   (Similar to orthogobs.)

function obs=sud(d)

obs=zeros(d,d,d^2-1);

obs=orthogobs(d);
obs=obs(:,:,2:end);
% Make it traceless
% The frist d-1 vecror have nonzero trace
for n=1:d-1
    obs(:,:,n)=obs(:,:,n)-trace(obs(:,:,n))*eye(d,d)/d;
end %for
diagpart=zeros(d,d-1);
for n=1:d-1
    diagpart(:,n)=diag(obs(:,:,n));
end %for
% Gram-Schmidt orthogonalization
% The last column of Q is the constant column
[Q,R]=qr(diagpart);
for n=1:d-1
    obs(:,:,n)=diag(Q(:,n));
end %for
for n=1:d^2-1
    obs(:,:,n)=obs(:,:,n)/sqrt(trace(obs(:,:,n)^2))*sqrt(2);
end %for
