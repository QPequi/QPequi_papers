function dydt = sch(t,y,alpha,j)

dim = 2*j+1;
lambda = (1/(4*j))*(tanh(alpha*(t-1)) + tanh(alpha));

dydt = zeros(dim,1);

% Equations of motion
dydt(1) = -1i*(3/2)*j.*y(1) - 1i*sqrt(j*(2*j-1))*(lambda/2).*y(3);

dydt(2) = -1i*((3*j-1)/2 + j - 1).*y(2)...
        -1i*sqrt((2*j+1)*(j-1)*3)*(lambda/2).*y(4);

ii = 3;
for n = (-j+2) : (j-2)

c1 = -1i*((j*(j+1) - n*n)/2 - n);
c2 = -1i*sqrt((j+n-1)*(j+n)*(j-n+1)*(j-n+2))*lambda/4;
c3 = -1i*sqrt((j+n+1)*(j-n)*(j-n-1)*(j+n+2))*lambda/4;

dydt(ii) = c1.*y(ii) + c2.*y(ii-2) + c3.*y(ii+2);    

ii = ii + 1;
end    

dydt(dim-1) = - 1i*((3*j-1)/2 - j + 1)*y(dim-1)...
            -  1i*sqrt((2*j-1)*(j-1)*3)*(lambda/2).*y(dim-3);

dydt(dim) = 1i*(j/2)*y(dim) - 1i*sqrt(j*(2*j-1))*(lambda/2).*y(dim-2);