% Partitions. 2^k, with k = 1 to 7
veccell = [2 4 8 16 32 64 128];

% matrix with the largest eigenvalues for each partition
num = 10;
mat_eigs = zeros(num,length(veccell));

for count = 1 : length(veccell)
    
    fprintf('\n count = %d',count);
    
% Control variable for wrong points    
var_test = 0; 

numcell = veccell(count); % number of cells for dividing the phase space
stepphi = 2*pi/numcell;
steptheta = pi/numcell; 

mat_stocr = zeros(numcell*numcell,numcell*numcell); % stochastic matrix

if (numcell*numcell < num)
    aux_num = numcell*numcell;
else
    aux_num = num;
end

fid = fopen('data_final30.dat','r');

while ~feof(fid)

line = fgets(fid); %# read line by line
points = str2num(line);

theta0 = points(1);
phi0 = points(2);
thetaf = points(3);
phif = points(4);

% Compute the cell in which the point lies.

celula0 = [(floor(theta0/steptheta)+1) (floor((phi0+pi)/stepphi)+1)];
celulaf = [(floor(thetaf/steptheta)+1) (floor((phif+pi)/stepphi)+1)];

ind0 = (celula0(2) - 1)*numcell + celula0(1);
indf = (celulaf(1) - 1)*numcell + celulaf(2);

    if (ind0 <= 0)
       ind0 = 1;
       var_test = var_test + 1;
    end
    if (ind0 >= (numcell*numcell+1))
        ind0 = numcell*numcell;
        var_test = var_test + 1;
    end

    if (indf <= 0)
        indf = 1;
        var_test = var_test + 1;
    end
    if (indf >= (numcell*numcell+1))
        indf = numcell*numcell;
        var_test = var_test + 1;
    end

mat_stocr(indf,ind0) = mat_stocr(indf,ind0) + 1;

end % End of while loop
fclose(fid);

% normalization of the main matrix
vec = sum(mat_stocr);
for kk = 1 : (numcell*numcell)
    mat_stocr(:,kk) = mat_stocr(:,kk)./vec(kk);
end

% eigenvalues and eigenvectors of the stochastic matrix
s = svds(mat_stocr,num);

for hh = 1 : num
mat_eigs(hh,count) = s(hh);
end

fprintf('    vartest = %d',var_test);

clear fid s;
end % end of the refinment loop
fclose(fid1);

fid = fopen('eigs.dat','w'); 
for ii = 1 : num
    for kk = 1 : length(veccell)
    fprintf(fid,'%10.6f     ',mat_eigs(ii,kk));
    end
    fprintf(fid,'\n');
end
fclose(fid);