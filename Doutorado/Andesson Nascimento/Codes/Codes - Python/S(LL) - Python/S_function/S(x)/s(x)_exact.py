import numpy as np
import math as mt
import matplotlib.pyplot as plt
import os

dim = 10
xvec = np.linspace(0, 0.99,dim)
e = 10**(-12)

s_exact_min = np.zeros(dim)
s_exact = np.zeros((dim,dim))

for xx in range(len(xvec)):
    rvec = np.linspace(xvec[xx]+e,1-e,dim)
    for rr in range(len(rvec)):
        s_exact[xx,rr] = (rvec[rr] - xvec[xx])*mt.log((rvec[rr]-xvec[xx])/rvec[rr]) + (1-rvec[rr]+xvec[xx])*mt.log((1-rvec[rr]+xvec[xx])/(1-rvec[rr]))
       
for ii in range(len(xvec)):
    s_exact_min[ii] = np.min(s_exact[ii,]) 
    
np.savetxt(os.path.join(os.getcwd(),'s_exact.txt'), s_exact_min)

plt.plot(xvec,s_exact_min)

