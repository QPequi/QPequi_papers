import numpy as np
import matplotlib.pyplot as plt
import math as mt
import os

var = np.linspace(0, mt.pi/2,100)
s_ub = np.zeros(len(var))
s_lb = np.zeros(len(var))
s = np.zeros(len(var))

for v in range(len(var)):  
    s_lb[v] = (8/mt.pi**2)*var[v]**2
    s_ub[v] = - np.log(1-(2/mt.pi)*var[v])        

plt.plot(var,s_ub,label='upper')
plt.plot(var,s_lb,label='lower')
# plt.ylabel(r'$\bar{\langle \Sigma \rangle}$', fontsize=14)
    # plt.xlabel('h', fontsize=14)
    # plt.axvline(0.5, color="red",linestyle='--')
plt.legend()

np.savetxt(os.path.join(os.getcwd(),'s_lb.txt'), s_lb)
np.savetxt(os.path.join(os.getcwd(),'s_ub.txt'), s_ub)

# aqui juntei os códigos

dim = 100
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

plt.plot(var,s_ub,label='upper')
plt.plot(var,s_lb,label='lower')
plt.plot(var,s_exact_min)



