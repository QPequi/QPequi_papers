

import numpy as np
import matplotlib.pyplot as plt


plt.rcParams["font.size"] = 17
plt.rcParams["text.usetex"] = True


vect = np.linspace(0, 30, 300)


##########################################
#
# Usual relative entropy
#
##########################################

# h=2,0
# tc_j50 = np.array([13,49,79,109,140,172,202,231,263,294])
# tm_j50 = np.array([5,21,36,51,66,81,96,111,126,141])

# tc_j100 = np.array([14,47,79,109,140,171,201,233,264,295])
# tm_j100 = np.array([10,26,42,58,74,90,106,121,137,153])

# tc_j200 = np.array([13,47,77,110,140,171,201,232,262,292])
# tm_j200 = np.array([11,27,43,59,75,91,107,122,138,154])

# tc_j300 = np.array([14,46,77,110,140,171,202,233,263,292])
# tm_j300 = np.array([11,27,43,59,77,91,107,123,139,154])

# h=1,6
# tc_j50 = np.array([17,60,98,135,173,209,251,288])
# tm_j50 = np.array([6,24,41,59,76,94,111,128])

# tc_j100 = np.array([14,60,98,136,174,213,248,286])
# tm_j100 = np.array([13,33,53,73,93,113,134,155])

# tc_j200 = np.array([16,59,98,136,174,211,250,290])
# tm_j200 = np.array([14,34,54,74,94,114,135,154])

# tc_j300 = np.array([16,57,97,136,174,212,249,288])
# tm_j300 = np.array([14,34,54,75,95,115,135,155])

# tc_j400 = np.array([18,57,96,136,174,212,251,288])
# tm_j400 = np.array([15,35,55,75,95,115,135,155])


#########################################
#
# Rényi relative entropy
#
#########################################

# h=1.4, \alpha = 0.8
# tc_j50 = np.array([17,60,98,135,173,209,251,288])
# tm_j50 = np.array([12,33,54,74,95,116,137,157])

# tc_j100 = np.array([14,60,98,136,174,213,248,286])
# tm_j100 = np.array([15,39,63,87,111,135,159,183])

# tc_j200 = np.array([16,59,98,136,174,211,250,290])
# tm_j200 = np.array([17,43,69,96,119,145,170,196])


# h=1.6, \alpha = 0.8
# tc_j50 = np.array([17,60,98,135,173,209,251,288])
# tm_j50 = np.array([7,25,44,61,80,98,106,116])

# tc_j100 = np.array([14,60,98,136,174,213,248,286])
# tm_j100 = np.array([12,33,55,76,98,119,140,162])

# tc_j200 = np.array([16,59,98,136,174,211,250,290])
# tm_j200 = np.array([15,37,61,84,105,128,148,169])



########################################
#
# Tsallis relative entropy
#
#########################################


# h=1.4, \alpha = 0.8
# tc_j50 = np.array([17,60,98,135,173,209,251,288])
# tm_j50 = np.array([12,33,54,74,95,116,137,157])

# tc_j100 = np.array([14,60,98,136,174,213,248,286])
# tm_j100 = np.array([14,38,62,87,111,135,159,183])

# tc_j200 = np.array([16,59,98,136,174,211,250,290])
# tm_j200 = np.array([14,39,64,90,116,142,168,194])

# h=1.6, \alpha = 0.8
# tc_j50 = np.array([17,60,98,135,173,209,251,288])
# tm_j50 = np.array([7,25,44,61,80,98,126,147])

# tc_j100 = np.array([14,60,98,136,174,213,248,286])
# tm_j100 = np.array([12,33,54,76,97,119,140,162])

# tc_j200 = np.array([16,59,98,136,174,211,250,290])
# tm_j200 = np.array([12,34,56,78,99,122,144,166])


plt.plot((vect[1:]+vect[0:-1])[tm_j50]/2, vect[tc_j50], linestyle='--', marker='o', label=r"$j=50$")
plt.plot((vect[1:]+vect[0:-1])[tm_j100]/2, vect[tc_j100], linestyle='--', marker='o', label=r"$j=100$")
plt.plot((vect[1:]+vect[0:-1])[tm_j200]/2, vect[tc_j200], linestyle='--', marker='o', label=r"$j=200$")
# plt.plot((vect[1:]+vect[0:-1])[tm_j300]/2, vect[tc_j300], linestyle='--', marker='o', label=r"$j=300$")
# plt.plot((vect[1:]+vect[0:-1])[tm_j400]/2, vect[tc_j400], linestyle='--', marker='o', label=r"$j=400$")
# plt.plot(np.linspace(0,20,50), 1.75*np.linspace(0,20,50), label=r"$t_c = 1.4 t_m$")
plt.legend(fontsize=12)
plt.xlabel(r"$t_m$")
plt.xlim([0, 19])
plt.ylabel(r"$t_c$")
plt.ylim([0,32])
# plt.title(r"Rényi entropy $(\alpha=0.8, h=1.6, \beta=1)$", fontsize=15)
plt.tight_layout()


# plt.savefig("tm_versus_tc_Renyi_h1.6_alfa0.8_beta1.pdf")