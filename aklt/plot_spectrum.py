import numpy as np
import matplotlib.pyplot as plt
import sys

L = 16 #int(sys.argv[1])
SS = 3 #int(sys.argv[2])

#ES0 = np.fromfile(f"data/ES_aklt_N_{2*L}_J1_{J1}_J2_{J2}_l_{L}.dat")
ES = np.fromfile(f"data/ES_1_mpo_aklt_S_{SS/2}_L_{L}_Sz_0.dat")
#ES = ES - ES[0]
k = np.fromfile(f"data/ES_k_mpo_aklt_S_{SS/2}_L_{L}_Sz_0.dat",int)
Stot = np.fromfile(f"data/ES_Stot_mpo_aklt_S_{SS/2}_L_{L}_Sz_0.dat",int)

#plt.plot(-np.log(ES0),'.')

#plt.figure()
for i in range(len(ES)):
    plt.plot(k[i],ES[i],'.',color=f"C{Stot[i]}",markersize=10)
    if k[i] == 0:
        plt.plot(L,ES[i],'.',color=f"C{Stot[i]}",markersize=10)

for i in range(max(Stot)+1):
    plt.plot([],[],'.',color=f"C{i}",label=f"$S={i}$",markersize=10)
plt.legend(loc="upper right",fontsize=18)

v = (ES[np.where(k==1)[0][0]] - ES[0]) / 2/np.pi*L

x = np.linspace(0,L,101)
y = np.abs(np.sin(x/L*2*np.pi)) * v + ES[0]
plt.plot(x,y,'--',color='grey')

plt.xlabel("$kL/2\pi$",fontsize=20)
plt.ylabel("$\lambda$",fontsize=20)
plt.xticks(range(0,L+1,4),fontsize=20)
#plt.yticks([6.2,6.6,7.0,7.4],fontsize=20)
plt.yticks([8.2,8.7,9.2],fontsize=20)
#plt.yticks([2,4,6,8],fontsize=20)
#plt.ylim([6.2,7.42])
plt.tight_layout()
plt.savefig(f"plot/ES_k_mpo_aklt_S_{SS/2}_L_{L}_Sz_0.pdf",format='pdf')
plt.show()

