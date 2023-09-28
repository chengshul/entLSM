import numpy as np
import matplotlib.pyplot as plt
import sys

L_l = [12,14,16]
marker = ['o','x','+','*','<']
SS = 3 #int(sys.argv[1])

for i in range(len(L_l)):
    L = L_l[i]
    SA = np.concatenate(
            (np.fromfile(f"data/SA_mpo1_aklt_S_{SS/2}_L_{L}.dat"),
             np.fromfile(f"data/SL_mpo1_aklt_S_{SS/2}_L_{L}.dat")))
    print(len(SA)-L)
    assert(len(SA) == L)
    S2A = np.fromfile(f"data/S2A_mpo1_aklt_S_{SS/2}_L_{L}.dat")

    lA_l = np.arange(1,(L+2)//3)
    I = []
    for lA in lA_l:
        I.append(SA[lA*3-1]+SA[L-lA-1]-S2A[lA-1]-SA[L-1])
    I_fit = np.poly1d(np.polyfit(2*lA_l,np.log(I),1))
    print(1/I_fit[1])

    plt.figure(1)
    plt.semilogy(2*lA_l,I,marker[i],label=f"$L={L}$",markersize=10,markeredgewidth=2)
    if i == 2:
        plt.semilogy(2*lA_l,np.exp(I_fit(2*lA_l)),'--',color='grey')#,label=f"$\\xi_\\rho={np.round(1/np.abs(I_fit[1]),2)}$")

plt.figure(1)
plt.legend(fontsize=26)#,loc='lower left')
plt.xticks([2,4,6,8,10],fontsize=30)
#plt.yticks([1e-3,1e-4,1e-5,1e-6,1e-7],fontsize=40)
plt.yticks([1e-3,1e-5,1e-7],fontsize=30)
plt.xlabel("$|B|$",fontsize=30)
plt.ylabel("$I(A:C|B)$",fontsize=30)
plt.tight_layout()
plt.savefig(f"plot/MI_mpo_aklt_S_{SS/2}_L_{L}.pdf",format='pdf')

plt.show()
