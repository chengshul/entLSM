import numpy as np
import matplotlib.pyplot as plt

L = 40

corr_H = np.abs(np.fromfile(f"data/mpo_corr_aklt_S_1.5_L_{L}.dat"))*3/4
corr_K = np.abs(np.fromfile(f"data/corr_dmrg_mpo_S_1.5_L_{L}_chi_400.dat"))*3/4

plt.semilogy(range(1,11),corr_H[:10],'o-',label=r"$|0\rangle_H$")
plt.semilogy(range(1,11),corr_K[:10],'s-',label=r"$|0\rangle_K$")

plt.xlabel("$|i-j|$",fontsize=30)
plt.ylabel(r"$|\langle\hat{\mathbf{S}}_{i,s}\cdot\hat{\mathbf{S}}_{j,s}\rangle|$",fontsize=30)
plt.xticks([2,4,6,8,10],fontsize=30)
plt.yticks([1e-1,1e-3,1e-5],fontsize=30)
plt.legend(fontsize=30)
plt.tight_layout()
plt.savefig(f"plot/corr_both_L_{L}.pdf",format='pdf')
plt.show()
