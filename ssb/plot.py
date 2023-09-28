import numpy as np
import matplotlib.pyplot as plt

L = 8
E0 = np.fromfile(f"data/E_mg_L_{L}.dat")[:40]
ES = -np.log(np.fromfile(f"data/ES_ssb_3h_z_N_{2*L}_J1_1.0_J2_1.0_D_1.0_l_{L}.dat"))[-1:-40:-1]

plt.plot(E0,'o')
plt.xlabel("$n$",fontsize=30)
plt.ylabel("$E$",fontsize=30)
plt.xticks([])
plt.yticks([-3,-2.5,-2,-1.5],fontsize=30)
plt.tight_layout()
plt.savefig("plot/E.pdf",format='pdf')

plt.figure()
plt.plot(ES,'o')
plt.xlabel("$n$",fontsize=30)
plt.ylabel("$\lambda$",fontsize=30)
plt.xticks([])
plt.yticks([4.4,4.6,4.8,5.0],fontsize=30)
plt.tight_layout()
plt.savefig("plot/lambda.pdf",format='pdf')
plt.show()
