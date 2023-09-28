import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.optimize import curve_fit as fit

N_l = np.array([8,10,12,14,16])
N_dmrg_l = np.array([8,10,12,14,16,20,24,28,32,36,40])
SS = 3 #int(sys.argv[1])
chi = 400

twist = []
twist_rho = []
E = []

for N in N_l:
    try:
        twist.append(np.real(np.fromfile(f"data/twist_mpo_aklt_S_{SS/2}_L_{N}.dat",complex)))
    except:
        twist.append(np.real(np.fromfile(f"data/twist_mpo_aklt_S_{SS/2}_L_{N}_Sz_0.dat",complex)))
for N in N_dmrg_l:
    E.append(np.fromfile(f"data/E_dmrg_mpo_S_{SS/2}_L_{N}_chi_{chi}.dat"))
    twist_rho.append(np.real(np.fromfile(f"data/twist_dmrg_mpo_S_{SS/2}_L_{N}_chi_{chi}_rho.dat",complex)))

twist = np.array(twist)
twist_rho = -np.log(np.array(twist_rho))
E = -np.log(-np.array(E))

i0 = -2
p1 = np.poly1d(np.polyfit(1/N_l[i0:],twist[i0:,1]-twist[i0:,0],1))
print(p1)
p2 = np.poly1d(np.polyfit(1/N_dmrg_l[:],twist_rho[:,1]-twist_rho[:,0],1))
print(p2)
p3 = np.poly1d(np.polyfit(1/N_dmrg_l[6:],E[6:,1]-E[6:,0],1))
print(p3)

def func(x,a,b,c):
    return a+b*x+c*x**3

popt1,pcov1 = fit(func,1/N_l,twist[:,1]-twist[:,0])
x1 = np.linspace(0,1.1/N_l[0],100)
y1 = [func(i,popt1[0],popt1[1],popt1[2]) for i in x1]

popt2,pcov2 = fit(func,1/N_dmrg_l,twist_rho[:,1]-twist_rho[:,0])
x2 = np.linspace(0,1.1/N_dmrg_l[0],100)
y2 = [func(i,popt2[0],popt2[1],popt2[2]) for i in x2]

#popt3,pcov3 = fit(func,1/N_dmrg_l[4:],E[4:,1]-E[4:,0])
#x3 = np.linspace(0,1.1/N_dmrg_l[0],100)
#y3 = [func(i,popt3[0],popt3[1],popt3[2]) for i in x3]

plt.plot(1/N_l,twist[:,1]-twist[:,0],
         'o',label="$\lambda_\mathrm{twist}-\lambda_0$")
#plt.plot([0,1.1/N_l[0]],p1([0,1.1/N_l[0]]),'--',color='grey')
plt.plot(x1,y1,':',color='grey')

plt.plot(1/N_dmrg_l,twist_rho[:,1]-twist_rho[:,0],
         's',label=r"$-\ln\langle\rho\rangle_\mathrm{twist}-\lambda_0$")
#plt.plot([0,1.1/N_dmrg_l[0]],p2([0,1.1/N_dmrg_l[0]]),'--',color='grey')
plt.plot(x2,y2,':',color='grey')

#plt.plot(1/N_dmrg_l,E[:,1]-E[:,0],'.')
#plt.plot([0,1.1/N_dmrg_l[0]],p3([0,1.1/N_dmrg_l[0]]),'--',color='grey')
#plt.plot(x3,y3,':',color='grey')

plt.xlim([0,1.1/N_l[0]])
plt.ylim([0,np.max(twist[:,1]-twist[:,0])*1.1])
#plt.xticks([0]+list(1/N_l),['$\infty$']+list(N_l),fontsize=20)
plt.xticks([0,1/40,1/24,1/16,1/12,1/10,1/8],['$\infty$',40,24,16,12,10,8],fontsize=20)
plt.yticks(fontsize=20)
plt.xlabel("$L$",fontsize=20)
#plt.ylabel("$\lambda_\mathrm{twist}-\lambda_0$",fontsize=20)
plt.ylabel("$\Delta\lambda$",fontsize=20)
plt.legend(fontsize=20,loc='lower right')
plt.tight_layout()
plt.savefig(f"plot/twist_mpo_aklt_S_{SS/2}_chi_{chi}.pdf",format='pdf')
plt.show()

