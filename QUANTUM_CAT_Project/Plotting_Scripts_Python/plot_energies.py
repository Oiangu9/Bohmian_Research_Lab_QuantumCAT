# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt
J10=5
J7=4
J5=3
K=10 #14
k0 = np.array([0.0,0.1,0.3,0.5,0.7,1.0,1.3,1.5,1.7,2.0,2.3,2.5,2.7,2.9])
E_GS = np.array([0.0374276,0.0424244,0.08239,0.162298,0.282076,0.536287,0.8795,1.15753,1.47461,2.02294,2.6576,3.1282,3.6362,4.18144])
E_GG = np.array([0.0590176,0.0640144,0.103983,0.183888,0.303666,0.557877,0.90112,1.17912,1.4962,2.04453,2.67923,3.14979,3.6578,4.20303])
E_eig_L5 = np.array([j**2*np.pi**2/25 for j in range(1,J5)])
E_eig_L7 = np.array([j**2*np.pi**2/49 for j in range(1,J7)])
E_eig_L85 = np.array([j**2*np.pi**2/(8.5**2) for j in range(1,J7)])
E_eig_L10 = np.array([j**2*np.pi**2/100 for j in range(1,J10)])


plt.plot(k0[:K],E_GS[:K], '-or', label="Gauss(x)*Phi^1(y)", linewidth=1.0, markersize=2.5)
plt.plot(k0[:K],E_GG[:K], '-ob', label="Gauss(x)*Gauss(y)", linewidth=1.0, markersize=2.5)
  

for j,E in enumerate(E_eig_L10):
    if j==1 or j==3:
        plt.plot(k0[K//2:K], [E for i in range(len(k0[K//2:K]))], '-', label=f"E^(j={j+1}) L=10", linewidth=1.0)        
    else:
        plt.plot(k0[:K], [E for i in range(len(k0[:K]))], '-', label=f"E^(j={j+1}) L=10", linewidth=1.0)
    
#for j,E in enumerate(E_eig_L7):
#    plt.plot(k0[:K], [E for i in range(len(k0[:K]))], ':', label=f"E^(j={j+1}) L=7", linewidth=1.0)
for j,E in enumerate(E_eig_L5):
    plt.plot(k0[:(K//2+1)], [E for i in range(len(k0[:(K//2+1)]))], '--', label=f"E^(j={j+1}) L=5", linewidth=1.0)
  

plt.legend(bbox_to_anchor=(1.01, 1.05))
plt.xlabel("k0 (a.u.) - Initial momentum")
plt.ylabel("Energy (a.u.)")
#plt.yscale('log')
#plt.show()
plt.savefig("Evsk_GS_GG_L5_.png", dpi=400,bbox_inches='tight')