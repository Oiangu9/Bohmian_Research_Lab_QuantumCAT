# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib


mode='GS'
L='5'
tmax=None
jmax=1
numIt=24000
outputEvery=130
dt=0.004
tIts=numIt//outputEvery
ks=['0.1', '0.3', '0.5', '0.7', '1.0']
Workers = ['W1','W2','W3','W4','W5']


ts=np.array([it*dt for it in range(numIt+1) if it%outputEvery==0])
if tmax==None:
    last_t = len(ts)
else:
    last_t = np.sum(ts<tmax)


matplotlib.rcParams.update({'font.size': 22})

fig1=plt.figure(figsize=(20,10)) # Transmissions
ax1=fig1.add_subplot(111)

fig2=plt.figure(figsize=(20,10)) # Absolute Errors
ax2=fig2.add_subplot(111)

fig3=plt.figure(figsize=(20,10)) # Relative Errors
ax3=fig3.add_subplot(111)

alg1_props={}
alg2_props={}
absolute_difs={}
relative_difs={}

ax1.set_xlabel("Time (a.u.)")
ax1.set_ylabel("Transmission through the slit")

ax2.set_xlabel("Time (a.u.)")

ax3.set_xlabel("Time (a.u.)")

try:
    os.mkdir(f"./{mode}_L_{L}_jmax_{jmax}/")
except:
    pass

for alg1, goodName1, shortName1 in zip(['CN_areaProps', 'XO_KA_trajProps'], ['Full Wavefunction 2D CN area', 'Truncated Born-Huang Approx (TBH)'], ['CN area', 'TBH']):
    for alg2, goodName2, shortName2 in zip(['XO_NoGJ_trajProps', 'CN_trajProps', 'XO_KA_trajProps'], ['Hermitian Approx (HA)', 'Full Wavefunction 2D CN traj. prop.', 'Truncated Born-Huang Approx (TBH)'], ['HA', 'CN trajs', 'TBH']):
        ax1.set_title(f"Transmission Comparison \n {goodName1} vs {goodName2}")
        ax2.set_title(f"Absolute difference between transmission (T) of \n {goodName1} vs {goodName2} ")
        ax3.set_title(f"% Relative transmission (T) difference of \n {goodName1} vs {goodName2} ")
        if alg1!=alg2:
            for k0, Wi in zip(ks, Workers):
                alg1_props[k0] = np.loadtxt(f"./{Wi}/DATA_{alg1}_k={k0}.txt")
                alg2_props[k0] = np.loadtxt(f"./{Wi}/DATA_{alg2}_k={k0}.txt")
                absolute_difs[k0] = np.abs(alg1_props[k0]-alg2_props[k0])
                relative_difs[k0] = 100*(alg1_props[k0]/(alg2_props[k0]+1e-5))
                ax2.set_ylabel(f"abs(T_({shortName1}) - T_({shortName2})")
                ax3.set_ylabel(f"100*(T_({shortName1}) / T({shortName2}))")
    
                ax1.plot(ts[:last_t], alg1_props[k0][:last_t], label=f"{shortName1} k0={k0}")
                ax1.plot(ts[:last_t], alg2_props[k0][:last_t], '--', label=f"{shortName2} k0={k0}")
                ax2.plot(ts[:last_t], absolute_difs[k0][:last_t], label=f"k0={k0}")
                ax3.plot(ts[:last_t], relative_difs[k0][:last_t], label=f"k0={k0}")
            ax1.legend()
            ax2.legend()
            ax3.legend()
            ax2.set_ylim((0,0.3))
            ax1.grid(True)
            ax2.grid(True)
            ax3.grid(True)
            fig1.savefig(f"./{mode}_L_{L}_jmax_{jmax}/transmission_{alg1}_vs_{alg2}", dpi=400)
            fig2.savefig(f"./{mode}_L_{L}_jmax_{jmax}/absolute_errors_{alg1}_vs_{alg2}", dpi=400)
            #fig3.savefig(f"./{mode}_L_{L}_jmax_{jmax}/relative_errors_{alg1}_vs_{alg2}", dpi=400)
    
            ax1.clear()
            ax2.clear()
            ax3.clear()
        
        
