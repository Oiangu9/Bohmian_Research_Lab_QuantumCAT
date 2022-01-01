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

fig1=plt.figure(figsize=(40,20))

# Transmission TBH CN
ax1=fig1.add_subplot(221)

# Transmission HA CN
ax2=fig1.add_subplot(222)


# Absolute Errors TBH
ax3=fig1.add_subplot(223)

# Absolute Errors HA
ax4=fig1.add_subplot(224)

CN_props={}
absolute_difs_CN_TBH={}
TBH_props={}
HA_props={}
absolute_difs_CN_HA={}

ax1.set_xlabel("Time (a.u.)")
ax1.set_ylabel("Transmission through the slit")

ax2.set_xlabel("Time (a.u.)")

ax3.set_xlabel("Time (a.u.)")
ax4.set_xlabel("Time (a.u.)")


try:
    os.mkdir(f"./{mode}_L_{L}_jmax_{jmax}/")
except:
    pass

ax1.set_title(f"Transmission Comparison \n Full Wavefunction 2D CN area vs Truncated Born-Huang Approx (TBH)")
ax2.set_title(f"Transmission Comparison \n Full Wavefunction 2D CN area vs Hermitian Approx (HA)")
ax3.set_title(f"Absolute difference between transmission (T) of \n Full Wavefunction 2D CN area vs Truncated Born-Huang Approx (TBH) ")
ax4.set_title(f"Absolute difference between transmission (T) of \n Full Wavefunction 2D CN area vs Hermitian Approx (HA) ")
for k0, Wi in zip(ks, Workers):
    CN_props[k0] = np.loadtxt(f"./{Wi}/DATA_CN_areaProps_k={k0}.txt")
    TBH_props[k0] = np.loadtxt(f"./{Wi}/DATA_XO_KA_trajProps_k={k0}.txt")
    HA_props[k0] = np.loadtxt(f"./{Wi}/DATA_XO_NoGJ_trajProps_k={k0}.txt")

    absolute_difs_CN_TBH[k0] = np.abs(CN_props[k0]-TBH_props[k0])
    absolute_difs_CN_HA[k0] = np.abs(CN_props[k0]-HA_props[k0])
    ax1.plot(ts[:last_t], CN_props[k0][:last_t], label=f"CN area k0={k0}")
    ax1.plot(ts[:last_t], TBH_props[k0][:last_t], '--', label=f"TBH k0={k0}")
    ax2.plot(ts[:last_t], CN_props[k0][:last_t], label=f"CN area k0={k0}")
    ax2.plot(ts[:last_t], HA_props[k0][:last_t], '--', label=f"HA k0={k0}")
    ax3.plot(ts[:last_t], absolute_difs_CN_TBH[k0][:last_t], label=f"k0={k0}")
    ax4.plot(ts[:last_t], absolute_difs_CN_HA[k0][:last_t], label=f"k0={k0}")
    
ax3.set_ylabel(f"abs(T_(CN area) - T_(TBH)")
ax4.set_ylabel(f"abs(T_(CN area) - T_(HA)")



ax1.legend()
ax2.legend()
ax3.legend()
ax4.legend()
ax3.set_ylim((0,0.3))
ax4.set_ylim((0,0.3))

ax1.grid(True)
ax2.grid(True)
ax3.grid(True)
ax4.grid(True)
fig1.tight_layout(pad=3.0)

matplotlib.rcParams.update({'font.size': 45})

fig1.suptitle(f"{mode} J={jmax}",fontweight='bold')

fig1.savefig(f"./{mode}_L_{L}_jmax_{jmax}/transmission_{mode}_L_{L}_jmax_{jmax}.png", dpi=400)

        
        
