#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

modes=['GG','GS']
Ltexts=['5','85']
Ls=[5,8.5]
tmax=70
jmax=11
numIt=23000
outputEvery=350
dt=0.004
tIts=numIt//outputEvery
ks=['0.1', '0.5', '1.0']


ts=np.array([it*dt for it in range(numIt+1) if it%outputEvery==0])
last_t = np.sum(ts<tmax)

fig1=plt.figure()
ax1=fig1.add_subplot(111)

fig2=plt.figure()
ax2=fig2.add_subplot(111)

for mode in modes:
    for Ltext, L in zip(Ltexts, Ls):
        for n, k0 in enumerate(ks):
            sumChis=np.loadtxt(f"./{mode}/L={Ltext}/DATA_sumChiInfo_CN_k0_{k0}.txt")
            
            sumChis_perj={}
            for j in range(jmax+1):
                sumChis_perj[j] = sumChis[j::(jmax+1)]
                for k in range(j):
                    sumChis_perj[j][:,1] -= sumChis_perj[k][:,1]
            
            for j in range(jmax+1):
                ax2.plot(ts[:last_t], sumChis_perj[j][:last_t,1], '-', label=f"State j={j}", linewidth=1.0, markersize=2 )
            if(n==0):
                ax2.legend(bbox_to_anchor=(1.01, 1.05))
            ax2.set_xlabel('Time (a.u.)')
            ax2.set_ylabel('Population in j-th TSEig')
            
            fig2.savefig(f"Population_{mode}_L{L}_k0_{k0}.png",dpi=400,bbox_inches='tight')
            ax2.clear()
            
            
            transmissions=np.loadtxt(f"./{mode}/L={Ltext}/DATA_CN_areaProps_k={k0}.txt")
            ax1.plot(ts[:last_t], transmissions[:last_t], '-', label=f"Cranck Nicolson k0={k0}")
        
        ax1.set_xlabel('Time (a.u.)')
        ax1.set_ylabel('Transmission trhough slit')
        ax1.legend()
        ax1.grid(True)
        ax1.set_ylim(0,0.62)
        fig1.savefig(f"Transmission_{mode}_CN_L{L}.png", dpi=400)
        ax1.clear()

