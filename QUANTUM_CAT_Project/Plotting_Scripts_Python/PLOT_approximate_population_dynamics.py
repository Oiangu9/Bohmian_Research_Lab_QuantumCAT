
import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib

mode='GS'
L='5'
tmax=None
J=1
jmax=7
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
times_computed=len(ts)

matplotlib.rcParams.update({'font.size': 22})

fig1=plt.figure(figsize=(20,10))

ax1 = []

ax1.append(fig1.add_subplot(111))

ax1[0].set_xlabel("Time (a.u.)")

ax1[0].grid(True)


pops_corr={}
pops_cwf={}
pops_leg={}

for k0, Wi in zip(ks, Workers):
    pops_corr[k0] = np.loadtxt(f"./{Wi}/DATA_population_correction.txt")
    pops_cwf[k0] = np.loadtxt(f"./{Wi}/DATA_populations_CWF.txt")
    pops_leg[k0] = np.loadtxt(f"./{Wi}/DATA_populations_legacy.txt")

numTrajs_pops_corr=[]
numTrajs_pops_leg=[]
numTrajs_pops_cwf=[]

for k0 in ks:
    numTrajs_pops_corr.append( pops_corr[k0].shape[0]//ts.shape[0])
    numTrajs_pops_leg.append( pops_leg[k0].shape[0]//ts.shape[0])
    numTrajs_pops_cwf.append( pops_cwf[k0].shape[0]//ts.shape[0])

print(numTrajs_pops_corr, numTrajs_pops_leg, numTrajs_pops_cwf)



# PLOT AVERAGE POPULATIONS of the ansatz GIVEN BY THE TRAJECTORIES
average_populations_corr={}
average_populations_leg={}

numTrajs=min(numTrajs_pops_corr)
for k0 in ks:
    average_populations_corr[k0]=np.zeros(( times_computed, next(iter(pops_corr.values())).shape[1]), dtype=np.float64)
    for i in range(numTrajs):
        beg=i*times_computed
        end=(i+1)*times_computed
        average_populations_corr[k0]+=pops_corr[k0][ beg:end,:]
    average_populations_corr[k0]/=numTrajs
        
numTrajs=min(numTrajs_pops_leg)
for k0 in ks:
    average_populations_leg[k0]=np.zeros(( times_computed, next(iter(pops_leg.values())).shape[1]), dtype=np.float64)
    for i in range(numTrajs):
        beg=i*times_computed
        end=(i+1)*times_computed
        average_populations_leg[k0]+=pops_leg[k0][ beg:end,:]
    average_populations_leg[k0]/=numTrajs

for k0 in ks:
    norm = np.sum(average_populations_corr[k0][:last_t,:], axis=1)
    for j in range(jmax+1):
        ax1[0].plot(ts[:last_t], average_populations_corr[k0][:last_t,j]/norm, '-', label=f"State j={j}", linewidth=1.0, markersize=2 )
    ax1[0].set_ylabel('Population in j-th TSEig')
    ax1[0].set_title(r"$\bf{"+mode+" L = "+L+" \ \ J = "+str(J)+" \ \ k0 ="+k0+"}$\nAverage accross trajectories\nof Corrected cwf product state populations", fontdict = {'fontsize' : 22})
    ax1[0].grid(True)
    ax1[0].legend()
    fig1.savefig(f"Average_Populations_Correction_{mode}_L{L}_k0_{k0}.png",dpi=400,bbox_inches='tight')
    ax1[0].clear()
    
for k0 in ks:
    #norm = np.sum(average_populations_leg[k0][:last_t,:], axis=1)
    norm=1
    for j in range(jmax+1):
        ax1[0].plot(ts[:last_t], average_populations_leg[k0][:last_t,j]/norm, '-', label=f"State j={j}", linewidth=1.0, markersize=2 )
    ax1[0].set_ylabel('Population in j-th TSEig')
    ax1[0].set_title(r"$\bf{"+mode+" L = "+L+" \ \  J = "+str(J)+" \ \ k0 ="+k0+"}$\nAverage accross trajectories\nof legacy cwf product state populations", fontdict = {'fontsize' : 22})
    ax1[0].grid(True)
    ax1[0].legend()
    fig1.savefig(f"Average_Populations_Legacy_{mode}_L{L}_k0_{k0}.png",dpi=400,bbox_inches='tight')
    ax1[0].clear()

# PLOT AVERAGE POPULATIONS computed with the approximated CWF definitions
numTrajs=min(numTrajs_pops_cwf)
current_time_diff_trajs=np.zeros((numTrajs, next(iter(pops_cwf.values())).shape[1]), dtype=np.float64)
populations_approx=np.zeros((times_computed, next(iter(pops_cwf.values())).shape[1]), dtype=np.float64)    

for k0 in ks:
    for t in range(last_t):
        for traj in range(numTrajs):
            current_time_diff_trajs[traj, :]= pops_cwf[k0][ traj*times_computed+t,:]
        current_time_diff_trajs=current_time_diff_trajs[np.argsort(current_time_diff_trajs[:,0]),:] # order them by x position
        for traj in range(numTrajs-1):
            dx=current_time_diff_trajs[traj+1,0]-current_time_diff_trajs[traj,0]
            populations_approx[t,:]+=dx*(current_time_diff_trajs[traj+1,:]+
                                                current_time_diff_trajs[traj,:])
        populations_approx[t,:]*=0.5
    norm = np.sum(populations_approx[:last_t, 1:], axis=1)
    for j in range(jmax+1):
        ax1[0].plot(ts[:last_t], populations_approx[:last_t,j+1]/norm, '-', label=f"State j={j}", linewidth=1.0, markersize=2 )
    ax1[0].set_ylabel('Population in j-th TSEig')
    ax1[0].set_title(r"$\bf{"+mode+" L = "+L+"\ \ J = "+str(J)+" \ \ k0 ="+k0+"}$\nIntegrated across different trajectory\ncwf-s in y used as a dynamic grid for full wf", fontdict = {'fontsize' : 22})
    ax1[0].grid(True)
    ax1[0].legend()
    fig1.savefig(f"Approx_Populations_using_CWF_definition_{mode}_L{L}_k0_{k0}.png",dpi=400,bbox_inches='tight')
    ax1[0].clear()
    

    









