"""Physical conservative coupler example, with unequal self-phase coefficients.
Run from any directory; writes the figure under josa-b-latex/plots/coupler.
Uses the audited elliptic reconstruction and leaves weierstrass_modified.py intact.
"""
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpmath import mp
from audit_coupler import make_solution

mp.dps=45
beta=[mp.mpf('.3'),mp.mpf('-.2')]
nonlinear=[[mp.mpf('1.2'),mp.mpf('.25')],[mp.mpf('.25'),mp.mpf('.7')]]
initial=[mp.mpc('.9','.2'),mp.mpc('.3','-.1')]
phase=mp.exp(-mp.j*mp.pi/4)
u0=[initial[0]*phase,mp.conj(initial[1])*phase]
v0=[-mp.conj(initial[0])/phase,initial[1]/phase]
aa=[-mp.j*beta[0],mp.j*beta[1]]
AA=[[mp.j*nonlinear[0][0],-mp.j*nonlinear[0][1]],[-mp.j*nonlinear[1][0],mp.j*nonlinear[1][1]]]
sol,R,_,_=make_solution(u0,v0,aa,AA)
z=np.linspace(0,8,501)
ratios=[R(mp.mpf(float(t))) for t in z]
angles=np.unwrap([float(mp.arg(r)) for r in ratios])
logs=[mp.log(abs(r/ratios[0]))+mp.j*mp.mpf(float(a-angles[0])) for r,a in zip(ratios,angles)]
fields=np.array([[complex(x) for x in sol(mp.mpf(float(t)),lr)] for t,lr in zip(z,logs)])
analytic=np.array([fields[:,0]/complex(phase),fields[:,1].conj()*complex(phase)]).T

def rhs(t,A):
 p=abs(A)**2
 return 1j*(np.array([float(b) for b in beta])+np.array(nonlinear,dtype=float)@p)*A+1j*A[::-1]
numeric=solve_ivp(rhs,(z[0],z[-1]),np.array(initial,dtype=complex),t_eval=z,method='DOP853',rtol=2e-13,atol=2e-14).y.T
error=float(np.max(abs(analytic-numeric)))
slice_error=float(max(np.max(abs(fields[:,2]+fields[:,0].conj())),np.max(abs(fields[:,3]-fields[:,1].conj()))))
power_error=float(np.max(abs(np.sum(abs(analytic)**2,axis=1)-sum(float(abs(x)**2) for x in initial))))
assert error<1e-10,error
assert slice_error<1e-11,slice_error
assert power_error<1e-11,power_error
# Match the existing paper figures: black curves, solid/dashed mode styles,
# circle/square numerical markers, dotted grey grids, and panel labels.
fig,axs=plt.subplots(1,2,figsize=(18,6))
linestyles=['-','--'];markers=['o','s']
indices=np.linspace(0,len(z)-1,20,dtype=int)
for j in range(2):
 axs[0].plot(z,abs(analytic[:,j])**2,color='black',linestyle=linestyles[j],linewidth=1.8,label=fr'Analytic $j={j+1}$')
 axs[0].plot(z[indices],abs(numeric[indices,j])**2,linestyle='None',marker=markers[j],markersize=7,color='black',label=fr'Numeric $j={j+1}$')
axs[0].set_ylabel(r'$|A_j|^2$',fontsize=18)
# Group analytic and numeric entries as in the other figures.
handles,labels=axs[0].get_legend_handles_labels();order=[0,2,1,3]
axs[0].legend([handles[i] for i in order],[labels[i] for i in order],loc='lower center',bbox_to_anchor=(.5,1.02),ncol=2,frameon=False,fontsize=14,handlelength=2.6,columnspacing=2)
phase_a=np.unwrap(np.angle(analytic[:,0]*analytic[:,1].conj()))
phase_n=np.unwrap(np.angle(numeric[:,0]*numeric[:,1].conj()))
axs[1].plot(z,phase_a,color='black',linestyle='-',linewidth=1.8,label='Analytic')
axs[1].plot(z[indices],phase_n[indices],linestyle='None',marker='o',markersize=7,color='black',label='Numeric')
axs[1].set_ylabel(r'$\arg(A_1A_2^*)$ (rad)',fontsize=18)
axs[1].legend(loc='lower center',bbox_to_anchor=(.5,1.02),ncol=2,frameon=False,fontsize=14,handlelength=2.6,columnspacing=2)
for ax,label in zip(axs,['A','B']):
 ax.set_xlabel(r'$z$',fontsize=18)
 ax.tick_params(labelsize=14)
 ax.grid(True,linestyle=':',color='gray',alpha=.5)
 ax.set_xlim(z[0],z[-1])
 ax.text(.02,1.02,label,transform=ax.transAxes,fontsize=20,fontweight='bold',va='bottom')
fig.subplots_adjust(left=.075,right=.98,bottom=.17,top=.82,wspace=.28)
out=Path(__file__).resolve().parents[1]/'josa-b-latex'/'plots'/'coupler'/'asymmetric_physical.pdf';out.parent.mkdir(parents=True,exist_ok=True);fig.savefig(out,bbox_inches="tight");plt.close(fig)
print(f'Asymmetric physical example: full-field max error {error:.6g}; slice residual {slice_error:.3g}; total-power residual {power_error:.3g}')
print(out)
