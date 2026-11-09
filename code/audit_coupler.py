"""Independent checks for the coupler manuscript; run with SymPy, mpmath and SciPy.
The historical Weierstrass helper is left unchanged. Period orientation and
high-precision refinement are applied here, at the call site.
"""
import sympy as s
u=s.symbols('u1:3');v=s.symbols('v1:3');I=[u[j]*v[j] for j in range(2)]
a=s.symbols('a1:3');A=s.Matrix([[s.Symbol('a11'),s.Symbol('a12')],[s.Symbol('a12'),s.Symbol('a22')]])
H=-sum(a[j]*I[j] for j in range(2))-sum(A[j,k]*I[j]*I[k]/2 for j in range(2) for k in range(2))+s.prod(u)+s.prod(v)
f=[s.diff(H,x) for x in v]+[-s.diff(H,x) for x in u]
def along(expr):return s.expand(sum(s.diff(expr,x)*dx for x,dx in zip(u+v,f)))
assert along(H)==0 and along(I[0]-I[1])==0
assert any(s.expand(s.diff(H,u[0])*s.diff(I[0]-I[1],x)-s.diff(H,x)*s.diff(I[0]-I[1],u[0]))!=0 for x in u+v)
rho,g,lam=s.symbols('rho g lam');gs=[g,-g];Is=[gj-rho for gj in gs]
B1=-sum(a)-sum(A[j,k]*gs[k] for j in range(2) for k in range(2));B2=sum(A)/2
phip=[a[j]-sum(a)/2-gs[j]*sum(A)/2+sum((A[j,k]-sum(A[l,k] for l in range(2))/2)*Is[k] for k in range(2)) for j in range(2)]
assert s.expand(sum(phip))==0
for j in range(2):
 assert s.expand(-a[j]-sum(A[j,k]*Is[k] for k in range(2))+phip[j]-(B1-2*B2*Is[j])/2)==0
 r=(sum(A[j,k] for k in range(2))-sum(A)/4)/B2
 assert s.factor(r-s.Rational(1,2)-(-1)**j*(A[0,0]-A[1,1])/sum(A))==0
print('PASS: Hamiltonian, commuting independent integrals, gauge and branching exponents')
b0,b1,b2=s.symbols('b0 b1 b2');q=lambda x:b0+b1*x+b2*x*x
g1,g2=s.symbols('g1 g2')
L=(rho-g2)/(g1-g2)*q(g1)+(rho-g1)/(g2-g1)*q(g2)+(g1-rho)*(g2-rho)*b2
assert s.factor(L-q(rho))==0
w,wp=s.symbols('w wp');wj=s.symbols('w1:4');yj=s.symbols('y1:4')
L=sum(yj[j]*s.prod((w-wj[k])/(wj[j]-wj[k]) for k in range(3) if k!=j) for j in range(3))
Delta=(wj[0]-wj[1])*(wj[0]-wj[2])*(wj[1]-wj[2])
M=s.Matrix([[1,w,wp,6*w*w]]+[[1,wj[j],-yj[j],6*wj[j]**2] for j in range(3)])
assert s.factor(M.det()+6*Delta*(wp+L))==0
print('PASS: corrected interpolation and exact FS determinant factor')
# The full projection, including the two conjugate equations and energy.
T=s.symbols('T',nonzero=True);Uh,Vh=s.symbols('Uh Vh')
Q=lambda x:T+b1*(x-lam)+b2*(x*x-lam*lam)
D=2*T*(b1+2*b2*lam)-8*lam
chi=8*T/D;J3=D/(4*(lam-rho));Js=[D*(gj-rho)/(8*(gj-lam)*(lam-rho)) for gj in gs]
# Reduce rational expressions using the root condition T^2=4*(lam^2-g^2).
def onroot(expr):
 n,d=s.fraction(s.factor(expr));return s.factor(s.rem(n,s.Poly(T*T-4*(lam*lam-g*g),T),T))
for j in range(2):
 J=Js[j];ih=gs[j]-rho;theta=(-1)**j*(-b2*g+2*g/T)
 logu=-2/chi+chi*J/2+chi*J3/4-chi*(Vh*J3**2/(T*4*J)+Uh*J3/(2*T))+theta
 logv=2/chi-chi*J/2-chi*J3/4+chi*(Uh*J3**2/(T*4*J)+Vh*J3/(2*T))-theta
 assert onroot((logu-(b1-2*b2*ih)/2-Vh/ih).subs(Uh,Q(rho)-Vh))==0
 assert onroot((logv+(b1-2*b2*ih)/2+Uh/ih).subs(Vh,Q(rho)-Uh))==0
 assert onroot(2*J-J3-D/(4*(gs[j]-lam)))==0
Hb=-chi*(Uh+Vh)*J3**2/(4*T)+chi*sum(J*J for J in Js)/4+chi*J3**2/8-(sum(Js)+J3)/chi
E=-b2*D/8-2/chi+64*(D+16*lam)*lam/(chi**3*D**2)
assert onroot(Hb.subs(Uh,Q(rho)-Vh)-E)==0
print('PASS: full degenerate-FWM projection, weighted powers and Hamiltonian value')

from mpmath import mp
from weierstrass_modified import Weierstrass
mp.dps=55;we=Weierstrass()
def make_solution(u0,v0,aa,AA):
 ip=[u0[j]*v0[j] for j in range(2)];r0=-sum(ip)/2;gam=[ij+r0 for ij in ip]
 H=-sum(aa[j]*ip[j] for j in range(2))-sum(AA[j][k]*ip[j]*ip[k]/2 for j in range(2) for k in range(2))+mp.fprod(u0)+mp.fprod(v0)
 b0=H+sum(aa[j]*gam[j] for j in range(2))+sum(AA[j][k]*gam[j]*gam[k]/2 for j in range(2) for k in range(2))
 b1=-sum(aa)-sum(AA[j][k]*gam[k] for j in range(2) for k in range(2));b2=sum(sum(row) for row in AA)/2
 ds=[b0*b0-4*mp.fprod(gam),2*b0*b1,2*b0*b2+b1*b1-4,2*b1*b2,b2*b2]
 root=mp.polyroots(list(reversed(ds)),maxsteps=500,error=False)[0]
 D=ds[1]+2*ds[2]*root+3*ds[3]*root**2+4*ds[4]*root**3
 a=ds[2]/12+ds[3]*root/4+ds[4]*root*root/2
 inv=[ds[0]*ds[4]-ds[1]*ds[3]/4+ds[2]**2/12,ds[0]*ds[2]*ds[4]/6-ds[0]*ds[3]**2/16-ds[1]**2*ds[4]/16+ds[1]*ds[2]*ds[3]/48-ds[2]**3/216]
 om=list(we.omega_from_g(*inv)[:2])
 if mp.im(om[1]/om[0])<0:om[1]=-om[1]
 def exact_invariants(x,y):
  g2,g3=we.g_from_omega(x,y);return g2*(mp.mpf(4)/3)/mp.mpf(4/3),g3*(mp.mpf(8)/27)/mp.mpf(8/27)
 om=list(mp.findroot(lambda x,y:tuple(v-w for v,w in zip(exact_invariants(x,y),inv)),tuple(om),tol=mp.mpf('1e-48')))
 rhop=mp.fprod(u0)-mp.fprod(v0)
 z0=we.invwp(a+D/(4*(r0-root)),om,D*rhop/(4*(r0-root)**2))
 z1=we.invwp(a,om,-D*b2/4)
 nus=[we.invwp(a+D/(4*(gj-root)),om,-(b0+b1*gj+b2*gj*gj)*D/(4*(gj-root)**2)) for gj in gam]
 sig=lambda x:we.wsigma(x,om);zet=lambda x:we.wzeta(x,om);p=lambda x:we.wp(x,om)
 l1=[sum(AA[j])-sum(sum(row) for row in AA)/4 for j in range(2)]
 l0=[-aa[j]-gam[j]*sum(sum(row) for row in AA)/4-sum(AA[j][l]*gam[l] for l in range(2))+sum(AA[k][k]*gam[k] for k in range(2))/2+sum(aa)/2 for j in range(2)]
 rr=[l1[j]/b2-mp.mpf('.5') for j in range(2)]
 ee=[l0[j]+l1[j]*root-2*zet(z1)*l1[j]/b2-(zet(nus[j]-z1)+zet(nus[j]+z1))/2 for j in range(2)]
 R=lambda t:sig(t-z0+z1)/sig(t-z0-z1)
 # For local complex tests use the logarithm continued from zero on this chart.
 def sol(t,logratio=None):
  logratio=mp.log(R(t)/R(0)) if logratio is None else logratio
  us=[u0[j]*sig(t-z0+nus[j])/sig(-z0+nus[j])*sig(-z0-z1)/sig(t-z0-z1)*mp.exp(t*ee[j]+rr[j]*logratio) for j in range(2)]
  rh=root-D/(4*(a-p(t-z0)))
  return us+[(gam[j]-rh)/us[j] for j in range(2)]
 return sol,R,(aa,AA), dict(om=om,z0=z0,z1=z1,nus=nus,b2=b2,b1=b1,b0=b0,root=root,D=D,gamma=gam)

cases=[([mp.mpc('.487','-.828'),mp.mpc('-.971','-.865')],[mp.mpc('.942','.279'),mp.mpc('.386','-.523')],[mp.mpf('-.161'),mp.mpf('-.897')],[[mp.mpf('.583'),mp.mpf('-.0754')],[mp.mpf('-.0754'),mp.mpf('-.0525')]])]
# A conservative physical coupler with unequal propagation and self-phase terms.
phase=mp.exp(-mp.j*mp.pi/4);A0=[mp.mpc('.9','.2'),mp.mpc('.3','-.1')]
u0=[A0[0]*phase,mp.conj(A0[1])*phase];v0=[-mp.conj(A0[0])/phase,A0[1]/phase]
cases.append((u0,v0,[-mp.j*mp.mpf('.3'),mp.j*mp.mpf('-.2')],[[mp.j*mp.mpf('1.2'),-mp.j*mp.mpf('.25')],[-mp.j*mp.mpf('.25'),mp.j*mp.mpf('.7')]]))
worst=mp.mpf(0)
for u0,v0,aa,AA in cases:
 sol,_,_,_=make_solution(u0,v0,aa,AA)
 assert max(abs(x-y) for x,y in zip(sol(0),u0+v0))<mp.mpf('1e-42')
 for t in [mp.mpf(0),mp.mpf('.025'),mp.mpc('.04','.01')]:
  vals=sol(t);us,vs=vals[:2],vals[2:];ip=[us[j]*vs[j] for j in range(2)]
  rhs=[-(aa[j]+sum(AA[j][k]*ip[k] for k in range(2)))*us[j]+vs[1-j] for j in range(2)]+[(aa[j]+sum(AA[j][k]*ip[k] for k in range(2)))*vs[j]-us[1-j] for j in range(2)]
  worst=max(worst,max(abs(mp.diff(lambda t:sol(t)[j],t)-rhs[j]) for j in range(4)))
assert worst<mp.mpf('1e-40'),mp.nstr(worst)
print('PASS: all four reconstructed fields, generic complex and unequal physical cases; max ODE residual',mp.nstr(worst,6))
# Generic FS identity and signed lattice shifts, tested independently of the
# quartic-to-cubic reduction and of the symbolic notebook derivation.
om=[mp.mpf(1),mp.mpc('.27','1.13')]
sig=lambda x:we.wsigma(x,om);zet=lambda x:we.wzeta(x,om);wp=lambda x:we.wp(x,om)
assert abs(zet(om[0])*om[1]-zet(om[1])*om[0]-mp.j*mp.pi/2)<mp.mpf('1e-45')
z=mp.mpc('.31','.14');worst_shift=mp.mpf(0)
for m in range(-2,3):
 for n in range(-2,3):
  shift=2*m*om[1]+2*n*om[0];eta=m*zet(om[1])+n*zet(om[0])
  rhs=(-1)**(m*n+m+n)*sig(z)*mp.exp((shift+2*z)*eta)
  worst_shift=max(worst_shift,abs(sig(z+shift)/rhs-1))
assert worst_shift<mp.mpf('1e-42'),mp.nstr(worst_shift)
nu=[mp.mpc('.37','.11'),mp.mpc('-.24','.29'),mp.mpc('.18','-.32')]
C=12*sig(nu[0]-nu[1])*sig(nu[0]-nu[2])*sig(nu[1]-nu[2])/mp.fprod(sig(t)**4 for t in nu)
mat=mp.matrix([[1,wp(z),mp.diff(wp,z),mp.diff(wp,z,2)]]+[[1,wp(t),-mp.diff(wp,t),mp.diff(wp,t,2)] for t in nu])
left=C*mp.fprod(sig(z+t) for t in nu)*sig(z-sum(nu))/sig(z)**4
assert abs(left/mp.det(mat)-1)<mp.mpf('1e-42')
# The kernel normalisations in the discussion, including the reciprocal sign.
x=mp.mpc('.3','.2');a=mp.mpc('.24','-.1')
Phi=lambda x,a:sig(x+a)/(sig(x)*sig(a))*mp.exp(-x*zet(a))
assert abs(Phi(x,a)*Phi(x,-a)-(wp(x)-wp(a)))<mp.mpf('1e-42')
assert abs(mp.diff(lambda x:Phi(x,a),x,2)-(2*wp(x)+wp(a))*Phi(x,a))<mp.mpf('1e-42')
print('PASS: oriented Legendre relation, 25 sigma shifts, generic FS prefactor, and BA kernel identities')
# Reconstruct the six printed lifted meromorphic modes, including their phase
# normalisations, and check their equations without using the projection proof.
u0,v0,aa,AA=cases[0];sol,_,_,data=make_solution(u0,v0,aa,AA)
om=data['om'];z0=data['z0'];z1=data['z1'];nus=data['nus']+[z1];root=data['root'];D=data['D'];gam=data['gamma'];b2=data['b2']
Qroot=data['b0']+data['b1']*root+b2*root**2;chi=8*Qroot/D
cs=[mp.sqrt(2*(gj-root)) for gj in gam];sign=mp.fprod(cs)/Qroot
J30=D/(4*(root+sum(u0[j]*v0[j] for j in range(2))/2))
sig=lambda x:we.wsigma(x,om);zet=lambda x:we.wzeta(x,om);wp=lambda x:we.wp(x,om)
bu30=sig(-z0+z1)/(sig(z1)*sig(-z0));bv30=J30/bu30
bu0=[u0[j]*bv30/cs[j] for j in range(2)]+[bu30]
bv0=[v0[j]*bu30/cs[j] for j in range(2)]+[bv30]
sj=[1,1,2];J0=[bu0[j]*bv0[j] for j in range(3)];rbar=-sum(J0)/4;gj=[J0[j]+sj[j]*rbar for j in range(3)]
def lifted(t):
 us=[bu0[j]*sig(t-z0+nus[j])/sig(-z0+nus[j])*sig(-z0)/sig(t-z0)*mp.exp(t*(chi*gj[j]/sj[j]-zet(nus[j])+1/chi)) for j in range(3)]
 return us+[sj[j]*(wp(nus[j])-wp(t-z0))/(2*us[j]) for j in range(3)]
worst=mp.mpf(0)
for t in [mp.mpf(0),mp.mpf('.03'),mp.mpc('.05','.01')]:
 vals=lifted(t);us,vs=vals[:3],vals[3:];Is=[us[j]*vs[j] for j in range(3)]
 Um=us[0]*us[1]*us[2]**2;Vm=vs[0]*vs[1]*vs[2]**2
 rhs=[(-1/chi+chi*Is[j]/(2*sj[j]))*us[j]-sign*chi*sj[j]*Vm/(4*vs[j]) for j in range(3)]+[(1/chi-chi*Is[j]/(2*sj[j]))*vs[j]+sign*chi*sj[j]*Um/(4*us[j]) for j in range(3)]
 worst=max(worst,max(abs(mp.diff(lambda t:lifted(t)[j],t)-rhs[j]) for j in range(6)))
assert worst<mp.mpf('1e-40'),mp.nstr(worst)
print('PASS: all six printed lifted sigma modes and their initial data; max ODE residual',mp.nstr(worst,6))
# Constant Wronskians follow directly from the lifted vector field after
# reciprocal exponential gauging; this precedes the Lame identification.
bu=s.symbols('bu1:4');bv=s.symbols('bv1:4');ch,sg=s.symbols('chi sign',nonzero=True)
J=[bu[j]*bv[j] for j in range(3)];weights=[1,1,2]
Um=bu[0]*bu[1]*bu[2]**2;Vm=bv[0]*bv[1]*bv[2]**2
Hb=-sg*ch*(Um+Vm)/4+ch*(J[0]**2+J[1]**2)/4+ch*J[2]**2/8-sum(J)/ch
fu=[s.diff(Hb,x) for x in bv];fv=[-s.diff(Hb,x) for x in bu]
gj=[J[j]-weights[j]*sum(J)/4 for j in range(3)]
for j in range(3):
 cj=1/ch+ch*gj[j]/weights[j]
 W=fu[j]*bv[j]-bu[j]*fv[j]-2*cj*J[j]
 expected=weights[j]*Hb-4*gj[j]/ch-ch*gj[j]**2/weights[j]-ch*weights[j]*sum(gj[k]**2/weights[k] for k in range(3))/4
 assert s.expand(W-expected)==0
 assert s.expand(sum(s.diff(W,x)*dx for x,dx in zip(bu+bv,fu+fv)))==0
print('PASS: all three constant gauged Wronskians directly from the lifted equations')
Egamma=sum(gj[k]**2/weights[k] for k in range(3))
def flow(expr):return s.expand(sum(s.diff(expr,x)*dx for x,dx in zip(bu+bv,fu+fv)))
for j in range(3):
 cj=1/ch+ch*gj[j]/weights[j]
 potential=-sum(J)+4/ch**2+ch**2*Egamma/4-ch*Hb/2+2*gj[j]/weights[j]
 for expr in [flow(fu[j])-2*cj*fu[j]+cj**2*bu[j]-potential*bu[j],flow(fv[j])+2*cj*fv[j]+cj**2*bv[j]-potential*bv[j]]:
  assert s.expand(expr.subs(sg,1))==0
  assert s.expand(expr.subs(sg,-1))==0
print('PASS: shared generated polynomial potentials for all six gauged lifted modes, both coupling signs')
