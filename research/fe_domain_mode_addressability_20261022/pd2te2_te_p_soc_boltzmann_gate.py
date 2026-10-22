#!/usr/bin/env python3
"""BiDB Pd2Te2 id 980 (or Pd2Se2 id 975): geometric Te/Se p-orbital
Slater-Koster plus atomic L.S, 24-band spinful illustrative Boltzmann model.

This is NOT a fitted Pd2Te2 band structure, a source parametrization, or
an experimental spin susceptibility. Pd d orbitals are intentionally
omitted. All numerical amplitudes, SOC, T and chemical potential are
chosen model values, not BiDB material quantities.
Input: original local ASE SQLite via --db; the DB is NEVER uploaded.
S = [sum (-f') <S> v] [sum (-f') v v]^{-1} j (units omitted).
C3v-orbit k quadrature avoids isolated degeneracies and makes the
integration symmetric without projecting the physical response itself.
Run: OPENBLAS_NUM_THREADS=1 python pd2te2_te_p_soc_boltzmann_gate.py \
      --db /path/to/bidb.db --id 980 --grid 8 12 16
"""
import argparse, itertools, json, sqlite3
import numpy as np
from scipy.linalg import expm
from scipy.special import expit

SX=np.array([[0,1],[1,0]],complex)
SY=np.array([[0,-1j],[1j,0]],complex)
SZ=np.diag([1.,-1.])
I2=np.eye(2)
LX=np.array([[0,0,0],[0,0,-1j],[0,1j,0]])
LY=np.array([[0,0,1j],[0,0,0],[-1j,0,0]])
LZ=np.array([[0,-1j,0],[1j,0,0],[0,0,0]])
LS=.5*(np.kron(LX,SX)+np.kron(LY,SY)+np.kron(LZ,SZ))
C3=np.array([[-.5,-np.sqrt(3)/2,0],[np.sqrt(3)/2,-.5,0],[0,0,1.]])
MX=np.diag([-1.,1.,1.])
MZ=np.diag([1.,1.,-1.])
ORBIT=(np.eye(3),C3,C3@C3,MX,C3@MX,C3@C3@MX)

class Model:
 def __init__(self,db,db_id=980):
  if db_id not in (980,975):
   raise ValueError("the Te/Se atom ordering is established only for ids 980/975")
  conn=sqlite3.connect("file:"+db+"?mode=ro",uri=True)
  conn.execute("pragma query_only=ON")
  row=conn.execute("SELECT numbers,positions,cell FROM systems WHERE id=?",
                   (db_id,)).fetchone()
  conn.close()
  if row is None: raise ValueError("BiDB id not found")
  nums=np.frombuffer(row[0],np.int32)
  pos=np.frombuffer(row[1],np.float64).reshape(-1,3)
  A=np.frombuffer(row[2],np.float64).reshape(3,3).copy()
  if not np.all(nums[[0,1,4,5]]==nums[0]):
   raise ValueError("chalcogen site ordering failed")
  self.db_id=db_id
  self.element_Z=int(nums[0])
  self.raw_A=A.copy()
  # Exact-hex projection: original BiDB relaxed cell differs by ~0.003 A.
  # This is a controlled MODEL ASSUMPTION, not an untouched database cell.
  A[1,0]=-A[0,0]/2
  A[1,1]=np.sqrt(3)*A[0,0]/2
  self.A=A
  xyz=pos[[0,1,4,5]].copy()
  frac=xyz@np.linalg.inv(self.raw_A)
  frac[:,:2]=np.array([[0,0],[1/3,2/3],[1/3,2/3],[0,0]])
  xyz[:,:2]=frac[:,:2]@A[:2,:2]
  center=(xyz[:,2].min()+xyz[:,2].max())/2
  outer=(xyz[3,2]-xyz[0,2])/2
  inner=(xyz[2,2]-xyz[1,2])/2
  xyz[:,2]=center+np.array([-outer,-inner,inner,outer])
  self.xyz=xyz
  self.cell_projection_A=float(np.linalg.norm(A-self.raw_A))
  self.matrices={"C3":C3,"mx":MX,"mz":MZ}
  self.permutations={"C3":range(4),"mx":range(4),"mz":(3,2,1,0)}
  self.spinU={"C3":expm(-1j*np.pi/3*SZ),
              "mx":-1j*SX,"mz":-1j*SZ}
  self.O=[]
  for odd in (True,False):
   operator=np.zeros((24,24),complex)
   for i in range(4):
    sgn=(-1 if i<2 else +1) if odd else 1
    operator[6*i:6*i+6,6*i:6*i+6]=sgn*np.kron(np.eye(3),SZ/2)
   self.O.append(operator)
  # Illustrative, adjustable p-p Slater–Koster amplitudes; NOT a DFT fit.
  self.pp_sigma=1.
  self.pp_pi=-.23
  self.decay=.8
  self.distance_ref=3.
  self.cutoff=7.2
  self.inner_onsite=.25
  self.outer_onsite=-.12

 def hoppings(self,u):
  u=np.asarray(u,float)
  xyz=self.xyz.copy()
  xyz[:2,:2]-=u/2
  xyz[2:,:2]+=u/2
  bonds=[]
  for i,j in itertools.product(range(4),repeat=2):
   for n1,n2 in itertools.product(range(-3,4),repeat=2):
    if i==j and n1==n2==0: continue
    d=xyz[j]+n1*self.A[0]+n2*self.A[1]-xyz[i]
    r=np.linalg.norm(d)
    if r>self.cutoff or r<1e-6: continue
    n=d/r
    orbital=(self.pp_pi*np.eye(3)+(self.pp_sigma-self.pp_pi)
             *np.outer(n,n))*np.exp(-(r-self.distance_ref)/self.decay)
    bonds.append((i,j,d,orbital))
  return bonds

 def h_v(self,k,bonds,soc):
  h=np.zeros((12,12),complex)
  v=np.zeros((2,12,12),complex)
  for i,j,d,t in bonds:
   phase=np.exp(1j*k@d[:2])
   sl=np.s_[3*i:3*i+3,3*j:3*j+3]
   h[sl]+=phase*t
   for a in (0,1): v[a][sl]+=1j*d[a]*phase*t
  H=np.kron(h,I2)
  V=np.array([np.kron(v[a],I2) for a in (0,1)])
  for i in range(4):
   sl=np.s_[6*i:6*i+6,6*i:6*i+6]
   H[sl]+=soc*LS+(self.inner_onsite if i in (1,2)
                 else self.outer_onsite)*np.eye(6)
  if np.linalg.norm(H-H.conj().T)>1e-9:
   raise AssertionError("H not Hermitian")
  return H,V

 def symU(self,op):
  U=np.zeros((24,24),complex)
  for i,j in enumerate(self.permutations[op]):
   U[6*j:6*j+6,6*i:6*i+6]=np.kron(self.matrices[op],
                                     self.spinU[op])
  return U

 def audit(self,soc=.48):
  rng=np.random.default_rng(123)
  errors={}
  for name in ("C3","mx","mz"):
   R=self.matrices[name]
   U=self.symU(name)
   worst=0.
   for _ in range(3):
    k=rng.uniform(-.9,.9,2)
    u=rng.uniform(-.1,.1,2)
    k2=(R@np.r_[k,0])[:2]
    u2=(R@np.r_[u,0])[:2]*(-1 if name=="mz" else 1)
    H,_=self.h_v(k,self.hoppings(u),soc)
    H2,_=self.h_v(k2,self.hoppings(u2),soc)
    worst=max(worst,float(np.linalg.norm(H2-U@H@U.conj().T)
                          /np.linalg.norm(H)))
   errors[name]=worst
  k=np.array([.31,.45])
  u=np.array([.11,-.07])
  H,V=self.h_v(k,self.hoppings(u),soc)
  Hm,Vm=self.h_v(-k,self.hoppings(u),soc)
  Theta=np.kron(np.eye(12),1j*SY)
  errors["time_reversal"]=float(np.linalg.norm(
    Hm-Theta@H.conj()@Theta.conj().T)/np.linalg.norm(H))
  errors["TR_velocity"]=max(float(np.linalg.norm(
    Vm[a]+Theta@V[a].conj()@Theta.conj().T)) for a in range(2))
  if max(errors.values())>1e-10: raise AssertionError(errors)
  return errors

 def chi_current(self,u=(0,0),nk=16,mu=1.5,T=.13,soc=.48):
  bonds=self.hoppings(u)
  conduct=np.zeros((2,2))
  source=np.zeros((2,2))
  reciprocal=np.linalg.inv(self.A[:2,:2].T)
  for n1,n2 in itertools.product(range(nk),repeat=2):
   # Shifted mesh prevents sampling exact high-symmetry degenerate bands.
   # Six C3v orbits make the quadrature point-group symmetric for any u.
   frac=np.array([(n1+.313)/nk,(n2+.127)/nk])
   k0=2*np.pi*frac@reciprocal
   for G in ORBIT:
    k=(G@np.r_[k0,0.])[:2]
    H,V=self.h_v(k,bonds,soc)
    ev,Q=np.linalg.eigh(H)
    vel=np.real(np.array([np.einsum("in,ij,jn->n",
                    Q.conj(),V[a],Q) for a in (0,1)]))
    spins=np.real(np.array([np.einsum("in,ij,jn->n",
                      Q.conj(),op,Q) for op in self.O]))
    w=expit((ev-mu)/T)*expit((mu-ev)/T)/T
    conduct+=(vel*w)@vel.T
    source+=(spins*w)@vel.T
  # Rows: S_-^z, S_+^z; columns: j_x, j_y.
  return source@np.linalg.inv(conduct)

def main():
 p=argparse.ArgumentParser()
 p.add_argument("--db",required=True)
 p.add_argument("--id",type=int,choices=[980,975],default=980)
 p.add_argument("--slide",type=float,default=.1)
 p.add_argument("--mu",type=float,default=1.5)
 p.add_argument("--T",type=float,default=.13)
 p.add_argument("--soc",type=float,default=.48)
 p.add_argument("--grid",nargs="+",type=int,default=[8,12,16])
 a=p.parse_args()
 model=Model(a.db,a.id)
 print("id:",a.id,"cell_symmetrization_A:",model.cell_projection_A)
 print("H/TR audit:",json.dumps(model.audit(a.soc),indent=2))
 for nk in a.grid:
  for soc in (0.,a.soc):
   neg=model.chi_current((-a.slide,0),nk,a.mu,a.T,soc)
   zero=model.chi_current((0,0),nk,a.mu,a.T,soc)
   pos=model.chi_current((a.slide,0),nk,a.mu,a.T,soc)
   odd_error=float(np.max(np.abs(neg[0]+pos[0])))
   even_error=float(np.max(np.abs(neg[1]-pos[1])))
   assert max(odd_error,even_error)<1e-9
   if soc==0:assert max(np.max(np.abs(x))
                        for x in (neg,zero,pos))<1e-8
   print(json.dumps({"nk":nk,"soc":soc,
     "chi_minus_u":np.round(neg,10).tolist(),
     "chi_zero_u":np.round(zero,10).tolist(),
     "chi_plus_u":np.round(pos,10).tolist(),
     "odd_err":odd_error,"even_err":even_error}))

if __name__=="__main__":
 main()
