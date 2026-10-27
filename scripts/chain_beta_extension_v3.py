#!/usr/bin/env python3
"""Chain Theory beta-function extension v3: actual finite-dimensional tests.

No result here derives the Dirac operator or 4D spacetime from the source axiom.
Requires numpy, sympy; optional import of v2 helper in same directory.
"""
from __future__ import annotations
import itertools
import math
import numpy as np
import sympy as sp


def gamma_matrices():
    sx=np.array([[0,1],[1,0]], complex)
    sy=np.array([[0,-1j],[1j,0]], complex)
    sz=np.diag([1.,-1.]).astype(complex)
    gam=[np.kron(sx, s) for s in (sx,sy,sz)] + [np.kron(sy,np.eye(2))]
    g5=gam[0]@gam[1]@gam[2]@gam[3]
    assert np.linalg.norm(g5-g5.conj().T)<1e-14
    assert np.linalg.norm(g5@g5-np.eye(4))<1e-14
    for i,j in itertools.product(range(4),repeat=2):
        assert np.linalg.norm(gam[i]@gam[j]+gam[j]@gam[i]-(2 if i==j else 0)*np.eye(4))<1e-14
    return gam,g5


def choi_orbit_test():
    """Different unitaries give distinct pure maximally entangled states,
    so entanglement alone does not distinguish a kinetic update."""
    d=4
    U0=np.eye(d,dtype=complex)
    U1=np.diag([1,1,1,-1]).astype(complex)
    results=[]
    vecs=[]
    for U in (U0,U1):
        v=sum(np.kron(np.eye(d)[:,k],U[:,k]) for k in range(d))/math.sqrt(d)
        rho=np.outer(v,v.conj())
        reshaped=v.reshape((d,d))
        rho_left=reshaped@reshaped.conj().T
        rho_right=reshaped.conj().T@reshaped
        assert np.linalg.norm(rho_left-np.eye(d)/d)<1e-14
        assert np.linalg.norm(rho_right-np.eye(d)/d)<1e-14
        assert abs(np.trace(rho@rho)-1)<1e-14
        vecs.append(v)
        results.append(dict(purity=float(np.trace(rho_left@rho_left).real),entropy=math.log2(d)))
    overlap=float(abs(np.vdot(*vecs)))
    assert abs(overlap-0.5)<1e-14
    return results,overlap


def wilson_corner_branches():
    """Exact massless corners for D_w(p) at tuned bare masses and r=1, a=1."""
    branches={}
    corners=list(itertools.product((0,1),repeat=4))
    for k in range(5):
        # m_bare=-2*k; scalar mass at a corner = m_bare+2*popcount(corner).
        massless=[p for p in corners if -2*k+2*sum(p)==0]
        assert len(massless)==math.comb(4,k)
        chirality_signs={(-1)**sum(p) for p in massless}
        assert chirality_signs=={(-1)**k}
        branches[k]=dict(bare_mass=-2*k, count=len(massless), corner_parity=(-1)**k)
    return branches


def overlap_free_symbol(p,m0,gam,g5):
    """Neuberger overlap D=I+gamma5 sign(gamma5*(D_w-m0)), r=1, a=1."""
    W=sum(1-math.cos(x) for x in p)
    sin2=sum(math.sin(x)**2 for x in p)
    x=(W-m0)*np.eye(4,dtype=complex)+1j*sum((math.sin(x)*g for x,g in zip(p,gam)),np.zeros((4,4),complex))
    return np.eye(4)+x/math.sqrt((W-m0)**2+sin2)


def overlap_free_spectrum():
    gam,g5=gamma_matrices()
    counts={}
    for m0 in (-1.,1.,3.,5.,7.,9.):
        zeros=[]
        for p in itertools.product((0.,math.pi),repeat=4):
            D=overlap_free_symbol(p,m0,gam,g5)
            if np.linalg.norm(D)<1e-12:
                zeros.append(p)
        count=sum(math.comb(4,k) for k in range(5) if 2*k<m0)
        assert count==len(zeros)
        counts[m0]=count
    assert list(counts.values())==[0,1,5,11,15,16]
    # for an arbitrary smooth noncorner momentum, verify GW and circle.
    p=(0.14,0.27,-0.32,0.55)
    D=overlap_free_symbol(p,1.,gam,g5)
    gw=np.linalg.norm(g5@D+D@g5-D@g5@D)
    herm=np.linalg.norm(D.conj().T-g5@D@g5)
    circ=max(abs(abs(l-1)-1) for l in np.linalg.eigvals(D))
    assert max(gw,herm,circ)<1e-12
    return counts,dict(gw=gw,gamma5_hermiticity=herm,overlap_circle=circ)


def lattice_wilson_matrix(K,L=2,seed=20261024):
    """4D Wilson-Dirac with discrete Z_K valued classical gauge BACKGROUND.
    Link gauge degrees of freedom are NOT quantized in this test.
    """
    gam,g5=gamma_matrices()
    sites=list(itertools.product(range(L),repeat=4)); ns=len(sites)
    ix={x:i for i,x in enumerate(sites)}
    rng=np.random.default_rng(seed)
    omega=np.exp(2j*np.pi/K)
    def shift(s,mu,step=1):
        out=list(s);out[mu]=(out[mu]+step)%L;return tuple(out)
    links={(s,mu):omega**int(rng.integers(K)) for s in sites for mu in range(4)}
    g={s:omega**int(rng.integers(K)) for s in sites}
    links_g={(s,mu):g[s]*links[(s,mu)]*np.conj(g[shift(s,mu)]) for s in sites for mu in range(4)}
    def wilson(lnks):
        out=np.zeros((4*ns,4*ns),dtype=complex)
        for s,idx in ix.items():
            sl=slice(4*idx,4*(idx+1))
            out[sl,sl]+=(4.*np.eye(4))
            for mu in range(4):
                fwd=shift(s,mu); bwd=shift(s,mu,-1)
                out[sl,slice(4*ix[fwd],4*(ix[fwd]+1))] += -0.5*(np.eye(4)-gam[mu])*lnks[(s,mu)]
                out[sl,slice(4*ix[bwd],4*(ix[bwd]+1))] += -0.5*(np.eye(4)+gam[mu])*np.conj(lnks[(bwd,mu)])
        return out
    DW=wilson(links); DW_g=wilson(links_g)
    BigG5=np.kron(np.eye(ns),g5)
    BigG=np.kron(np.diag([g[s] for s in sites]),np.eye(4))
    relative_cov=np.linalg.norm(DW_g-BigG@DW@BigG.conj().T)/np.linalg.norm(DW)
    gamma5herm=np.linalg.norm(DW.conj().T-BigG5@DW@BigG5)/np.linalg.norm(DW)
    assert max(relative_cov,gamma5herm)<1e-12
    return DW,DW_g,BigG,BigG5,relative_cov,gamma5herm


def matrix_overlap(DW,m0,G5):
    H=G5@(DW-m0*np.eye(DW.shape[0]))
    herm=np.linalg.norm(H-H.conj().T)
    assert herm<1e-11
    vals,vecs=np.linalg.eigh(H)
    min_gap=float(min(abs(vals)))
    assert min_gap>1e-9, 'Kernel has zero mode, overlap sign undefined.'
    signH=(vecs*np.sign(vals))@vecs.conj().T
    return np.eye(DW.shape[0])+G5@signH,min_gap


def overlap_gauge_test(K=7):
    DW,DWg,G,G5,cov,herm=lattice_wilson_matrix(K)
    D,gap=matrix_overlap(DW,1.05,G5)
    Dg,gapg=matrix_overlap(DWg,1.05,G5)
    relcov=np.linalg.norm(Dg-G@D@G.conj().T)/np.linalg.norm(D)
    relgw=np.linalg.norm(G5@D+D@G5-D@G5@D)/np.linalg.norm(D)
    relg5=np.linalg.norm(D.conj().T-G5@D@G5)/np.linalg.norm(D)
    circle=max(abs(abs(v-1)-1) for v in np.linalg.eigvals(D))
    assert max(relcov,relgw,relg5,circle)<1e-11
    return dict(K=K, wilson_cov=cov, wilson_gamma5_herm=herm, kernel_gap=gap,
                overlap_cov=relcov,gw=relgw,overlap_gamma5_herm=relg5,circle=circle)




def one_dimensional_qubit_walk():
    """Construct Dirac-like propagation from additional unitary LOCAL shift+coin rules.

    U(k,m)=exp(-i*m*sigma_x) exp(-i*k*sigma_z).
    Its eigenphases satisfy cos(omega)=cos(m)cos(k), with
    omega=sqrt(k^2+m^2)+O((|k|+|m|)^3) near the origin.
    Only an illustrative added QCA model, not forced by Chain Theory.
    """
    sx=np.array([[0,1],[1,0]],complex)
    sz=np.diag([1.,-1.]).astype(complex)
    eye=np.eye(2,dtype=complex)
    rows=[]
    for eps in (0.2,0.1,0.05,0.025):
        m=0.7*eps; k=0.9*eps
        coin=math.cos(m)*eye-1j*math.sin(m)*sx
        shift=np.diag([np.exp(-1j*k),np.exp(1j*k)])
        U=coin@shift
        unitary=np.linalg.norm(U.conj().T@U-eye)
        omega=math.acos(math.cos(m)*math.cos(k))
        dirac=math.hypot(m,k)
        err=abs(omega-dirac)
        materr=np.linalg.norm(U-(eye-1j*(m*sx+k*sz)))
        assert unitary<1e-13
        rows.append(dict(epsilon=eps,phase=omega,dirac_energy=dirac,
                         phase_error=err,linear_operator_error=materr))
    assert rows[1]['phase_error'] < rows[0]['phase_error']/6
    assert rows[2]['phase_error'] < rows[1]['phase_error']/6
    assert rows[3]['phase_error'] < rows[2]['phase_error']/6
    return rows


def anomalies_and_vectorlike():
    Yq=sp.Rational(1,6)
    su2su2y_Q=3*Yq*sp.Rational(1,2)
    SU2_weak_doublets_Q=3
    assert su2su2y_Q==sp.Rational(1,4)
    assert SU2_weak_doublets_Q%2==1
    YL=sp.Rational(-1,2)
    assert su2su2y_Q+YL*sp.Rational(1,2)==0
    assert (SU2_weak_doublets_Q+1)%2==0
    # Vectorlike U1 pair, two left Weyl fields charges +q and -q
    q=sp.Integer(1)
    anomaly_cubic=q**3+(-q)**3
    anomaly_grav=q+(-q)
    sum_sq=q**2+(-q)**2
    db=sp.Rational(2,3)*sum_sq  # U(1) beta_e coefficient multiplying e^3/(16*pi^2)
    assert (anomaly_cubic,anomaly_grav,sum_sq,db)==(0,0,2,sp.Rational(4,3))
    # A FULL SM generation is a separate supplied matter roster, NOT a qudit derivation.
    # All fermions listed as left-handed Weyl fields, with conjugates for RH matter.
    generation=[('Q_L',3,2,sp.Rational(1,6)),
                ('u_c',3,1,sp.Rational(-2,3)),
                ('d_c',3,1,sp.Rational(1,3)),
                ('L_L',1,2,sp.Rational(-1,2)),
                ('e_c',1,1,sp.Integer(1))]
    grav=sum(nc*nw*y for _,nc,nw,y in generation)
    cubic=sum(nc*nw*y**3 for _,nc,nw,y in generation)
    su2sqY=(3*sp.Rational(1,6)+sp.Rational(-1,2))/2
    su3sqY=(2*sp.Rational(1,6)+sp.Rational(-2,3)+sp.Rational(1,3))/2
    su3cube=2-1-1
    n_doublets=3+1
    assert (grav,cubic,su2sqY,su3sqY,su3cube,n_doublets%2)==(0,0,0,0,0,0)
    return dict(mixed_su2su2y=su2su2y_Q,witten_odd=SU2_weak_doublets_Q%2,
                full_SM_generation_anomalies=(grav,cubic,su2sqY,su3sqY,su3cube),
                full_SM_witten_even=(n_doublets%2==0),
                vectorlike_y3=anomaly_cubic,vectorlike_grav=anomaly_grav,
                vectorlike_charge2=sum_sq, vectorlike_delta_b=db)


def main():
    print('STUDY: Chain Theory qudit beta extension v3; no original-axiom derivation claimed.')
    choi,ov=choi_orbit_test()
    print('Maximally entangled Choi alternatives:',choi,'mutual overlap=',ov)
    print('Wilson tuned critical branches (k=0..4):',wilson_corner_branches())
    cs,errs=overlap_free_spectrum()
    print('Free overlap m0 -> massless species:',cs,'smooth-momentum residuals:',errs)
    for K in (5,7):print('4D background-field overlap test:',overlap_gauge_test(K))
    print('Constructive 1D qubit-walk/Dirac small-k test:',one_dimensional_qubit_walk())
    print('Chiral anomaly and vectorlike counterexample:',anomalies_and_vectorlike())
    print('NEW v3 TESTS PASSED')

if __name__=='__main__': main()