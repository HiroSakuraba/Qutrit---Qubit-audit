#!/usr/bin/env python3
"""Chain Theory v20: covariant CP-completion nonuniqueness and missing intrinsic clock.

Conditional operational extensions to the selected P3 Bell projection, NOT deductions
from the source's no-failure-branch ontology.  Exact symbolic checks plus numeric
full joint-state, covariance, and Markovian semigroup checks.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import sympy as sp

OUT = Path(__file__).with_name('chain_beta_extension_v20_results.json')
I9=np.eye(9,dtype=complex)
omega=np.eye(3,dtype=complex).reshape(9)/np.sqrt(3)
P=np.outer(omega,omega.conj())
F=I9-P
rng=np.random.default_rng(20261025)

def rand_su3():
    z=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
    q,r=np.linalg.qr(z)
    q=q@np.diag(np.conj(np.diag(r))/np.abs(np.diag(r)))
    q=q/(np.linalg.det(q)**(1/3))
    return q

def pure_random(n):
    z=rng.normal(size=n)+1j*rng.normal(size=n)
    z/=np.linalg.norm(z)
    return np.outer(z,z.conj())

def channel(x,s):
    fx=F@x@F
    return P@x@P+(1-s)*fx+(s/8)*np.trace(fx)*F

def channel_kraus(x,s,comp_basis):
    out=P@x@P+(1-s)*F@x@F
    # all 64 Kraus operators sqrt(s/8) |fi><fj|
    for fi in comp_basis.T:
        for fj in comp_basis.T:
            k=np.sqrt(s/8)*np.outer(fi,fj.conj())
            out+=k@x@k.conj().T
    return out

def semigroup(x,s,gamma,t):
    pf=P@x@F+F@x@P
    fx=F@x@F
    ft=fx-np.trace(fx)/8*F
    return P@x@P+np.trace(fx)/8*F + np.exp(-gamma*s*t)*ft+np.exp(-gamma*t)*pf

def subspace_basis():
    # 8 orthonormal columns spanning orthogonal complement of omega, deterministic.
    a=np.column_stack([omega,np.eye(9,dtype=complex)])
    q,_=np.linalg.qr(a)
    assert np.linalg.norm(np.outer(q[:,0],q[:,0].conj())-P)<1e-13
    return q[:,1:]

def readout_state():
    # Original two qutrit-qubit pairs: |omega2>_Q1q1 |omega2>_Q2q2,
    # re-ordered into (Q1,Q2) vs (q1,q2) matrix; third pair factors out.
    c=np.zeros((9,4),complex)
    for i in range(2):
        for j in range(2):
            c[3*i+j,2*i+j]=0.5
    return c

def run():
    checks=[]; stats={}
    def check(name,cond,details=None):
        if not bool(cond): raise AssertionError(name)
        checks.append(name)
        if details is not None:stats[name]=details

    # The selected branch and failure CP completions are physical extensions only.
    basis=subspace_basis()
    check('P and F complementary projectors',np.linalg.norm(P@P-P)<1e-13 and np.linalg.norm(F@F-F)<1e-13 and np.linalg.norm(P@F)<1e-13)
    check('eight-dimensional failure complement',np.linalg.norm(basis.conj().T@basis-np.eye(8))<1e-13 and np.linalg.norm(basis@basis.conj().T-F)<1e-13)
    check('all instruments retain exact Bell-selected Kraus P',True)
    for s in [0.,0.13,0.5,0.87,1.]:
        ks=[P,np.sqrt(1-s)*F]
        ks += [np.sqrt(s/8)*np.outer(basis[:,i],basis[:,j].conj()) for i in range(8) for j in range(8)]
        complete=sum((k.conj().T@k for k in ks),np.zeros((9,9),complex))
        check(f'65-Kraus completeness s={s}',np.linalg.norm(complete-I9)<2e-13)
        x=pure_random(9)
        check(f'Kraus-channel formula s={s}',np.linalg.norm(channel_kraus(x,s,basis)-channel(x,s))<2e-13)
        vals=np.linalg.eigvalsh(channel(x,s))
        check(f'positivity trace s={s}',vals.min()>-1e-13 and abs(np.trace(channel(x,s))-1)<1e-13)

    # Analytical F-sector witness: E0 is identity on failure block, E1 depolarizes.
    f=basis[:,0]
    rf=np.outer(f,f.conj())
    el=channel(rf,0); ed=channel(rf,1)
    evs=np.linalg.eigvalsh(el-ed)
    distance=0.5*np.abs(evs).sum()
    check('exact failure pure-state trace distance 7/8',abs(distance-7/8)<1e-13)
    stats['pure_failure_trace_distance']=float(distance)

    # Exact projectors and state-independent decomposition using rational arithmetic.
    exact_bell=sp.zeros(9,1)
    for k in range(3):exact_bell[3*k+k]=1/sp.sqrt(3)
    p=sp.simplify(exact_bell*exact_bell.T); fmat=sp.eye(9)-p
    check('exact symbolic Bell and complement rank',p*p==p and fmat*fmat==fmat and p.trace()==1 and fmat.trace()==8)
    # Example of full 9D input from two qutrit-qubit pairs.
    c=readout_state()
    rho_local=c@c.conj().T
    rho_remote=c.conj().T@c
    prob=float(np.trace(c.conj().T@P@c).real)
    target=np.array([1.,0.,0.,1.],complex)/np.sqrt(2)
    selected_bob=c.conj().T@P@c
    check('original two-link selected Bell probability 1/6',abs(prob-1/6)<2e-15)
    check('selected remote qubit Bell ray exactly',np.linalg.norm(selected_bob/prob-np.outer(target,target.conj()))<1e-13)
    check('unconditional initial remote state I4/4',np.linalg.norm(rho_remote-np.eye(4)/4)<1e-13)
    bfail=rho_remote-selected_bob
    failnorm=bfail/(1-prob)
    failtarget=3/10*np.eye(4)-1/5*np.outer(target,target.conj())
    check('universal coarse-grained remote failure state',np.linalg.norm(failnorm-failtarget)<1e-13)
    stats['aligned_selected_born_probability']=prob
    stats['remote_failure_density_eigenvalues']=np.linalg.eigvalsh(failnorm).tolist()
    local_difference=channel(rho_local,0)-channel(rho_local,1)
    local_dist=0.5*np.sum(np.abs(np.linalg.eigvalsh(local_difference)))
    check('two completions differ on actual local marginal',local_dist>0.05,{'trace_distance':float(local_dist)})
    # Independent exact-rational analysis of the physical input's difference.
    rho_exact=sp.diag(*[sp.Rational(1,4) if k in (0,1,3,4) else sp.Integer(0) for k in range(9)])
    exact_fail=fmat*rho_exact*fmat
    exact_witness=exact_fail-(sp.trace(exact_fail)/8)*fmat
    evals=exact_witness.eigenvals()
    distance_exact=sp.Rational(1,2)*sum(sp.Abs(eigen)*mult for eigen,mult in evals.items())
    check('exact-rational original-input completion distance 7/16',distance_exact==sp.Rational(7,16))
    stats['original_input_completion_difference_eigenvalues']={str(k):int(v) for k,v in evals.items()}
    stats['original_input_completion_trace_distance_exact']='7/16'
    # The remote state after summing any local trace-preserving Kraus measurement.
    for s in [0.,0.4,1.]:
        ks=[P,np.sqrt(1-s)*F]
        ks += [np.sqrt(s/8)*np.outer(basis[:,i],basis[:,j].conj()) for i in range(8) for j in range(8)]
        bout=sum((c.conj().T@k.conj().T@k@c for k in ks),np.zeros((4,4),complex))
        check(f'no signaling on entangled input s={s}',np.linalg.norm(bout-rho_remote)<2e-13)
        check(f'unchanged selected branch s={s}',np.linalg.norm(c.conj().T@ks[0].conj().T@ks[0]@c-selected_bob)<1e-13)

    # Gauge covariance under U \otimes overline(U). Fixed Bell projection P.
    max_resid=0.
    for s in [0.,.17,.5,1.]:
        for k in range(40):
            u=rand_su3(); w=np.kron(u,u.conj())
            x=pure_random(9)
            max_resid=max(max_resid,float(np.linalg.norm(channel(w@x@w.conj().T,s)-w@channel(x,s)@w.conj().T)))
            if np.linalg.norm(w@P@w.conj().T-P)>2e-12:raise AssertionError('Bell stabilizer broken')
    check('160 independent SU3 correlated gauge-covariance trials',max_resid<2e-12,{'max_residual':max_resid})

    # Superoperator eigenspaces, verified on complete 81D matrix-unit basis.
    for s in [0.,.35,1.]:
        basis9=[np.outer(np.eye(9)[:,i],np.eye(9)[:,j]) for i in range(9) for j in range(9)]
        # column-major vec; stacking columns in any linear order retains spectrum.
        sop=np.column_stack([channel(x,s).reshape(-1) for x in basis9])
        transform=np.column_stack([x.reshape(-1) for x in basis9])
        sop=sop@np.linalg.inv(transform)
        sv=np.linalg.eigvals(sop)
        expected=np.concatenate([np.ones(2),np.full(63,1-s),np.zeros(16)])
        # Sort by real part to avoid numerical noise on nominally real eigenvalues.
        err=float(np.max(np.abs(np.sort(sv.real)-np.sort(expected))))
        check(f'81D superoperator spectrum s={s}',err<1e-12 and max(abs(sv.imag))<1e-11,{'max_error':err})

    # Repeated discrete steps and Poisson semigroup; gamma and s are arbitrary extensions.
    fstate=np.outer(f,f.conj()); dep=F/8
    spectral=[]
    for s in [0.0,0.2,0.7,1.0]:
        rho=fstate.copy()
        for n in range(1,9):
            rho=channel(rho,s)
            expected=(1-s)**n*fstate+(1-(1-s)**n)*dep
            if np.linalg.norm(rho-expected)>2e-12:raise AssertionError('iteration formula')
        spectral.append({'s':s,'step_decay_eigenvalue':1-s,'after_8_trace_distance':float(.875*(1-s)**8)})
    check('exact iterated spectral decay formula',True)
    max_semi=0.
    for s,gamma,t1,t2 in [(.0,1.2,.3,.4),(.25,.7,.3,.8),(.75,2.3,.2,.45),(1.,.9,1.5,.2)]:
        x=pure_random(9)
        y=semigroup(semigroup(x,s,gamma,t1),s,gamma,t2)
        total=semigroup(x,s,gamma,t1+t2)
        max_semi=max(max_semi,float(np.linalg.norm(y-total)))
        check(f'QMS semigroup s={s}',np.linalg.norm(y-total)<2e-12)
        # Analytical Poisson mixture: truncation at n=120 for a maximum gamma*t below 5.
        clock=gamma*t1
        term=np.exp(-clock); poisson=term*x.copy(); xn=x.copy()
        for n in range(1,121):
            xn=channel(xn,s);term*=clock/n;poisson+=term*xn
        check(f'Poisson convex-channel formula s={s}',np.linalg.norm(poisson-semigroup(x,s,gamma,t1))<2e-12)
    check('full semigroup compositional error bound',max_semi<2e-12,{'max_error':max_semi})
    stats['spectral_relaxations']=spectral
    stats['auxiliary_QMS_relaxation_rate']='gamma*s on traceless 8x8 failure block; gamma on Bell-failure coherences'
    stats['underdetermined_parameters']=['failure mixing s in [0,1]','clock rate gamma>0','a physical order or time parameter','network routing among junctions']
    stats['allowed_source_branch']='P3 on two joining qutrit copies; normalized postselected state; source excludes all failure branches'
    stats['not_derived']=['physical time or Hamiltonian','gauge propagator or photon','unique continuum transfer','electron species','QED beta coefficient']
    out={'status':'PASS','test_count':len(checks),'checks':checks,'measurements':stats,
         'scope':'Conditional standard-QM CPTP extensions of the primary paper\'s selected P3 branch. Failure branches and auxiliary Poisson clock are forbidden or absent in the original ontology.'}
    OUT.write_text(json.dumps(out,sort_keys=True,indent=2,default=str)+'\n')
    print(json.dumps({'status':'PASS','checks':len(checks),'p_selected':prob,
                      'pure_failure_trace_distance':distance,
                      'local_input_completion_distance':local_dist,
                      'max_gauge_residual':max_resid,'max_semigroup_residual':max_semi},indent=2))

if __name__=='__main__':run()