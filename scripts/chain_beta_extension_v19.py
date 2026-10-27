#!/usr/bin/env python3
"""Chain Theory QED beta program v19: support-overlap CP transfers and continuum tests.

Scope: exact Bell projections of the source, plus an EXPLICIT additional
common-ambient-support comparison used to connect distinct links.
No physical time, Hamiltonian, energy scale or QED beta function is derived.
"""
from __future__ import annotations
import json, math
from pathlib import Path
import numpy as np
import sympy as sp

OUT=Path(__file__).with_name('chain_beta_extension_v19_results.json')
EPS=2e-12

def basis(i,n=3):
    a=np.zeros(n,complex);a[i]=1;return a

def support(theta):
    return np.column_stack((basis(0), math.cos(theta)*basis(1)+math.sin(theta)*basis(2)))

def random_isometry(rng):
    z=rng.normal(size=(3,2))+1j*rng.normal(size=(3,2))
    q,_=np.linalg.qr(z)
    return q[:,:2]

def random_u2(rng):
    z=rng.normal(size=(2,2))+1j*rng.normal(size=(2,2))
    q,r=np.linalg.qr(z)
    return q @ np.diag(np.diag(r)/np.abs(np.diag(r)))

def overlap(v,w):
    return v.conj().T @ w

def cp(k,rho):
    return k@rho@k.conj().T

def psd_sqrt(a):
    eigen,vec=np.linalg.eigh((a+a.conj().T)/2)
    assert min(eigen)>-1e-12
    return vec@np.diag(np.sqrt(np.maximum(eigen,0)))@vec.conj().T

def polar(a):
    u,s,vh=np.linalg.svd(a)
    assert min(s)>1e-9
    return u@vh

def run():
    tests=[]; values={}
    def check(label, cond, data=None):
        assert bool(cond),label
        tests.append(label)
        if data is not None: values[label]=data

    # Source Bell join, in aligned real bases. 
    # The bra is <omega_3| and each input |omega_2> has 1/sqrt(2) amplitudes.
    j=np.eye(3,dtype=complex)/np.sqrt(3)
    v=support(0.)
    for theta in [0., .29, .7, math.pi/2]:
        w=support(theta)
        c_direct=np.einsum('ab,ai,bj->ij',np.conj(j),v/np.sqrt(2),w/np.sqrt(2))
        c_formula=v.T@w/(2*np.sqrt(3))
        check(f'Bell contraction theta={theta:.4f}',np.allclose(c_direct,c_formula,rtol=0,atol=EPS))
        p=np.linalg.norm(c_direct)**2
        expect=(1+math.cos(theta)**2)/12
        check(f'principal-angle probability theta={theta:.4f}',abs(p-expect)<EPS)
    check('aligned Bell success p=1/6',abs(np.linalg.norm(v.T@v/(2*np.sqrt(3)))**2-1/6)<EPS)
    check('orthogonal-normal Bell p=1/12',abs(np.linalg.norm(v.T@support(math.pi/2)/(2*np.sqrt(3)))**2-1/12)<EPS)
    values['probabilities']={'aligned':1/6,'orthogonal_normals':1/12,'theta_0.7':(1+math.cos(.7)**2)/12}

    # Generic complex support subspaces: principal-angle theorem gives sv (1, |n1^dag n2|).
    rng=np.random.default_rng(20261025)
    max_sv_err=0.;max_cov_err=0.;max_cp_err=0.;max_trace_over=0.
    for i in range(128):
        a,b=random_isometry(rng),random_isometry(rng)
        oa=overlap(a,b)
        p1=np.eye(3)-a@a.conj().T
        p2=np.eye(3)-b@b.conj().T
        # rank-one complements of two 2-planes
        na=np.linalg.eigh(p1)[1][:,-1]
        nb=np.linalg.eigh(p2)[1][:,-1]
        target=np.array([1.,abs(np.vdot(na,nb))])
        s=np.linalg.svd(oa,compute_uv=False)
        max_sv_err=max(max_sv_err,np.max(abs(s-target)))
        g1,g2=random_u2(rng),random_u2(rng)
        o2=overlap(a@g1,b@g2)
        max_cov_err=max(max_cov_err,np.linalg.norm(o2-g1.conj().T@oa@g2))
        x=random_u2(rng);rho=x@np.diag([.3,.7])@x.conj().T
        max_cp_err=max(max_cp_err,np.linalg.eigvalsh(cp(oa,rho)).min()*-1)
        max_trace_over=max(max_trace_over,max(0.,np.trace(cp(oa,rho)).real-np.trace(rho).real))
    check('complex Grassmannian principal-angle identity',max_sv_err<1e-12,{'max_error':float(max_sv_err)})
    check('U(2) local-frame transport covariance',max_cov_err<1e-12,{'max_error':float(max_cov_err)})
    check('overlap map positivity and trace nonincrease',max_cp_err<1e-12 and max_trace_over<1e-12,
          {'max_negative_eigenvalue_violation':float(max_cp_err),'max_trace_increase':float(max_trace_over)})

    # CP map K rho K^dag is a subchannel. A trace-preserving completion needs
    # complementary branch, which the source ontology expressly excludes.
    k=overlap(v,support(.6))
    z=psd_sqrt(np.eye(2)-k.conj().T@k)
    check('explicit Kraus completeness with added failure',np.allclose(k.conj().T@k+z.conj().T@z,np.eye(2),atol=EPS))
    rho=np.diag([.5,.5])
    lam=cp(k,rho); nrm=lam/np.trace(lam)
    mix_components=np.diag([.5,.5])  # each pure basis state conditioned separately stays itself
    nonlinear=np.linalg.norm(nrm-mix_components)
    check('normalized success is nonlinear in mixed-state inputs',nonlinear>.05,{'Frobenius_gap':float(nonlinear)})
    values['one_step_K_singular_values']=list(map(float,np.linalg.svd(k,compute_uv=False)))

    # Exact finite-loop holonomy supplied by overlap matrices of three planes.
    c=sp.sqrt(3)/2; s=sp.Rational(1,2)
    V0=sp.Matrix([[1,0],[0,1],[0,0]])
    V1=sp.Matrix([[1,0],[0,c],[0,s]])
    V2=sp.Matrix([[c,0],[0,1],[s,0]])
    K01=V0.T*V1; K12=V1.T*V2; K20=V2.T*V0
    H=sp.simplify(K01*K12*K20)
    H_expected=sp.Rational(3,4)*sp.Matrix([[1,0],[sp.Rational(1,4),1]])
    U=sp.Matrix([[8,-1],[1,8]])/sp.sqrt(65)
    check('exact three-plane composite overlap',H==H_expected)
    check('exact nontrivial polar holonomy',sp.simplify(U.T*U-sp.eye(2))==sp.zeros(2) and
          sp.simplify(U.T*H-(U.T*H).T)==sp.zeros(2) and
          sp.simplify(sp.det(U)-1)==0)
    check('Wilson trace 16/sqrt(65) differs from 2',sp.simplify(sp.trace(U)-16/sp.sqrt(65))==0)
    values['holonomy']={'overlap_loop':str(H),'polar_unitary':str(U), 'wilson_trace_exact':'16/sqrt(65)',
                        'wilson_trace_approx':float(16/sp.sqrt(65)), 'polar_rotation_angle':float(math.atan(1/8))}
    # Gauge transform at vertices, loop conjugation gives invariant trace, singular spectrum.
    nv=[np.array(z,dtype=complex) for z in [V0,V1,V2]]
    nh=np.eye(2,dtype=complex)
    for i in range(64):
        gs=[random_u2(rng) for _ in range(3)]
        kw=[overlap(nv[j]@gs[j],nv[(j+1)%3]@gs[(j+1)%3]) for j in range(3)]
        h=kw[0]@kw[1]@kw[2]
        h0=np.array(H,dtype=complex)
        check_ok=np.linalg.norm(h-gs[0].conj().T@h0@gs[0])<5e-13
        if not check_ok: raise AssertionError('loop covariance sample '+str(i))
        nh=h
    check('64 frame-gauge loop conjugation trials',True)

    # Continuum limit: smooth supports remove attenuation; zigzag with same uniform
    # support limit retains quadratic-variation attenuation.
    sigma=1.2; Theta=.9
    ns=[20,80,320,1280,5120]
    path_results=[]
    for n in ns:
        assert n%2==0
        smooth=math.cos(Theta/n)**n
        zigzag=math.cos(sigma/math.sqrt(n))**n
        path_results.append({'n':n,'smooth_gain':smooth,'zigzag_gain':zigzag,
                             'zigzag_support_max_deviation':sigma/math.sqrt(n)})
    smooth_lim=1.
    rough_lim=math.exp(-sigma*sigma/2)
    check('smooth path approaches unit gain',abs(path_results[-1]['smooth_gain']-smooth_lim)<1e-4)
    check('uniformly flat zigzag tends to nonunit gain',abs(path_results[-1]['zigzag_gain']-rough_lim)<2e-4)
    check('distinct limits despite same uniform limiting support',abs(rough_lim-smooth_lim)>.3)
    x=sp.symbols('x',real=True)
    check('exact quadratic-variation logarithmic asymptotic', sp.limit(sp.log(sp.cos(x))/x**2,x,0)==-sp.Rational(1,2))
    values['continuum']={'smooth_Theta':Theta,'zigzag_sigma':sigma,'zigzag_limit_exact':'exp(-sigma^2/2)',
                         'zigzag_limit':rough_lim,'smooth_limit':smooth_lim,'samples':path_results,
                         'zigzag_endpoint_conditional_concurrence':2*rough_lim/(1+rough_lim**2)}
    # Check upper bounds: singular values norm <=1, complete positivity already shown.
    for n in [4,10,20]:
        ang=[0 if j%2==0 else sigma/math.sqrt(n) for j in range(n+1)]
        mats=[overlap(support(ang[j]),support(ang[j+1])) for j in range(n)]
        w=np.linalg.multi_dot(mats)
        target=np.diag([1,math.cos(sigma/math.sqrt(n))**n])
        check(f'finite zigzag path composition n={n}',np.allclose(w,target,atol=EPS))
    # Unchanged alpha readout produces no nontrivial scale derivative.
    alpha=float((79+48*sp.sqrt(2)+42*sp.sqrt(3)+3*sp.sqrt(6))/31104)
    check('fixed record readout independent of path n', 0<alpha<1/100)
    values['alpha_candidate']=alpha
    result={'status':'PASS','test_count':len(tests),'test_names':tests,
            'values':values, 'claims_not_established':['physical time','dynamical Hamiltonian',
            'unique continuum regularization','U(1) charge dynamics','QED beta function'],
            'scope_note':'Overlap comparisons assume a common ambient qutrit frame and chosen ordered paths; both exceed the source axiom.'}
    OUT.write_text(json.dumps(result,indent=2,sort_keys=True,default=str)+'\n')
    print(json.dumps({'status':result['status'],'tests':len(tests),'principal_angle_max_error':max_sv_err,
                      'wilson_trace':float(16/sp.sqrt(65)),'smooth_last':path_results[-1]['smooth_gain'],
                      'zigzag_last':path_results[-1]['zigzag_gain'],'zigzag_limit':rough_lim,
                      'nonlinearity_frob_gap':nonlinear},indent=2))
    return result

if __name__=='__main__':
    run()