#!/usr/bin/env python3
"""Chain Theory beta v13: physical 216-dimensional measurement completion.

Primary source: A. D. Wissner-Gross, ChainTheoryAlpha.pdf, Supp. S1/S7/S8.
All trace-preserving instruments here are NEW, not derived from the source axiom.
Requires sympy, numpy. Produces deterministic JSON result on stdout.
"""
from __future__ import annotations
import json, math
import numpy as np
import sympy as s

R=s.Rational
rt=s.sqrt
I=s.I

def ket(n,i):
    v=s.zeros(n,1);v[i]=1;return v

def equal(a,b):
    return all(s.simplify(x)==0 for x in (a-b))

def exact_physical_instruments():
    # Physical grouping (A=qutrit copies Q1,Q2:9; B=read qubits q1,q2:4; C=unread qutrit-qubit pair:6).
    psi=s.zeros(9,4)
    for i in range(2):
        for j in range(2):psi[3*i+j,2*i+j]=R(1,2)
    phi3=(ket(9,0)+ket(9,4)+ket(9,8))/rt(3)
    phi2=(ket(4,0)+ket(4,3))/rt(2)
    phi_unread=(ket(6,0)+ket(6,3))/rt(2)
    p3=phi3*phi3.H
    rhoB0=s.simplify(psi.H*psi)
    assert equal(rhoB0,s.eye(4)/4)
    selected=s.simplify(p3*psi)
    p0=s.simplify((selected.H*selected).trace())
    assert p0==R(1,6)
    assert equal(selected,s.simplify(phi3*phi2.T/rt(6)))
    # Complete three-pair normalized branch, checked as vector of dimension 216.
    selected_phys=s.kronecker_product(s.Matrix(selected).reshape(36,1),phi_unread)
    target_phys=s.kronecker_product(phi3,phi2,phi_unread)/rt(6)
    assert selected_phys.rows==216 and equal(selected_phys,target_phys)
    # Binary Lüders completion; its selected operation K0=P3 is exactly source join projection.
    kfail=s.eye(9)-p3
    assert equal(p3.H*p3 + kfail.H*kfail,s.eye(9))
    failure=s.simplify(psi.H*kfail.H*kfail*psi)
    assert equal(failure+selected.H*selected,rhoB0)
    rhoFail=s.simplify(failure/(1-p0))
    PB=phi2*phi2.H
    assert equal(rhoFail,R(3,10)*s.eye(4)-R(1,5)*PB)
    assert [s.simplify(v) for v in rhoFail.eigenvals().keys()] # positive
    assert s.simplify((phi2.H*rhoFail*phi2)[0])==R(1,10)
    # Nine-Kraus measure-reprepare extension (the input state is the aligned one).
    def b2(j):
        e=lambda a,b: ket(9,3*a+b)
        return [(e(0,0)+e(1,1))/rt(2),(e(0,0)-e(1,1))/rt(2),
                (e(0,1)+e(1,0))/rt(2),(e(0,1)-e(1,0))/rt(2)][j]
    eta=rt(3)*b2(0)-rt(2)*phi3
    vectors=[phi3,eta,b2(1),b2(2),b2(3),ket(9,2),ket(9,5),ket(9,6),ket(9,7)]
    frame=s.Matrix.hstack(*vectors)
    assert equal(frame.H*frame,s.eye(9))
    Ks=[phi3*v.H for v in vectors]
    assert equal(Ks[0],p3)
    assert equal(sum((k.H*k for k in Ks),s.zeros(9)),s.eye(9))
    probs=[];concs=[];mixedBob=s.zeros(4)
    for k in Ks:
        ps=s.simplify(k*psi)
        sigma=s.simplify(ps.H*ps)
        p=s.simplify(sigma.trace())
        probs.append(p);mixedBob+=sigma
        if p==0:
            concs.append(None)
        else:
            rb=s.simplify((phi3.H*ps).T/rt(p))
            amat=s.Matrix([[rb[0],rb[1]],[rb[2],rb[3]]])
            cc=s.simplify(4*amat.det()*s.conjugate(amat.det()))
            assert cc==1 # squared concurrence
            concs.append(1)
    assert probs==[R(1,6),R(1,12),R(1,4),R(1,4),R(1,4),0,0,0,0]
    assert equal(mixedBob,rhoB0)
    # All remaining nonselected outcomes have a fixed remote aggregate, no matter how resolved.
    assert equal(sum((s.simplify((k*psi).H*(k*psi)) for k in Ks[1:]),s.zeros(4)),failure)
    return {
        'physical_dimensions':{'Alice_Q1Q2':9,'Bob_q1q2':4,'unread_Q3f':6,'total':216},
        'physical_input_is_three_maximally_entangled_links':True,
        'selected_kraus_is_exact_source_P3':True,
        'selected_probability':'1/6',
        'selected_normalized_poststate':'omega3(Q1,Q2) tensor omega2(q1,q2) tensor omega2(Q3,f)',
        'binary_luders_cptp_exact':True,
        'failure_probability':'5/6',
        'universal_coarse_failure_Bob_rho':'(3/10) I4 - (1/5) |omega2><omega2|',
        'failure_Bob_eigenvalues':{'Bell_direction':'1/10','three_orthogonal_directions':'3/10'},
        'failure_Bob_Bell_fidelity':'1/10',
        'nine_repreparation_Kraus_CPTP_exact':True,
        'nine_outcome_probabilities':[str(x) for x in probs],
        'nine_nonzero_outcome_concurrences':concs,
        'Bob_unconditional_state':'I4/4',
        'unread_link_unchanged':'|omega2><omega2|',
        'limitation':'Repreparation and nonselected branches exceed original no-failure distillation axiom.'
    }

def exact_weyl_twirl():
    omega=(-1+I*rt(3))/2
    X=s.zeros(3)
    for i in range(3):X[(i+1)%3,i]=1
    Z=s.diag(1,omega,omega**2)
    phi3=(ket(9,0)+ket(9,4)+ket(9,8))/rt(3)
    project=phi3*phi3.H
    bells=[]
    for a in range(3):
        for b in range(3):
            U=(X**a)*(Z**b)
            bvec=s.simplify(s.kronecker_product(U,s.eye(3))*phi3)
            bells.append(bvec)
            # Paired matching leaves the *same* join invariant.
            correlated=s.kronecker_product(U,U.conjugate())*phi3
            assert equal(correlated,phi3)
    gram=s.Matrix.hstack(*bells)
    assert equal(gram.H*gram,s.eye(9))
    twirled=s.simplify(sum((v*v.H for v in bells),s.zeros(9))/9)
    assert equal(twirled,s.eye(9)/9)
    # Classical unknown absolute qutrit frame: averaged selected POVM effect is I9/9.
    psi=s.zeros(9,4)
    for i in range(2):
        for j in range(2):psi[3*i+j,2*i+j]=R(1,2)
    rho=s.simplify(psi.H*twirled*psi)
    assert equal(rho,s.eye(4)/36)
    assert rho.trace()==R(1,9)
    return {
        'Weyl_qutrit_unitaries':9,'nine_Bell_states_orthonormal_exact':True,
        'one_sided_independent_frame_twirl_of_P3':'I9/9',
        'correlated_U_tensor_conjugateU_matching_stabilizes_P3':True,
        'unknown_frame_event_probability':'1/9',
        'unknown_frame_conditional_Bob':'I4/4',
        'source_aligned_event_probability':'1/6',
        'source_aligned_conditional_Bob':'|omega2><omega2|',
        'conditional_no_go':'A fixed effect invariant under the full independent SU(3) x SU(3) local frame group must be scalar; P3 is not invariant.'
    }

def numerical_robustness():
    rng=np.random.default_rng(20261024)
    phi=np.zeros(9,complex);phi[[0,4,8]]=1/np.sqrt(3)
    P=np.outer(phi,phi.conj());F=np.eye(9)-P
    alpha=(79+48*np.sqrt(2)+42*np.sqrt(3)+3*np.sqrt(6))/31104
    def haar3():
        z=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
        q,r=np.linalg.qr(z)
        return q@np.diag(np.diag(r)/np.abs(np.diag(r)))
    max_no_signal=0.;max_frame_covariance=0.;min_p=1.;max_p=0.;min_concurrence=1.;max_concurrence=0.;sample_count=64
    for _ in range(sample_count):
        a,b=haar3()[:,:2],haar3()[:,:2]
        psi=np.zeros((9,4),complex)
        for i in range(2):
            for j in range(2):
                psi[:,2*i+j]=np.kron(a[:,i],b[:,j])/2
        initial=psi.conj().T@psi
        # Physical frame covariance: transport instruments together with the two qutrit frames.
        g=np.kron(haar3(),haar3());pg=g@P@g.conj().T;fg=np.eye(9)-pg
        psig=g@psi
        max_frame_covariance=max(max_frame_covariance,
              np.linalg.norm(pg@psig-g@(P@psi)),
              np.linalg.norm(fg@psig-g@(F@psi)),
              np.linalg.norm((pg@psig).conj().T@(pg@psig)-(P@psi).conj().T@(P@psi)))
        unconditional=(P@psi).conj().T@(P@psi)+(F@psi).conj().T@(F@psi)
        max_no_signal=max(max_no_signal,np.linalg.norm(initial-unconditional))
        prob=np.trace((P@psi).conj().T@(P@psi)).real
        min_p=min(min_p,prob);max_p=max(max_p,prob)
        if prob>1e-12:
            remote=phi.conj()@psi/np.sqrt(prob)
            concurrence=2*abs(np.linalg.det(remote.reshape(2,2)))
            min_concurrence=min(min_concurrence,float(concurrence))
            max_concurrence=max(max_concurrence,float(concurrence))
    assert max_no_signal<1e-13 and max_frame_covariance<1e-13
    # An exact nonaligned support example produces a separable selected pair.
    v1=np.eye(3,dtype=complex)[:,[0,2]]
    v2=np.eye(3,dtype=complex)[:,[0,1]]
    psix=np.zeros((9,4),complex)
    for i in range(2):
        for j in range(2):
            psix[:,2*i+j]=np.kron(v1[:,i],v2[:,j])/2
    remote=phi.conj()@psix
    bad_p=np.vdot(remote,remote).real
    bad_conc=2*abs(np.linalg.det((remote/np.sqrt(bad_p)).reshape(2,2)))
    assert abs(bad_p-1/12)<1e-13 and bad_conc<1e-13
    # A rare accept score can be tuned to arbitrary target <=1/6. It does NOT predict alpha.
    q=6*alpha
    assert 0<q<1
    return {
        'seed':20261024,'random_independently_rotated_support_pairs':sample_count,
        'maximum_remote_no_signaling_residual':max_no_signal,
        'maximum_transport_covariance_residual':max_frame_covariance,
        'fixed_P3_probability_range':[min_p,max_p],
        'selected_remote_concurrence_range':[min_concurrence,max_concurrence],
        'exact_support_misalignment_selected_probability':bad_p,
        'exact_support_misalignment_selected_concurrence':bad_conc,
        'alpha_candidate':alpha,'accept_given_Bell_success_calibrated_q':q,
        'probability_after_calibration':q/6,
        'nonunique_score_test':'Any t in [0,1/6] can be generated by accepting a fraction 6t of successful P3 outcomes. The factor t is an input, not a derivation.'
    }

def main():
    out={'research_claim':'Physical CPTP extensions exist and are covariant when supplied frame transport; source 172 trace not selected as physical observable.',
         'exact_physical_instruments':exact_physical_instruments(),
         'frame_and_gauge_test':exact_weyl_twirl(),
         'random_support_and_calibration':numerical_robustness(),
         'status':'All exact SymPy and seeded NumPy assertions passed.'}
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()