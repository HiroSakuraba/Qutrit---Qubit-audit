#!/usr/bin/env python3
"""Chain Theory v15: local POVM and sharp Bell-join refinement.

This tests a chosen v11/v14 embedding, not a uniquely axiomatic measurement.
Uses SymPy for exact results and the preexisting v14 source construction
for an independent numerical verification of all 172 sector contributions.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import sympy as s
from chain_beta_extension_v14 import construct_effects, WEIGHTS, ALPHA

ROOT=Path(__file__).resolve().parent
R=s.Rational
I=s.I
S=s.zeros(9)
for a in range(3):
    for b in range(3):
        S[3*a+b,3*b+a]=1
ONE=s.eye(9)
J=s.zeros(9,1)
for a in range(3):J[4*a]=1/s.sqrt(3)
P=s.simplify(J*J.H)
F=ONE-P
S_NUM=np.asarray(S,dtype=np.complex128)
P_NUM=np.asarray(P,dtype=np.complex128)

EXACT={
    '1':(55*ONE-7*S)/20736,
    'sqrt2':ONE/648,
    'sqrt3':(77*ONE+21*S)/62208,
    'sqrt6':ONE/10368,
}
EXPECTED_MULT={'1':79,'sqrt2':48,'sqrt3':42,'sqrt6':3}
W_EXACT={'1':s.Integer(1),'sqrt2':s.sqrt(2),'sqrt3':s.sqrt(3),'sqrt6':s.sqrt(6)}
ALPHA_EXACT=(79+48*s.sqrt(2)+42*s.sqrt(3)+3*s.sqrt(6))/31104
RHO_IN=s.zeros(9)
for a in range(2):
    for b in range(2):RHO_IN[3*a+b,3*a+b]=R(1,4)
RHO_MIX=ONE/9


def trace(m):return s.trace(m)

def numeric_localization():
    """Trace out q1,q2 (d=4), Q3,f (d=6) from v14 physical effects."""
    data=construct_effects()
    sectors=data['sector']
    found={}
    residual=0.0
    for name, read_effect in sectors.items():
        assert read_effect.shape==(36,36)
        # Ordering Q1 q1 Q2 q2 (each index on a read link).
        reduced=np.einsum('abcd ebfd -> acef',
            read_effect.reshape(3,2,3,2,3,2,3,2)).reshape(9,9)
        effect=reduced/(4*6*144)
        exact=np.asarray(EXACT[name],dtype=np.complex128)
        err=float(np.max(np.abs(effect-exact)))
        residual=max(residual,err)
        assert err<1e-14,(name,err)
        found[name]=effect
    assert np.linalg.norm(sum(found.values())-float(R(172,31104))*np.eye(9))<1e-13
    return found,residual


def exact_tests():
    assert S.T==S and S*S==ONE
    Pplus=(ONE+S)/2;Pminus=(ONE-S)/2
    assert Pplus*Pplus==Pplus and Pminus*Pminus==Pminus
    assert trace(Pplus)==6 and trace(Pminus)==3
    assert S*J==J and F*F==F
    assert s.simplify(sum(EXACT.values(),s.zeros(9))-R(172,31104)*ONE)==s.zeros(9)
    Ezero=(ONE-sum(EXACT.values(),s.zeros(9)))
    assert Ezero==R(30932,31104)*ONE
    eigs={}
    for name,E in EXACT.items():
        assert s.simplify(P*E-E*P)==s.zeros(9)
        eplus=s.simplify(trace(Pplus*E)/6)
        eminus=s.simplify(trace(Pminus*E)/3)
        assert eplus>0 and eminus>0
        assert s.simplify(E-eplus*Pplus-eminus*Pminus)==s.zeros(9)
        assert s.simplify(trace(E)-R(EXPECTED_MULT[name],3456))==0
        eigs[name]={'symmetric_multiplicity_6':str(eplus),'antisymmetric_multiplicity_3':str(eminus)}
    # The swap expectation gives the entire score on arbitrary qutrit pair state.
    t=s.symbols('t',real=True)
    f=(55-7*t)/20736+s.sqrt(2)/648+(77+21*t)*s.sqrt(3)/62208+s.sqrt(6)/10368
    vals={}
    for name,rho in [('maximally_mixed',RHO_MIX),('original_unjoined',RHO_IN),('selected_Bell',P)]:
        sv=s.simplify(trace(rho*S))
        probs={k:s.simplify(trace(rho*E)) for k,E in EXACT.items()}
        calc=s.simplify(sum(W_EXACT[k]*probs[k] for k in probs))
        assert s.simplify(calc-f.subs(t,sv))==0
        vals[name]={'swap_expectation_exact':str(sv),
                    'sector_probabilities_exact':{k:str(v) for k,v in probs.items()},
                    'weighted_score_exact':str(calc),
                    'weighted_score_numeric':float(calc.evalf(17))}
    assert s.simplify(f.subs(t,R(1,3))-ALPHA_EXACT)==0
    delta=s.simplify(f.subs(t,R(1,2))-ALPHA_EXACT)
    assert delta==7*(s.sqrt(3)-1)/124416
    # exact local Bell-join refinement: selected K0=P3; failure K_w=sqrt(E_w)F.
    # E_w commute with P3, so the outcome effects are F E_w F, summing to F.
    assert s.simplify(P+F*(Ezero+sum(EXACT.values(),s.zeros(9)))*F-ONE)==s.zeros(9)
    seq={name:{k:s.simplify(trace(rho*F*E*F)) for k,E in EXACT.items()}
         for name,rho in [('original_unjoined',RHO_IN),('selected_Bell',P),('maximally_mixed',RHO_MIX)]}
    seq_score={name:s.simplify(sum(W_EXACT[k]*p for k,p in probs.items())) for name,probs in seq.items()}
    assert seq_score['selected_Bell']==0
    # The source P3 outcome has its literal probability 1/6.
    assert trace(RHO_IN*P)==R(1,6)
    # Gauge test: W=(U tensor conjugate(U)) fixes P3 but not the swap-dependent effects.
    U=s.diag(I,-I,1)
    assert s.det(U)==1
    W=s.kronecker_product(U,U.conjugate())
    assert W.H*W==ONE
    assert s.simplify(W*P*W.H-P)==s.zeros(9)
    swapdiff=s.simplify(W*S*W.H-S)
    norm_squared=s.simplify(trace(swapdiff.H*swapdiff))
    assert norm_squared==16
    changes={k:s.simplify(trace((W*E*W.H-E).H*(W*E*W.H-E))) for k,E in EXACT.items()}
    assert changes['1']==R(784,20736**2) # 16*(7/20736)**2 = 784/20736**2
    assert changes['sqrt3']==changes['1']
    assert changes['sqrt2']==0 and changes['sqrt6']==0
    # Reassure: transported effect and state have unchanged probabilities.
    for name,rho in [('initial',RHO_IN),('joined',P)]:
        for E in EXACT.values():
            assert s.simplify(trace(rho*E)-trace(W*rho*W.H*W*E*W.H))==0
    # An entire one-parameter family of 2D support orientations returns different scores.
    # P1=diag(1,1,0); P2 projects onto span{e0, cos(theta)e1+sin(theta)e2}.
    # s=Tr[(P1/2)(P2/2)] = (1+cos(theta)^2)/4.
    span=s.Matrix.hstack(s.Matrix([1,0,0]),s.Matrix([0,1/s.sqrt(3),s.sqrt(R(2,3))]))
    P2=s.simplify(span*span.H)
    assert P2*P2==P2
    rotate_rho=s.kronecker_product(s.diag(1,1,0)/2,P2/2)
    assert trace(rotate_rho)==1
    assert trace(rotate_rho*S)==R(1,3)
    return {'local_effect_spectra':eigs,'coarse_nonzero_effect':str(R(172,31104)),
        'zero_effect':str(R(30932,31104)),
        'source_state_scores':vals,
        'delta_initial_minus_alpha_exact':str(delta),
        'delta_initial_minus_alpha_numeric':float(delta.evalf(17)),
        'sequential_join_plus_records':{
            name:{'sector_probabilities_exact':{k:str(v) for k,v in probs.items()},
                  'score_exact':str(seq_score[name]),'score_numeric':float(seq_score[name].evalf(17))}
            for name,probs in seq.items()},
        'matched_frame_witness':{
            'U':'diag(i,-i,1)',
            'P3_invariant_exact':True,
            'squared_swap_difference_HS':'16',
            'sector_effect_difference_HS_squared':{k:str(v) for k,v in changes.items()},
            'fixed_effects_invariant_under_this_matched_frame':False,
            'jointly_transported_effects_are_covariant':True},
        'misaligned_support_state_with_exact_original_alpha':{
            'cos_theta_squared':'1/3','swap_expectation':'1/3',
            'score_matches_alpha_exact':True,
            'caution':'The relative support orientation is an added state choice, not an axiom-selected state.'}}


def random_unitary(rng,d):
    z=rng.normal(size=(d,d))+1j*rng.normal(size=(d,d))
    q,r=np.linalg.qr(z)
    q=q@np.diag(np.conjugate(np.diag(r)/np.abs(np.diag(r))))
    q=q/np.linalg.det(q)**(1/d)
    return q


def sqrt_effect(E):
    vals,vecs=np.linalg.eigh(E)
    assert vals.min()>-1e-13
    return (vecs*np.sqrt(vals.clip(min=0)))@vecs.conj().T


def numeric_instrument_tests(effects):
    rng=np.random.default_rng(20261025)
    p=P_NUM;ff=np.eye(9)-p
    zero=np.eye(9)-sum(effects.values())
    kraus={'Bell':p}
    for k,E in effects.items():kraus[k]=sqrt_effect(E)@ff
    kraus['record_zero']=sqrt_effect(zero)@ff
    completeness=sum((K.conj().T@K for K in kraus.values()),np.zeros((9,9),complex))
    assert np.max(np.abs(completeness-np.eye(9)))<1e-13
    # The original pre-join two-link state in 9 x 4 reshaping (unread factor remains spectator).
    psi=np.zeros((9,4),complex)
    for a in range(2):
        for b in range(2):psi[3*a+b,2*a+b]=1/2
    rhoB=psi.conj().T@psi
    assert np.max(np.abs(rhoB-np.eye(4)/4))<1e-14
    selected=kraus['Bell']@psi
    psucc=float(np.vdot(selected,selected).real)
    assert abs(psucc-1/6)<1e-14
    bell2=np.zeros(4,complex);bell2[[0,3]]=1/np.sqrt(2)
    assert np.linalg.norm(selected-np.asarray(J,dtype=complex)@bell2[np.newaxis,:]/np.sqrt(6))<1e-13
    sum_remote=sum(((K@psi).conj().T@(K@psi) for K in kraus.values()),np.zeros((4,4),complex))
    no_signal=float(np.max(np.abs(sum_remote-rhoB)))
    assert no_signal<1e-13
    # Frame transport covariance; contrasts with holding E fixed.
    max_cov=0.0;max_comp=0.0;max_nosig=0.0
    for trial in range(32):
        u=random_unitary(rng,3)
        W=np.kron(u,u.conj())
        transformed={k:W@K@W.conj().T for k,K in kraus.items()}
        new_psi=W@psi
        for k,K in kraus.items():
            a=float(np.vdot(K@psi,K@psi).real)
            b=float(np.vdot(transformed[k]@new_psi,transformed[k]@new_psi).real)
            max_cov=max(max_cov,abs(a-b))
        c=sum((K.conj().T@K for K in transformed.values()),np.zeros((9,9),complex))
        max_comp=max(max_comp,float(np.linalg.norm(c-np.eye(9))))
        rb=sum(((K@new_psi).conj().T@(K@new_psi) for K in transformed.values()),np.zeros((4,4),complex))
        max_nosig=max(max_nosig,float(np.linalg.norm(rb-rhoB)))
    assert max_cov<1e-13 and max_comp<1e-12 and max_nosig<1e-12
    return {'physical_join_system_dim':9,'remote_read_qubits_dim':4,'unread_pair_dim':6,
            'instrument_outcomes':list(kraus),'all_kraus_act_only_on_joining_qutrits':True,
            'selected_Kraus_is_exact_P3':True,'selected_probability':psucc,
            'trace_preserving_completeness_max_error':float(np.max(np.abs(completeness-np.eye(9)))),
            'unconditional_no_signaling_max_error':no_signal,
            'random_correlated_matching_trials':32,
            'max_transport_probability_error':max_cov,
            'max_transport_completeness_error':max_comp,
            'max_transport_no_signaling_error':max_nosig}


def main():
    effects,residual=numeric_localization()
    exact=exact_tests()
    physical=numeric_instrument_tests(effects)
    report={
        'program':'Chain Theory QED beta, v15 — local instrument and gauge-descent audit',
        'status':'Source-motivated local CPTP instrument constructed; exact original alpha not recovered on source pure state; instrument not uniquely selected.',
        'readout_model':'v11 172 record embedding and v14 Stinespring operator-frame bridge, then normalized trace over remote/read spectator factors',
        'record_counts':EXPECTED_MULT,
        'alpha_candidate':float(ALPHA_EXACT.evalf(18)),
        'max_residual_localization_vs_exact_formula':residual,
        'symbolic_audit':exact,
        'physical_instrument':physical,
        'beta_status':'No momentum or scale dependence derived, hence no independently determined QED beta function.'}
    dest=ROOT/'chain_beta_extension_v15_results.json'
    dest.write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))
if __name__=='__main__':main()