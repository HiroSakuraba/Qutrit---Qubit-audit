#!/usr/bin/env python3
"""Chain Theory / QED beta program v14.

A canonical, unital CP Heisenberg bridge from the source's 31,104D
operator-HS space into physical 216D measurement effects. We use the
particular v11 embedding and distinguish abstract CP realizability
from locality and uniqueness. No QED beta-function coefficient is derived.

Requires NumPy, SymPy and the companion v9/v11 scripts.
"""
from __future__ import annotations
import json
from collections import Counter
from fractions import Fraction
from pathlib import Path
import numpy as np
import sympy as sp
from chain_beta_extension_v9 import make_basis, gellmann, pauli
from chain_beta_extension_v11 import build_canonical_records, exact_connected

ROOT=Path(__file__).resolve().parent
SQRT=np.sqrt
N=24*36*36
D=216
FRAME_CONSTANT=N//D
WEIGHTS={'1':1.,'sqrt2':SQRT(2.),'sqrt3':SQRT(3.),'sqrt6':SQRT(6.)}
EXPECTED={'1':79,'sqrt2':48,'sqrt3':42,'sqrt6':3}
ALPHA=(79+48*SQRT(2.)+42*SQRT(3.)+3*SQRT(6.))/N


def verify_frame_constant():
    """Hilbert--Schmidt orthonormal bases and their quadratic Casimirs."""
    ls3=gellmann(); ls2=pauli()
    assert len(ls3)==8 and len(ls2)==3
    c3=sum((a.conj().T@a for a in ls3),np.zeros((3,3),complex))
    c2=sum((a.conj().T@a for a in ls2),np.zeros((2,2),complex))
    assert np.linalg.norm(c3-(8/3)*np.eye(3))<1e-14
    assert np.linalg.norm(c2-(3/2)*np.eye(2))<1e-14
    # Each End(C6) factor has sum B_a^dagger B_a = 6 I6.
    # The unread su(3) tensor su(2) contributes (8/3)(3/2) I6=4I6.
    assert FRAME_CONSTANT==6*6*4==144
    # Thus V=(1/sqrt(144)) sum |a> tensor A_a is an isometry:
    # V^dagger V=I216. Phi(X)=V^dagger(X tensor I216)V is unital CP.
    return {'full_readout_dimension':N,'physical_dimension':D,
            'frame_constant_exact':FRAME_CONSTANT,
            'qutrit_traceless_casimir':'8/3',
            'qubit_traceless_casimir':'3/2',
            'unread_connected_casimir':'4',
            'complete_positivity':'Proven by isometric Stinespring compression',
            'unitality':'Phi(I_R)=I_phys (144I/144)'}


def construct_effects():
    """No 31104x31104 dense arrays; exploit the 172-column low-rank embedding."""
    arms,blocks,singles,_=make_basis()
    vec,vals,tags=build_canonical_records(blocks,singles)
    basis=np.stack(arms,axis=0)  # 175 orthonormal operators in End(C36)
    read_ops=np.einsum('jr,jab->rab',vec,basis)
    gram=np.einsum('rab,sab->rs',read_ops.conj(),read_ops)
    assert np.linalg.norm(gram-np.eye(172))<2e-12
    assert Counter({k:int(sum(abs(vals-v)<1e-12)) for k,v in WEIGHTS.items()})==Counter(EXPECTED)
    sector={k:np.zeros((36,36),complex) for k in WEIGHTS}
    for j,w in enumerate(vals):
        key=min(WEIGHTS,key=lambda k:abs(w-WEIGHTS[k]))
        sector[key]+=read_ops[j].conj().T@read_ops[j]
    _,cs=exact_connected()
    c=sp.simplify(2*cs/sp.sqrt(3))  # HS-normalized unread connected operator
    assert sp.trace(c*c)==1
    c2=sp.simplify(c*c)
    C=np.array(c,dtype=complex)
    CtC=C.conj().T@C
    assert np.linalg.norm(CtC-np.array(c2,dtype=complex))<2e-15
    score_read=sum(WEIGHTS[k]*S for k,S in sector.items())
    assert abs(np.trace(score_read).real-sum(EXPECTED[k]*WEIGHTS[k] for k in WEIGHTS))<2e-12
    # Phi(P_w)=sum_{record in w} A_record^dagger A_record / 144.
    # Each A_record=R_record tensor C; sector effects factor exactly.
    # All effects are positive; complement is positive since they are
    # contributions to the resolution of the identity from an orthonormal basis.
    ev={k:np.linalg.eigvalsh(S) for k,S in sector.items()}
    max_nonzero=max(np.linalg.eigvalsh(sum(sector.values())))*max(np.linalg.eigvalsh(CtC))/FRAME_CONSTANT
    assert max_nonzero<1-1e-10
    assert min(float(x.min()) for x in ev.values())>-1e-13
    assert np.linalg.norm(score_read-score_read.conj().T)<2e-13
    # Sector effect traces: tr(Phi(P_w))=rank(P_w)/144.
    for k,S in sector.items():
        assert abs(np.trace(S).real-EXPECTED[k])<1e-11
    assert abs(np.trace(CtC)-1)<1e-14
    return {'sector':sector,'read_ops':read_ops,'connected':C,'connected2':CtC,
            'connected2_exact':c2,'score_read':score_read,
            'max_eigenvalue_nonzero_effect':float(max_nonzero),
            'record_counts':EXPECTED.copy(),
            'effect_trace_checks':{k:float(np.trace(S).real/144) for k,S in sector.items()}}


def get_pair():
    b=np.zeros(6,complex)
    b[0]=b[3]=1/SQRT(2.)
    return b


def physical_states():
    b=get_pair()
    pre=np.kron(np.kron(b,b),b)
    joined_read=np.zeros((3,2,3,2),complex)
    for q in range(3):
        for z in range(2):joined_read[q,z,q,z]=1/SQRT(6.)
    post=np.kron(joined_read.reshape(-1),b)
    assert abs(np.linalg.norm(pre)-1)<1e-13 and abs(np.linalg.norm(post)-1)<1e-13
    # Verify the paper's explicit qutrit Bell-selected 216D state and p=1/6.
    psi=pre.reshape(3,2,3,2,3,2)
    j=np.eye(3,dtype=complex)/SQRT(3)
    contracted=np.einsum('ac,ajclmn->jlmn',j.conj(),psi)
    projected=np.einsum('ik,jlmn->ijklmn',j,contracted).reshape(-1)
    prob=float(np.vdot(projected,projected).real)
    assert abs(prob-1/6)<1e-14
    assert np.linalg.norm(projected-post/SQRT(6))<2e-14
    return pre,post


def expectation(Ar,CtC,vec):
    """State-vector physical probabilities under the CP-mapped five-outcome POVM."""
    resh=vec.reshape(36,6)
    p={k:float(np.einsum('ai,ab,ij,bj->',resh.conj(),S,CtC,resh).real/144)
       for k,S in Ar.items()}
    assert all(v>=-1e-14 for v in p.values())
    p['0']=float(1-sum(p.values()))
    assert 0<p['0']<=1
    score=float(sum(WEIGHTS[k]*p[k] for k in WEIGHTS))
    return p,score


def mixture_nonuniqueness(probs,physical_score):
    """Phi_t=t Phi+(1-t) Tr(X)/N *I are all unital CP and equivariant."""
    mixed={k:EXPECTED[k]/N for k in WEIGHTS}
    mixed['0']=1-sum(mixed.values())
    assert abs(sum(WEIGHTS[k]*mixed[k] for k in WEIGHTS)-ALPHA)<1e-15
    return {'maximally_mixed_effect_probabilities':mixed,
            'alpha_at_t_0':ALPHA,'alpha_at_t_half':(ALPHA+physical_score)/2,
            'alpha_at_t_1':physical_score,
            'family':'Phi_t(X)=t Phi(X)+(1-t) tr_R(X)/31104 I_phys, 0<=t<=1'}


def random_unitary(rng,d):
    z=rng.normal(size=(d,d))+1j*rng.normal(size=(d,d))
    q,r=np.linalg.qr(z)
    diag=np.diag(r)
    return q@np.diag(np.conjugate(diag/np.abs(diag)))


def gauge_covariance_checks(sector,CtC,states):
    rng=np.random.default_rng(20261025)
    delta=0.; matrix_delta=0.
    for _ in range(12):
        u=[np.kron(random_unitary(rng,3),random_unitary(rng,2)) for _ in range(3)]
        u_read=np.kron(u[0],u[1]); u_unread=u[2]
        Bnew={k:u_read@S@u_read.conj().T for k,S in sector.items()}
        Cnew=u_unread@CtC@u_unread.conj().T
        U=np.kron(u_read,u_unread)
        matrix_delta=max(matrix_delta,max(np.linalg.norm(Bnew[k]-u_read@sector[k]@u_read.conj().T) for k in sector))
        for st in states:
            original,_=expectation(sector,CtC,st)
            rotated,_=expectation(Bnew,Cnew,U@st)
            delta=max(delta,max(abs(original[k]-rotated[k]) for k in original))
    assert delta<5e-13
    return {'trials':12,'max_probability_covariance_residual':delta,
            'transported_effects_transform_by_conjugation':True}


def nonlocality_and_exact_checks(CtC):
    _,C0=exact_connected()
    Cn=sp.simplify(2*C0/sp.sqrt(3))
    x=sp.kronecker_product(sp.eye(3),sp.Matrix([[0,1],[1,0]]))
    comm=sp.simplify((Cn*Cn)*x-x*(Cn*Cn))
    nsq=sp.simplify(sp.trace(comm.H*comm))
    assert nsq==sp.Rational(8,9)
    # Alice's original selected Bell effect acts on Q1,Q2 alone and
    # MUST commute with the Pauli X of Bob's unread qubit f.
    # Phi(P_w) contains CtC on unread Q3f and generally does not.
    assert np.linalg.norm(CtC@np.array(x,dtype=complex)-np.array(x,dtype=complex)@CtC)>0.9
    return {'exact_unread_remote_qubit_X_commutator_HS_squared':str(nsq),
            'local_Alice_only_effect_possible_for_this_POVM':False,
            'reason':'The positive nonzero-record effects act nontrivially on the unread remote qubit f.',
            'normalized_unread_connected_eigenvalues':{'sqrt(3)/2':1,'-sqrt(3)/6':3,'0':2}}


def main():
    frame=verify_frame_constant()
    effects=construct_effects()
    pre,post=physical_states()
    p_pre,a_pre=expectation(effects['sector'],effects['connected2'],pre)
    p_post,a_post=expectation(effects['sector'],effects['connected2'],post)
    assert abs(a_pre-0.03302342926119737)<2e-14
    assert abs(a_post-0.0342898130578542)<2e-14
    assert abs(a_pre-ALPHA)>0.02 and abs(a_post-ALPHA)>0.02
    # Exact sector-probability fractions from the explicit record basis.
    expected_pre={'1':sp.Rational(103,9216),'sqrt2':sp.Rational(1,144),
        'sqrt3':sp.Rational(175,27648),'sqrt6':sp.Rational(1,2304)}
    expected_post={'1':sp.Rational(55,6912),'sqrt2':sp.Rational(1,108),
        'sqrt3':sp.Rational(25,3456),'sqrt6':sp.Rational(1,3456)}
    for k in WEIGHTS:
        assert abs(p_pre[k]-float(expected_pre[k]))<2e-15
        assert abs(p_post[k]-float(expected_post[k]))<2e-15
    orig=mix=mixture_nonuniqueness(p_pre,a_pre)
    cov=gauge_covariance_checks(effects['sector'],effects['connected2'],[pre,post])
    exact=nonlocality_and_exact_checks(effects['connected2'])
    out={
      'research_question':'Does the source record readout canonically become a physical measurement on the 216D three-link system?',
      'status':'Constructive unital CP bridge exists conditional on v11 embedding; it is not uniquely selected and does not recover alpha for the source pure states.',
      'source':'ChainTheoryAlpha.pdf supplemental S6-S8; v11 record embedding',
      'frame_calculation':frame,
      'alpha_candidate':ALPHA,
      'record_weight_multiplicities':EXPECTED,
      'physical_state_expectations':{
          'initial_unjoined_three_pairs':{'sector_probabilities':p_pre,'score':a_pre,
                'exact_nonzero_probabilities':{k:str(v) for k,v in expected_pre.items()}},
          'postselected_joined_three_pairs':{'sector_probabilities':p_post,'score':a_post,
                'exact_nonzero_probabilities':{k:str(v) for k,v in expected_post.items()}},
          'maximally_mixed_physical_216D':{'score':ALPHA}},
      'maximum_eigenvalue_of_nonzero_effect':effects['max_eigenvalue_nonzero_effect'],
      'gauge_covariance':cov,
      'covariant_CP_nonuniqueness':mix,
      'locality':exact,
      'beta_function_status':'Still not defined by original source; no spacetime, propagator or scale-dependent flow derived.'
    }
    print(json.dumps(out,indent=2))
    with (ROOT/'chain_beta_extension_v14_results.json').open('w') as f:json.dump(out,f,indent=2)

if __name__=='__main__':main()