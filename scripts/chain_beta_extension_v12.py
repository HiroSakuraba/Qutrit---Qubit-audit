#!/usr/bin/env python3
"""Chain Theory v12 — instrument-valued gauge descent, conditional and explicitly scoped.

Constructs a mathematically valid Lüders instrument on the ABSTRACT readout
operator Hilbert space, not automatically a CP instrument on the underlying
216D physical Hilbert space. Tests gauge-covariance, nonuniqueness and a
precise operational obstruction. Requires NumPy, SymPy and v9/v11 scripts.
"""
import json
import math
from pathlib import Path
import numpy as np
import sympy as sp
from chain_beta_extension_v9 import make_basis,record_matrix,transport
from chain_beta_extension_v11 import build_canonical_records,qutrit_match,bell3

N=31104
sqrt=np.sqrt
S=79+48*sp.sqrt(2)+42*sp.sqrt(3)+3*sp.sqrt(6)
ALPHA=S/N
ROOT=Path(__file__).resolve().parent

def commutator_twirl_exact():
    I=sp.eye(3)
    Z=sp.diag(0,0,1)
    P=sp.Rational(1,2)*sp.BlockMatrix([[I,-I],[-I,I]]).as_explicit()
    total=sp.zeros(6)
    mats=[]
    # theta_j=j*pi/8, qubit rotations u(theta)=diag(exp(i theta),exp(-i theta)).
    # Their adjoint SO(3) rotations are R(+2 theta), R(-2 theta).
    for j in range(4):
        theta=sp.pi*sp.Rational(j,8)
        a=2*theta
        R1=sp.Matrix([[sp.cos(a),-sp.sin(a),0],[sp.sin(a),sp.cos(a),0],[0,0,1]])
        R2=R1.T
        W=sp.diag(R1,R2)
        total+=W*P*W.T
        mats.append(W)
    Tw=sp.simplify(total/4)
    expected=sp.Rational(1,2)*sp.BlockMatrix([[I,-Z],[-Z,I]]).as_explicit()
    assert sp.simplify(Tw-expected)==sp.zeros(6)
    assert sp.simplify(sp.trace(Tw))==3
    assert Tw.rank()==5
    eig={str(k):int(v) for k,v in Tw.eigenvals().items()}
    assert eig=={'1/2':4,'1':1,'0':1},eig
    assert all(sp.simplify(W*Tw*W.T-Tw)==sp.zeros(6) for W in mats)
    assert sp.simplify(sp.trace((Tw-P)**2))==1 # 4 entries (1/4)^2 etc
    return {'projector_rank_before':3,'twirled_effect_rank':5,
            'eigenvalues_and_multiplicities':eig,'trace_before_and_after':'3',
            'HS_squared_change_exact':str(sp.trace((Tw-P)**2))}

def physical_CP_counterexample():
    A=sp.diag(1,-1)/sp.sqrt(2)
    X=sp.diag(1,0)
    # |A)(A| is an orthogonal projector on Liouville / HS space.
    out=sp.simplify(A*sp.trace(A*X))
    assert out==sp.diag(sp.Rational(1,2),-sp.Rational(1,2))
    assert sp.trace(out)==0 and out[1,1]<0
    # This map is not even positive; it cannot be a quantum instrument on physical states.
    return {'A':'sigma_z/sqrt(2)','input':'|0><0|',
            'output':'diag(1/2,-1/2)','output_positive':False}

def main():
    assert sp.simplify(S-(79+48*sp.sqrt(2)+42*sp.sqrt(3)+3*sp.sqrt(6)))==0
    assert 226 < S.evalf() < 227 and 172 < S.evalf()
    basis,blocks,singles,_=make_basis()
    M,P=record_matrix(blocks,singles)
    E,w,tags=build_canonical_records(blocks,singles)
    assert np.linalg.norm(E.T@E-np.eye(172))<3e-13
    group={'1':[],'sqrt2':[],'sqrt3':[],'sqrt6':[]}
    for j,x in enumerate(w):
        group[min(group,key=lambda k:abs({'1':1,'sqrt2':sqrt(2),'sqrt3':sqrt(3),'sqrt6':sqrt(6)}[k]-x))].append(j)
    assert {k:len(v) for k,v in group.items()}=={'1':79,'sqrt2':48,'sqrt3':42,'sqrt6':3}
    proj={k:E[:,ids]@E[:,ids].T for k,ids in group.items()}
    proj['0']=np.eye(175)-P
    ident=np.eye(175)
    for k,A in proj.items():
        assert np.linalg.norm(A@A-A)<2e-13,k
        for l,B in proj.items():
            if k!=l: assert np.linalg.norm(A@B)<2e-13,(k,l)
    assert np.linalg.norm(sum(proj.values())-ident)<2e-13
    assert np.linalg.norm(sum({'1':1,'sqrt2':sqrt(2),'sqrt3':sqrt(3),'sqrt6':sqrt(6)}[k]*proj[k]
                              for k in group)-M)<2e-13
    # E_w are projections in the FULL 31104 operator-HS Hilbert space. The
    # 175D computation is valid because these projectors factor through one
    # unread connected vector; the rest of the full space is zero-weight.
    full_ranks={'0':N-172,'1':79,'sqrt2':48,'sqrt3':42,'sqrt6':3}
    probs={k:float(v/N) for k,v in full_ranks.items()}
    assert abs(sum(probs.values())-1)<1e-15
    observable_mean=sum(float({'1':sp.Integer(1),'sqrt2':sp.sqrt(2),
                               'sqrt3':sp.sqrt(3),'sqrt6':sp.sqrt(6)}.get(k,0))*probs[k]
                        for k in full_ranks)
    assert abs(observable_mean-float(ALPHA))<1e-16
    # Concrete covariant instrument under a matched-frame gauge
    # transformation that stabilizes the normalized three-pair state.
    u=np.diag(np.exp(1j*np.array([np.pi/4,-np.pi/4])))
    uc=u.conj()
    A=qutrit_match(u); B=qutrit_match(uc)
    Ujoin=np.kron(A,B)
    assert np.linalg.norm(Ujoin@bell3()-bell3())<3e-15
    W=transport(blocks,singles,uq1=u,uQ1=A,uq2=uc,uQ2=B)
    assert np.linalg.norm(W@W.T-ident)<3e-13
    ket_comm=E[:,169]
    assert tags[169]=='comm'
    psi=ket_comm.astype(complex)
    moved=W@psi
    pre_probs={k:float(np.vdot(psi,mat@psi).real) for k,mat in proj.items()}
    fixed_frame_probs={k:float(np.vdot(moved,mat@moved).real) for k,mat in proj.items()}
    transported_probs={k:float(np.vdot(moved,(W@mat@W.T)@moved).real) for k,mat in proj.items()}
    err_cov=max(abs(pre_probs[k]-transported_probs[k]) for k in proj)
    assert err_cov<1e-13
    assert abs(pre_probs['sqrt6']-1)<1e-13
    assert abs(fixed_frame_probs['sqrt6'])<1e-13
    assert abs(fixed_frame_probs['0']-1)<1e-13
    # Test branch STATE covariance on generic coherent superposition.
    missing=(np.eye(175)[:,singles['q1'][0]]+np.eye(175)[:,singles['q2'][0]])/sqrt(2)
    random_state=(2*E[:,0].astype(complex)+1j*E[:,90]+1.3*E[:,171]+0.7j*missing)
    random_state/=np.linalg.norm(random_state)
    random_moved=W@random_state
    generic_probs=[]
    max_branch_err=0.0
    for k,mat in proj.items():
        source=mat@random_state
        target=(W@mat@W.T)@random_moved
        max_branch_err=max(max_branch_err,float(np.linalg.norm(target-W@source)))
        generic_probs.append(float(np.vdot(random_state,mat@random_state).real))
    assert max_branch_err<2e-13
    assert abs(sum(generic_probs)-1)<2e-13
    # The full physical space is dimension 216, its operator space dim 46656;
    # R of dimension 31104 is a subspace of it, not automatically a Hilbert space
    # on which physically prepared 31104-dimensional quantum states reside.
    assert 6**3==216 and 216**2==46656
    # Probability obstruction, and an infinite family of binary POVMs all
    # giving the same candidate probability on tau=I/N.
    assert int(sp.ceiling(S))==227
    must_add=int(sp.ceiling(S))-172
    binary_effects=[]
    for t in [0.0,0.25,0.5,float(ALPHA),0.75,1.0]:
        a=float((S-172*sp.Float(t))/(N-172))
        assert 0<=t<=1 and 0<=a<=1
        p=(172*t+(N-172)*a)/N
        assert abs(p-float(ALPHA))<2e-16
        binary_effects.append({'t_records':t,'a_complement':a,'probability_at_tau':p})
    assert abs(binary_effects[3]['a_complement']-float(ALPHA))<2e-14
    data={
      'research_scope':'v12 mathematical instrument on operator-HS readout space; not yet a physical 216D CP instrument',
      'source_dimensions':{'physical_three_arm_state_space':216,'physical_operator_space':46656,
                           'specified_readout_operator_space':N,'record_rank':172},
      'spectral_instrument':{'outcome_multiplicities':full_ranks,
                             'outcome_probabilities_at_tau':probs,
                             'weighted_mean_at_tau':observable_mean,
                             'alpha_exact':str(ALPHA),
                             'alpha_numeric':float(ALPHA),
                             'kraus_completeness':'sum_w P_w=I_31104',
                             'note':'These CP-TP Kraus operators act on an abstract HS Hilbert space, not automatically on physical 216D density matrices.'},
      'matched_frame_covariance':{'three_pair_join_stabilized':True,
                    'max_outcome_probability_covariance_error':err_cov,
                    'max_postmeasurement_branch_covariance_error':max_branch_err,
                    'commutator_selected_outcome_probability_before':pre_probs['sqrt6'],
                    'after_frame_transform_without_transport':fixed_frame_probs['sqrt6'],
                    'after_covariant_instrument_transport':transported_probs['sqrt6'],
                    'fixed_frame_missing_sector_probability':fixed_frame_probs['0']},
      'no_go_probability':{'tr_M':float(S),'max_trace_yes_effect_supported_on_172':172,
                           'minimum_rank_of_binary_effect_matching_alpha_at_tau':227,
                           'minimum_new_support_directions_if_original_172_all_retained':must_add,
                           'nonunique_binary_POVMs':binary_effects},
      'commutator_sector_twirling':commutator_twirl_exact(),
      'physical_CP_failure_of_naive_superoperator':physical_CP_counterexample(),
      'open':['No canonical map from the physical 216D state to a readout-space density operator was specified.',
              'The original model does not define a physical CP instrument for the abstract 172 weighted operators.',
              'Gauge covariance of transported measurement operators does not select them uniquely.',
              'No QED momentum dependence, beta coefficient, or spacetime dynamics is derived.']
    }
    output=ROOT/'chain_beta_extension_v12_results.json'
    output.write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(data,indent=2))
    print('PASS v12 — mathematical spectral instrument, matched covariance, POVM bound, exact twirl, CP counterexample')

if __name__=='__main__': main()