#!/usr/bin/env python3
"""Chain Theory beta v6: explicit CPTP completions and 172-record audit.

Primary source: A. D. Wissner-Gross (2026), ChainTheoryAlpha.pdf,
Supplement §§S1, S7, S8. All CPTP instruments defined here are NEW
operational hypotheses, not consequences of its 'no-failure' axiom.

Python >= 3.10, sympy, numpy. Run: python chain_beta_extension_v6.py
"""
from __future__ import annotations
import json, math
import numpy as np
import sympy as s

SQ=s.sqrt
R=s.Rational
I=s.I

def basis(n,k):
    v=s.zeros(n,1);v[k,0]=1;return v

def embed_bell(j):
    # qutrit A1,A2, but each has rank-2 support selected by the original links
    e=lambda a,b:basis(9,3*a+b)
    return ((e(0,0)+e(1,1)), (e(0,0)-e(1,1)),
            (e(0,1)+e(1,0)), (e(0,1)-e(1,0)))[j]/SQ(2)

def exact_instruments():
    # Rows: qutrit A1,A2; columns: remote B1,B2 (qubits).
    psi=s.zeros(9,4)
    for i in range(2):
        for j in range(2): psi[3*i+j,2*i+j]=R(1,2)
    assert (psi.T*psi)==s.eye(4)/4
    omega=(basis(9,0)+basis(9,4)+basis(9,8))/SQ(3)
    b=[embed_bell(k) for k in range(4)]
    eta=SQ(3)*b[0]-SQ(2)*omega
    c=[basis(9,2),basis(9,5),basis(9,6),basis(9,7)]
    # Orthogonal measurement in 9d; each outcome prepares the same maximally
    # entangled qutrit pair. First Kraus is the SOURCE'S literal P3 projector.
    v=[omega,eta,b[1],b[2],b[3]]+c
    basis_mat=s.Matrix.hstack(*v)
    assert s.simplify(basis_mat.H*basis_mat)==s.eye(9)
    kraus=[omega*vec.H for vec in v]
    assert kraus[0]==omega*omega.H  # P_3, not just on the input subspace
    assert s.simplify(sum((k.H*k for k in kraus),s.zeros(9)))==s.eye(9)
    sigmas=[]; p=[]; concurrence=[]
    for k in kraus:
        # k*psi = |omega3> (bra_v psi), pure qutrit Bell x remote state
        projected=k*psi
        sigma=s.simplify(projected.T*projected)  # here all vectors are real
        sigmas.append(sigma)
        prob=s.simplify(sigma.trace());p.append(prob)
        if prob==0:
            concurrence.append(None)
        else:
            # remote ket coefficients are projected onto output qutrit Bell
            amplitude=s.simplify((omega.T*projected).T/SQ(prob))
            A=s.Matrix([[amplitude[0],amplitude[1]],[amplitude[2],amplitude[3]]])
            concurrence.append(s.simplify(2*abs(A.det())))
            assert concurrence[-1]==1
            assert s.simplify(sigma*sigma-prob*sigma)==s.zeros(4)
    assert p==[R(1,6),R(1,12),R(1,4),R(1,4),R(1,4),0,0,0,0],p
    assert s.simplify(sum(sigmas,s.zeros(4))-s.eye(4)/4)==s.zeros(4)
    pb=basis(4,0)/SQ(2)+basis(4,3)/SQ(2)
    P_B=pb*pb.T
    assert s.simplify(sigmas[0]-R(1,6)*P_B)==s.zeros(4)
    # Binary Lüders instrument: (P3,I-P3). Failure is outside the pair-only
    # source ontology, although standard CPTP nonselective measurement works.
    P3=kraus[0]; F=s.eye(9)-P3
    assert P3*P3==P3 and F*F==F
    assert s.simplify(P3.H*P3+F.H*F)==s.eye(9)
    other=s.simplify((F*psi).T*(F*psi))
    assert s.simplify(other - (s.eye(4)/4-R(1,6)*P_B))==s.zeros(4)
    failure=s.simplify(other/R(5,6))
    assert s.simplify(failure-(R(3,10)*s.eye(4)-R(1,5)*P_B))==s.zeros(4)
    assert s.simplify((P_B*failure).trace())==R(1,10)
    return {'dimension_A':9,'dimension_B':4,
      'source_Bell_projector_is_exactly_first_Kraus':True,
      'CPTP_completeness_all_9_A_dimensions':True,
      '9_Kraus_outcome_probabilities':[str(q) for q in p],
      '9_Kraus_remote_pair_concurrence':[str(q) if q is not None else None for q in concurrence],
      'all_5_nonzero_reprepared_outcomes_are_pure_Bell3_x_Bell2':True,
      'after_averaging_remote_state':'I_4/4',
      'binary_Luders_fail_remote_density':'(3/10)*I_4 - (1/5)*P_B',
      'binary_Luders_failure_probability':'5/6',
      'binary_Luders_failure_Bell_fidelity':'1/10',
      'limitation':'Extra Kraus outcomes reset/reprepare qutrit Bell pairs: excluded by the source\'s distillation-only rule; fixed (3,2)x(3,2) source input only.'}

def full_qutrit_bell_measurement():
    # |Phi_(t,s)> = sum_j exp(2*pi*i*t*j/3)|j,j+s mod 3>/sqrt(3)
    # Use NumPy only for roots of unity; probability formula is verified against
    # rational expected results for all nine outcomes.
    psi=np.zeros((9,4),complex)
    for a in range(2):
        for b in range(2): psi[3*a+b,2*a+b]=0.5
    betas=[]; probs=[]; conc=[]; rho=np.zeros((4,4),complex)
    for shift in range(3):
      for phase in range(3):
        x=np.zeros(9,complex)
        for j in range(3):
            x[3*j+(j+shift)%3]=np.exp(2j*np.pi*phase*j/3)/np.sqrt(3)
        betas.append(x)
        remote=x.conj()@psi
        pp=float(np.vdot(remote,remote).real)
        expected=1/6 if shift==0 else 1/12
        assert abs(pp-expected)<2e-14
        probs.append(pp)
        A=(remote/np.sqrt(pp)).reshape((2,2))
        C=2*abs(np.linalg.det(A))
        assert abs(C-(1 if shift==0 else 0))<2e-14
        conc.append(C)
        rho+=np.outer(remote,remote.conj())
    B=np.vstack(betas)
    assert np.linalg.norm(B@B.conj().T-np.eye(9))<2e-14
    assert np.linalg.norm(rho-np.eye(4)/4)<2e-14
    assert abs(sum(probs)-1)<1e-14
    return {'qutrit_Bell_basis':9, 'probability_3_shift_zero_each':'1/6',
      'probability_6_nonzero_shift_each':'1/12',
      'total_probability_Bell_pair_remote':'1/2',
      'total_probability_product_remote':'1/2',
      'all_outcomes_remote_average_I4_over_4':True,
      'remote_concurrences':[round(x,12) for x in conc],
      'caveat':'Six outcomes have unentangled remote qubits, outside the source pure-pair-only ontology.'}



def independent_support_rotation_test():
    """Counterexample to universal pair-state preservation for fixed P3.

    Rotate only the first link's qutrit support span{0,1} to span{0,2}
    via the unitary permutation 1<->2. Each input link remains a normalized
    maximally entangled (3,2) state. Keep the same physical P3 join fixed.
    Whether such independently oriented supports belong to a particular
    interpretation of the source's matching rules must be decided explicitly.
    """
    psi=s.zeros(9,4)
    a1={0:0,1:2}
    a2={0:0,1:1}
    for i in range(2):
        for j in range(2): psi[3*a1[i]+a2[j],2*i+j]=R(1,2)
    assert psi.T*psi==s.eye(4)/4
    omega=(basis(9,0)+basis(9,4)+basis(9,8))/SQ(3)
    selected=s.simplify(omega.T*psi)
    prob=s.simplify((selected*selected.T)[0])
    assert prob==R(1,12)
    normalized=s.simplify(selected.T/SQ(prob))
    assert normalized==basis(4,0)
    # Bob qubits are product |0,0>, not a maximally entangled pair.
    A=s.Matrix([[normalized[0],normalized[1]],[normalized[2],normalized[3]]])
    assert A.det()==0
    return {'first_support':'span(|0>,|2>)', 'second_support':'span(|0>,|1>)',
       'each_initial_3x2_link_is_maximally_entangled':True,
       'fixed_P3_outcome_probability':'1/12',
       'normalized_remote_output':'|0,0> (separable product)',
       'remote_pure_state_concurrence':'0',
       'meaning':'The particular fixed P3 selected branch is not pair-preserving for all independently oriented two-dimensional qutrit supports. A covariant support-dependent join or an admissibility restriction is necessary. The source does not settle that operational choice.'}


def record_spectrum():
    # Table I: conditional readout eigenvalues; source fixes NO failure weights.
    count=[79,48,42,3]
    vals=[s.Integer(1),SQ(2),SQ(3),SQ(6)]
    n=24*36*36
    assert sum(count)==172
    a=s.simplify(sum(c*v for c,v in zip(count,vals))/n)
    assert s.simplify(a-(79+48*SQ(2)+42*SQ(3)+3*SQ(6))/31104)==0
    p=R(1,6)
    branch_examples={'M_failure_zero':p*a,'M_failure_same_as_success':a,
                     'M_failure_identity':p*a+(1-p)}
    # All examples are explicit added branch extensions, not predicted couplings.
    return {'conditional_sector_counts':count, 'readout_dimension':n,
            'conditional_alpha_exact':str(a),
            'conditional_alpha_inverse_decimal':str(s.N(1/a,25)),
            'branch_readout_examples':{k: str(s.simplify(v)) for k,v in branch_examples.items()},
            'branch_readout_approximations':{k: float(s.N(v,16)) for k,v in branch_examples.items()},
            'inference':'Only the selected-outcome count is specified. CPTP completion does not uniquely determine the unread or failed-branch readout; unconditional alpha is undefined without a new observable and an extension rule.'}

def main():
    output={
        'scope':'CPTP instrument and 172-record compatibility audit on the original pair input; NO QED beta function derived.',
        'binary_and_pair_preserving_instruments':exact_instruments(),
        'standard_complete_qutrit_Bell_measurement':full_qutrit_bell_measurement(),
        'independently_rotated_support_stress_test':independent_support_rotation_test(),
        '172_record_alpha_status':record_spectrum(),
        'status':'All exact SymPy and finite complex NumPy assertions passed.'
    }
    print(json.dumps(output,indent=2))
if __name__=='__main__': main()