#!/usr/bin/env python3
"""Chain Theory QED beta v5: exact operational-completion tests.

SOURCE: A. D. Wissner-Gross (2026), ChainTheoryAlpha.pdf, §§ S1, S7, S8.
The finite pair and distillation formulae are source-derived. The possibility
of choosing whether a local operation occurs, standard subsystem causality,
and complete quantum instruments are EXTERNAL operational assumptions.

Dependencies: Python 3.10+, NumPy, SymPy.
No beta function or physical spacetime is claimed to emerge.
"""
from __future__ import annotations

import itertools
import json
from fractions import Fraction

import numpy as np
import sympy as sp

TOL = 2e-12


def exact_qutrit_qubit_join():
    # The qutrit-qubit link occupies two of the qutrit's three basis states.
    # Contraction: ( <omega_3|_{A1,A2} )
    #              ( |omega_2>_{A1,B1}|omega_2>_{A2,B2} ).
    # The unnormalized remote state is (|00>+|11>)/(2*sqrt(3)),
    # equivalently |omega_2>/sqrt(6).
    ket_b = sp.Matrix([sp.sqrt(sp.Rational(1, 2)), 0, 0,
                       sp.sqrt(sp.Rational(1, 2))])
    remote_vec = ket_b / sp.sqrt(6)
    p = sp.simplify((remote_vec.T * remote_vec)[0])
    assert p == sp.Rational(1, 6)
    post = sp.simplify((remote_vec*remote_vec.T) / p)
    before = sp.eye(4)/4
    diff = post-before
    spectrum = sorted(sp.Matrix(diff).eigenvals().items(), key=lambda kv: float(kv[0]))
    dtr = sp.simplify(sp.Rational(1, 2)*sum(abs(lam)*mult for lam, mult in spectrum))
    assert dtr == sp.Rational(3,4)
    assert (post*post-post) == sp.zeros(4)

    # A trace-preserving, nonselective Bell measurement on Alice (i.e. P, I-P)
    # leaves the remote marginal equal to before.
    unnormalized_success = p*post
    unnormalized_failure = sp.simplify(before-unnormalized_success)
    assert unnormalized_failure.trace()==sp.Rational(5,6)
    failure = sp.simplify(unnormalized_failure / sp.Rational(5,6))
    assert sp.simplify(unnormalized_success + unnormalized_failure - before)==sp.zeros(4)
    assert sp.simplify((post*before).trace())==sp.Rational(1,4)
    assert sp.simplify((post*failure).trace())==sp.Rational(1,10)
    assert sp.simplify((post*post).trace())==1

    # If Alice can choose to enforce only the normalized Bell branch, Bob can
    # distinguish it from 'do nothing' by a Bell projector on B1,B2.
    return {
        'link_type': '(3,2) crossed with (3,2)',
        'joint_success_probability_under_Born': str(p),
        'remote_before': sp.sstr(before),
        'remote_after_conditioned': sp.sstr(post),
        'trace_distance_conditional_vs_unconditioned': str(dtr),
        'bob_bell_test_probability_if_no_join': '1/4',
        'bob_bell_test_probability_if_conditioned_join': '1',
        'bob_bell_test_probability_in_other_measurement_branch': '1/10',
        'bob_bell_test_probability_if_unconditioned_complete_measurement': '1/4',
        'success_plus_failure_identity': True,
        'interpretation': 'Operational signalling only if the normalized branch is freely selectable and deterministic with independent A/B subsystems; neither locality nor freedom of choice is defined by source.'
    }


def bell(d: int):
    return np.eye(d, dtype=np.complex128)/np.sqrt(d)


def contraction_given_order(d: int, n: int, order: tuple[int,...]):
    """n independent Bell links, join neighbouring internal copies.

    Each link j is a tensor with endpoints (2j,2j+1). Join j contracts
    endpoints (2j+1,2j+2). Their tensor factors are disjoint and all
    join bras commute; contraction order should be immaterial.
    """
    state = bell(d)
    for j in range(1, n):
        state = np.tensordot(state, bell(d), axes=0)
    labels = list(range(2*n))
    for j in order:
        left,right = 2*j+1,2*j+2
        i1,i2=labels.index(left),labels.index(right)
        state=np.tensordot(bell(d).conj(),state,axes=([0,1],[i1,i2]))
        labels=[x for x in labels if x not in (left,right)]
    assert set(labels)=={0,2*n-1}, labels
    if labels != [0,2*n-1]:
        state=state.T
    return state


def graph_history_commutation():
    tests=[]
    for d,n in ((2,3),(2,4),(2,5),(3,3),(3,4),(3,5),(4,3),(4,4)):
        expected = bell(d)/d**(n-1)
        perms=list(itertools.permutations(range(n-1)))
        worst=0.0
        for order in perms:
            x=contraction_given_order(d,n,order)
            worst=max(worst,float(np.linalg.norm(x-expected)))
        assert worst<TOL
        tests.append({'d':d,'links':n,'number_of_orderings':len(perms),
                      'all_orderings_yield_same_unnormalized_endpoint_state':True,
                      'max_Frobenius_error':worst,
                      'amplitude_factor':f'1/{d**(n-1)}',
                      'joint_Born_branch_probability':f'1/{d**(2*(n-1))}'})
    return {'results':tests,'interpretation':
      'For disjoint virtual-copy Bell contractions, a final tensor-network state provides no intrinsic order of the joins; this does not prove all possible overlapping read operations commute.'}


def matching_invariance():
    # The (U tensor conjugate(U)) matchings leave a maximally entangled
    # state unchanged: this is a basis symmetry, not an evolution clock.
    cases=[]
    for d in (2,3,4):
        for theta in (0.19, 0.51, 1.17):
            rng = np.random.default_rng(d*300 + round(theta*100))
            z=rng.normal(size=(d,d)) + 1j*rng.normal(size=(d,d))
            q,r=np.linalg.qr(z)
            phases=np.diag(r)
            q=q * np.conj(phases/np.abs(phases))[None,:]
            om=bell(d).reshape(-1)
            changed=np.kron(q,q.conj())@om
            err=float(np.linalg.norm(changed-om))
            assert err < TOL
            cases.append({'d':d,'test_seed':d*300+round(theta*100),'residual':err})
    return {'cases':cases,'interpretation':
            'Matchings U tensor U-bar are gauge/basis transformations on the Bell state, not a nontrivial propagation law.'}


def finite_kinematics_no_canonical_time():
    # Demonstrate a zero-Hamiltonian completion and a (local Bell stabilizer)
    # Hamiltonian acting on copies. Both preserve every input product of Bell
    # link kets, yet their action on arbitrary nonaxiomatic disturbances differs.
    d=2
    p=np.outer(bell(d).reshape(-1), bell(d).reshape(-1).conj())
    i4=np.eye(4)
    h0=np.zeros((4,4))
    h1=i4-p
    assert np.linalg.norm(h1@bell(d).reshape(-1))<TOL
    # e^{-it(I-P)} = P + exp(-it)(I-P)
    t=0.7
    u0=i4
    u1=p+np.exp(-1j*t)*(i4-p)
    assert np.linalg.norm(u1.conj().T@u1-i4)<TOL
    assert np.linalg.norm(u1@bell(d).reshape(-1)-bell(d).reshape(-1))<TOL
    assert np.linalg.norm(u0@bell(d).reshape(-1)-bell(d).reshape(-1))<TOL
    diff=float(np.linalg.norm(u1-u0))
    assert diff>0.1
    return {'two_bell_state_preserving_unitary_extensions_are_distinct':True,
            'norm_difference':diff,
            'caveat':'Hamiltonian and off-axiom disturbances are ADDITIONS, not source-derived states; on the axiom pair state these completions are indistinguishable.'}


def main():
    out={
        'scope':'Exact representation/operational-completion tests; not a QED beta calculation.',
        'source':'Wissner-Gross (2026), ChainTheoryAlpha.pdf, esp. S1 and S7.',
        'exact_qutrit_qubit_operational_probe':exact_qutrit_qubit_join(),
        'all_join_orderings':graph_history_commutation(),
        'matched_basis_invariance':matching_invariance(),
        'nonunique_extended_hamiltonians':finite_kinematics_no_canonical_time(),
        'status':'All finite analytic and numeric checks passed. Conditional locality/choice assumption remains external.'
    }
    print(json.dumps(out,indent=2))


if __name__=='__main__':
    main()