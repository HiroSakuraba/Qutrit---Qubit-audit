#!/usr/bin/env python3
"""Chain Theory QED beta v4: structural checks, not a derivation of QED.

Source basis: Wissner-Gross (2026), A Parameter-Free Candidate for the
Fine-Structure Constant from Entangled Qudits, especially S1, S3, S7, S8.

DEPENDENCIES: Python >=3.10, numpy, networkx.
All additional spacetime/lattice/time rules are explicitly extra axioms.
"""
from __future__ import annotations

import json
import math
from fractions import Fraction

import networkx as nx
import numpy as np

ATOL = 2e-12


def bell(d: int) -> np.ndarray:
    return np.eye(d, dtype=np.complex128)/math.sqrt(d)


def swapped_amplitude(big: int, small: int) -> tuple[np.ndarray,float,float]:
    """Contract (Q1,q1)(Q2,q2) with <omega_big| on (Q1,Q2).

    Both (Q,q) source states are max entangled on support min(big,small).
    Output is an *unnormalized* amplitude tensor on q1,q2.
    """
    assert big >= small >= 2
    omega = np.zeros((big,small),dtype=np.complex128)
    omega[np.arange(small),np.arange(small)] = 1/math.sqrt(small)
    psi = np.einsum('ab,cd->abcd',omega,omega)
    contracted = np.einsum('ac,abcd->bd',np.conjugate(bell(big)),psi)
    expected = bell(small)/math.sqrt(big*small)
    error = float(np.linalg.norm(contracted-expected))
    p = float(np.sum(np.abs(contracted)**2))
    assert error < ATOL
    assert abs(p - 1/(big*small)) < ATOL
    return contracted,p,error


def source_projection_checks() -> dict:
    results=[]
    for big,small in ((2,2),(3,2),(3,3),(4,2),(5,3)):
        out,p,error=swapped_amplitude(big,small)
        results.append({'type':[big,small], 'conditional_success_probability':p,
                        'exact_success_probability':str(Fraction(1,big*small)),
                        'unnormalized_output_error':error})
    ket=bell(2).reshape(-1)
    projector=np.outer(ket,ket.conj())
    own_read_err=float(np.linalg.norm(projector@ket-ket))
    assert own_read_err<ATOL
    # Overlapping pair-state inputs are not all left invariant by a new read.
    return {'existing_link_read_invariance':own_read_err,
            'joins':results}


def overlapping_read_order_test() -> dict:
    """Three ordinary tensor factors; proposed event-conflict diagnostic only.

    Source sharing uses one virtual factor per incident link. Identifying two
    projectors as acting on the SAME middle tensor factor is additional.
    """
    ret=[]
    for d in (2,3,4):
        ket=bell(d).reshape(-1)
        p=np.outer(ket,ket.conj())
        a=np.kron(p,np.eye(d))
        b=np.kron(np.eye(d),p)
        comm=a@b-b@a
        observed=float(np.linalg.norm(comm)**2)
        expected=2*(d*d-1)/(d**3)
        sandwich=float(np.linalg.norm(a@b@a-a/(d*d)))
        assert abs(observed-expected)<ATOL
        assert sandwich<ATOL
        ret.append({'d':d,'overlap_projector_commutator_HS_squared':observed,
                    'exact_formula':f'2*(d^2-1)/d^3 = {Fraction(2*(d*d-1),d**3)}',
                    'ABA_sandwich_error':sandwich})
    return {'overlapping_projectors':ret,
            'limitation':'Tests ordinary overlapping factors; the source uses distinct virtual copies for links at a shared node. Noncommutation does not by itself define a time orientation.'}


def finite_chain_distillation() -> dict:
    """Exact contractions of n Bell links with n-1 Bell joins.

    No probability calculus is ascribed to the source ontology: we report
    ordinary Hilbert-space norm factors for comparison.
    """
    rows=[]
    for d in (2,3,4):
        link=bell(d)
        state=link.copy()
        for n in range(1,6):
            expected=bell(d)/(d**(n-1))
            err=float(np.linalg.norm(state-expected))
            ordinary_branch_p=float(np.linalg.norm(state)**2)
            assert err < ATOL
            assert abs(ordinary_branch_p-d**(-2*(n-1)))<ATOL
            rows.append({'d':d,'links':n,'joins':n-1,
                         'unnormalized_amplitude_factor':float(d**(-(n-1))),
                         'ordinary_Born_success_product':ordinary_branch_p,
                         'amplitude_error':err})
            state=np.einsum('ab,cd,bc->ad',state,link,link.conj())
    return {'ordinary_Hilbert_contraction_of_chains':rows,
            'interpretation':'A path-length amplitude exists before branch normalization; the normalized endpoint Bell state does not determine number of joins.'}


def postselection_signalling_if_promoted_to_deterministic() -> dict:
    """A warning for *expanded* state spaces, not a contradiction within source axiom.

    A,B share a Bell pair. C is independently prepared |c>. Normalized
    Bell projection on B,C deterministically would set A to |c>.
    This is a conditional channel, not a CPTP operation.
    """
    b=bell(2)
    out=[]
    for c in (0,1):
        ket=np.zeros(2,dtype=np.complex128);ket[c]=1
        psi=np.einsum('ab,c->abc',b,ket)
        # <omega|_BC, leftover A
        contracted=np.einsum('bc,abc->a',b.conj(),psi)
        p=float(np.vdot(contracted,contracted).real)
        state=contracted/math.sqrt(p)
        expect=np.zeros(2,dtype=np.complex128);expect[c]=1
        assert np.linalg.norm(state-expect)<ATOL
        assert abs(p-0.25)<ATOL
        out.append({'remote_input_C':int(c), 'output_A_probabilities':np.abs(state).tolist(),
                    'branch_probability':p})
    return {'warning':'Only applies if arbitrary state preparation and deterministic normalized projection are allowed; those are not guaranteed by the source axiom.', 'cases':out}


def circular_bipartite_graph(n:int, offsets:tuple[int,...]) -> nx.Graph:
    g=nx.Graph()
    for i in range(n):
        g.add_node(('Q',i),dimension=3)
        g.add_node(('q',i),dimension=2)
        for o in offsets:
            g.add_edge(('Q',i),('q',(i+o)%n),link_type=(3,2))
    return g


def girth(g:nx.Graph) -> int:
    minimum=math.inf
    # Given each edge (u,v), shortest path with that edge removed closes a cycle.
    for u,v in list(g.edges()):
        h=g.copy();h.remove_edge(u,v)
        if nx.has_path(h,u,v):
            minimum=min(minimum,nx.shortest_path_length(h,u,v)+1)
    return int(minimum) if minimum!=math.inf else -1


def lattice_nonuniqueness() -> dict:
    # Two simple connected 3-regular bipartite graphs, same number of nodes,
    # all edges readable (3,2), each qutrit shares three links.
    configs=[(0,1,2),(0,1,3)]
    rows=[]
    for offsets in configs:
        g=circular_bipartite_graph(7,offsets)
        assert nx.is_connected(g)
        assert nx.is_bipartite(g)
        assert all(g.degree(v)==3 for v in g.nodes())
        assert len(g.nodes())==14 and len(g.edges())==21
        assert set(g[u][v]['link_type'] for u,v in g.edges())=={(3,2)}
        adjacency=nx.to_numpy_array(g,nodelist=sorted(g.nodes()),dtype=float)
        spectral=np.sort(np.linalg.eigvalsh(adjacency))
        qutrits=[('Q',i) for i in range(7)]
        fourcycles=sum(math.comb(len(set(g.neighbors(a))&set(g.neighbors(b))),2)
                       for j,a in enumerate(qutrits) for b in qutrits[j+1:])
        rows.append({'offsets':list(offsets), 'vertices':14,'edges':21,
                     'distinct_4_cycles':fourcycles,
                     'qutrit_degree':3,'qubit_degree':3,'girth':girth(g),
                     'trace_adjacency_4':int(round(np.trace(np.linalg.matrix_power(adjacency.astype(int),4)))) ,
                     'max_adjacency_eigenvalue':float(max(spectral))})
    assert [r['girth'] for r in rows]==[4,6]
    assert [r['distinct_4_cycles'] for r in rows]==[7,0]
    assert rows[0]['trace_adjacency_4'] != rows[1]['trace_adjacency_4']
    return {'distinct_graphs_with_identical_local_3_to_2_star':rows,
            'note':'Per-edge copies represent shared qudits following the paper; this does not prove a unique physical Hilbert-space identification of those copies.'}


def coin(m:float) -> np.ndarray:
    sx=np.array([[0,1],[1,0]],dtype=np.complex128)
    return math.cos(m)*np.eye(2)-1j*math.sin(m)*sx


def step_on_states(psi:np.ndarray,m:float) -> np.ndarray:
    """Local coined walk: coin then spin-dependent nearest-neighbor shifts.

    psi[n,0] shifts to n+1 and psi[n,1] shifts to n-1.
    An explicit 1+1D spatial lattice, clock and coin are extra assumptions.
    """
    new=psi@coin(m).T
    out=np.zeros_like(new)
    out[:,0]=np.roll(new[:,0],1)
    out[:,1]=np.roll(new[:,1],-1)
    return out


def full_step(n:int,m:float) -> np.ndarray:
    total=2*n
    u=np.zeros((total,total),dtype=np.complex128)
    for j in range(total):
        psi=np.zeros((n,2),dtype=np.complex128)
        psi.reshape(-1)[j]=1
        u[:,j]=step_on_states(psi,m).reshape(-1)
    return u


def lattice_translation(n:int) -> np.ndarray:
    t=np.zeros((2*n,2*n),dtype=np.complex128)
    for x in range(n):
        for s in range(2):
            t[2*((x+1)%n)+s,2*x+s]=1
    return t


def lattice_parity(n:int) -> np.ndarray:
    p=np.zeros((2*n,2*n),dtype=np.complex128)
    for x in range(n):
        for s in range(2):
            p[2*((-x)%n)+(1-s),2*x+s]=1
    return p


def local_walk_nonuniqueness() -> dict:
    n=25
    eye=np.eye(2*n,dtype=np.complex128)
    t=lattice_translation(n)
    p=lattice_parity(n)
    rows=[];steps=[]
    for m in (0.,0.21,0.61):
        u=full_step(n,m)
        unitary=float(np.linalg.norm(u.conj().T@u-eye))
        homogeneity=float(np.linalg.norm(u@t-t@u))
        parity=float(np.linalg.norm(u@p-p@u))
        assert max(unitary,homogeneity,parity)<ATOL
        # For a localized initial one-component state, support after step k
        # lies exactly inside the radius-k graph-distance light cone.
        psi=np.zeros((n,2),dtype=np.complex128);psi[0,0]=1
        violations=[];prob_pairs=[]
        for k in range(1,7):
            psi=step_on_states(psi,m)
            distances=np.minimum(np.arange(n),n-np.arange(n))
            tail=float(np.sum(np.abs(psi[distances>k])**2))
            assert tail < ATOL
            violations.append(tail)
            if k==1:
                prob_pairs=[float(np.sum(abs(psi[x])**2)) for x in (1,n-1)]
        rows.append({'m':m,'unitarity_defect':unitary,
                     'translation_covariance_defect':homogeneity,
                     'parity_covariance_defect':parity,
                     'largest_superluminal_tail':max(violations),
                     'step1_probability_at_plusminus_one':prob_pairs})
        steps.append(u)
    inequivalent=[float(np.linalg.norm(steps[0]-steps[i])) for i in (1,2)]
    assert min(inequivalent)>1
    return {'properties':'unitary, radius 1, translation-invariant, parity-symmetric',
            'all_work_even_with_different_coin_parameters':rows,
            'difference_from_massless_update':inequivalent,
            'note':'Coin angle controls dynamics and is not selected by maximally entangled pairs.'}


def blocking_nonuniqueness() -> dict:
    ret=[]
    for n in (6,8,10):
        a=frozenset(tuple(sorted((2*i,(2*i+1)%n))) for i in range(n//2))
        b=frozenset(tuple(sorted(((2*i+1)%n,(2*i+2)%n))) for i in range(n//2))
        rotated=frozenset(tuple(sorted(((x+1)%n,(y+1)%n))) for x,y in a)
        assert a!=b and rotated==b
        ret.append({'ring_sites':n,'translation_related_nearest_neighbor_blockings':2,
                    'individual_blocking_invariant_under_single_shift':False})
    return {'blockings':ret,'note':'This addresses pairwise decimation only, not every possible RG construction.'}


def qudit_charges() -> dict:
    # Source link transforms as (3,2)_{1/6}. Identifying Q_em=T3+Y is extra.
    q=np.kron(np.eye(3),np.diag([2/3,-1/3]))
    charge_square=float(np.trace(q@q).real)
    assert abs(charge_square-5/3)<ATOL
    return {'source_representation':'(3,2)_{1/6}',
            'additional_physics':'Q_em=T3+Y',
            'Tr(Q_em^2)':str(Fraction(5,3)),
            'numeric_trace':charge_square,
            'beta_function_from_original_axiom':'undetermined'}


def main() -> None:
    findings={'study':'Chain Theory beta extension v4; no beta-function derivation claimed',
              'source_projection':source_projection_checks(),
              'chain_distillation':finite_chain_distillation(),
              'overlapping_projection_order':overlapping_read_order_test(),
              'postselection_conditional_signalling_warning':postselection_signalling_if_promoted_to_deterministic(),
              'graph_topology':lattice_nonuniqueness(),
              'unitary_local_walk':local_walk_nonuniqueness(),
              'RG_blocking':blocking_nonuniqueness(),
              'charge_trace':qudit_charges()}
    print(json.dumps(findings,indent=2))
    print('ALL v4 ALGEBRAIC AND FINITE-SIZE TESTS PASSED')


if __name__=='__main__':main()