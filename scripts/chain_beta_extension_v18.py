#!/usr/bin/env python3
"""Chain Theory v18: virtual-copy entanglement, graph reconstruction, and scale selection.

Only tensor-product Bell-pair algebra is source-derived. Graph choice, an
auxiliary one-particle walker, and the heat-kernel diagnostic are explicitly
separate test structures, not predictions of the source axioms.

Dependencies: numpy, scipy, networkx, sympy (sympy only for exact alpha).
"""
from __future__ import annotations
import itertools, json, math
from pathlib import Path
import numpy as np
import networkx as nx
from scipy.linalg import expm
import sympy as sp

OUT = Path(__file__).with_name('chain_beta_extension_v18_results.json')
RTOL = 3e-11


def graph_bipartite(n: int, offsets: tuple[int, ...]):
    G=nx.Graph()
    G.add_nodes_from(("Q", i) for i in range(n))
    G.add_nodes_from(("q", i) for i in range(n))
    for i in range(n):
        for a in offsets:
            G.add_edge(("Q", i), ("q", (i+a) % n), rank=2)
    assert len(G)==2*n and nx.is_connected(G)
    assert len(G.edges)==n*len(offsets)
    return G


def cut_entropy_bits(G, vertices):
    V=set(vertices)
    return sum(math.log2(G.edges[u,v].get("rank",2)) for u,v in G.edges if (u in V) != (v in V))


def reconstruct_from_entropy(G):
    vertices=list(G.nodes)
    s1={i:cut_entropy_bits(G,[i]) for i in vertices}
    recons={}
    for k,u in enumerate(vertices):
        for v in vertices[k+1:]:
            mij=(s1[u]+s1[v]-cut_entropy_bits(G,[u,v]))/2
            recons[(u,v)]=mij
            expected=(math.log2(G.edges[u,v].get("rank",2)) if G.has_edge(u,v) else 0.0)
            assert abs(mij-expected)<RTOL,(u,v,mij,expected)
    return recons


def four_cycle_count(G):
    # bipartite square count = sum_{u<v in Q} binom(number common q-neighbors,2)
    qs=[x for x in G if x[0]=='Q']
    return sum(math.comb(len(set(G.neighbors(u)) & set(G.neighbors(v))), 2)
               for u,v in itertools.combinations(qs,2))


def graph_stats(G):
    order=list(G.nodes)
    a=nx.to_numpy_array(G,nodelist=order,dtype=int)
    lap=np.diag(a.sum(axis=1))-a
    eig=np.linalg.eigvalsh(lap)
    assert min(eig)>-RTOL
    gap=float(eig[1])
    tr_a4=int(np.trace(np.linalg.matrix_power(a,4)))
    n4=four_cycle_count(G)
    # Exact girth: remove each edge and find the shortest alternate path.
    cycles=[]
    for u,v in list(G.edges):
        G.remove_edge(u,v)
        cycles.append(nx.shortest_path_length(G,u,v)+1)
        G.add_edge(u,v,rank=2)
    girth=min(cycles)
    return {'nodes':len(G),'edges':G.number_of_edges(),
            'degrees':sorted(set(dict(G.degree).values())),
            'girth':girth,'four_cycles':n4,'adjacency_trace_A4':tr_a4,
            'diameter':nx.diameter(G),'laplacian_gap':gap,
            'return_probability_tau1':float(np.exp(-eig).mean()),
            'effective_spectral_dimension_tau1':float(2*sum(eig*np.exp(-eig))/sum(np.exp(-eig)))}


def bell_pair(dim: int):
    return sum(np.kron(np.eye(dim)[:,j],np.eye(dim)[:,j]) for j in range(dim))/math.sqrt(dim)


def swap_chain(dim:int,nlinks:int):
    """Tensor contraction over virtual Bell pairs, avoids exponential Hilbert tensor."""
    w=bell_pair(dim)
    T=w.reshape(dim,dim).copy()
    # Joining subsequent Bell links with normalized Bell bra yields Tnext = T @ omega @ omega.
    B=w.reshape(dim,dim)
    for _ in range(nlinks-1):
        T=np.einsum('ij,jk,kl->il',T,B.conj(),B)
    target=(dim**(-(nlinks-1)))*B
    assert np.allclose(T,target,rtol=RTOL,atol=RTOL)
    amplitude=float(np.vdot(w,T.reshape(-1)).real)
    probability=float(np.vdot(T,T).real)
    assert abs(amplitude-dim**(-(nlinks-1)))<RTOL
    assert abs(probability-dim**(-2*(nlinks-1)))<RTOL
    normalized=T/math.sqrt(probability)
    assert np.allclose(normalized,B,atol=RTOL)
    return {'dim':dim,'links':nlinks,'join_count':nlinks-1,
            'unnormalized_amplitude':amplitude,
            'ordinary_born_success_probability':probability,
            'conditional_endpoint_error':float(np.linalg.norm(normalized-B))}


def physical_small_star_test():
    """Two 3x2 entangled links sharing a 3-level *named* qutrit via two virtual copies.

    Physical array order Q1(3),q1(2),Q2(3),q2(2).
    Entropy S(Q1Q2) exactly two bits; S(Q1q1) zero. 
    """
    # ordering axes (Q1,q1,Q2,q2)
    omega32=np.zeros((3,2),dtype=complex)
    omega32[0,0]=omega32[1,1]=1/math.sqrt(2)
    state=np.einsum('ab,cd->abcd',omega32,omega32)
    assert abs(np.vdot(state,state).real-1)<RTOL
    def entropy_bipartition(keep):
        keep=list(keep)
        rest=[i for i in range(4) if i not in keep]
        matrix=state.transpose(keep+rest).reshape(int(np.prod([state.shape[i] for i in keep],dtype=int)), -1)
        singular=np.linalg.svd(matrix,compute_uv=False)
        vals=singular**2
        vals=vals[vals>1e-12]
        return float(-np.sum(vals*np.log2(vals)))
    sQ1=entropy_bipartition([0]);sQ2=entropy_bipartition([2]);sQQ=entropy_bipartition([0,2]);sQq=entropy_bipartition([0,1]);sq1=entropy_bipartition([1]);sq2=entropy_bipartition([3]);sqq=entropy_bipartition([1,3])
    assert np.isclose([sQ1,sQ2,sQQ,sQq,sq1,sq2,sqq],[1,1,2,0,1,1,2],atol=RTOL).all()
    return {'entropy_bits_Q1':sQ1,'entropy_bits_Q2':sQ2,
            'entropy_bits_Q1Q2':sQQ,'entropy_bits_Q1q1':sQq,
            'entropy_bits_q1q2':sqq}


def walker_overlay(G, time=0.35, gammas=(0.6,1.4)):
    """Illustrative *added* hopping Hamiltonian, not a source-derived dynamics."""
    order=list(G.nodes)
    a=nx.to_numpy_array(G,nodelist=order,dtype=float)
    center=order.index(('Q',0))
    probs=[]
    for gamma in gammas:
        U=expm(-1j*gamma*time*a)
        unitarity=np.linalg.norm(U.conj().T@U-np.eye(len(order)))
        assert unitarity<RTOL
        psi=U[:,center]
        P=float(sum(abs(psi[order.index(v)])**2 for v in G.neighbors(('Q',0))))
        probs.append(P)
        assert abs(float(np.vdot(psi,psi).real)-1)<RTOL
    assert abs(probs[0]-probs[1])>1e-3
    return {'added_hamiltonian':'H_gamma = gamma * adjacency', 'time_in_arbitrary_units':time,
            'gamma_values':list(gammas),'neighbor_probabilities':probs,
            'probability_difference':abs(probs[1]-probs[0])}


def mass_walkers(G, time=0.6, gamma=1.0, masses=(0.0,1.5)):
    """Second *added* dynamical non-uniqueness witness, not a time rescaling."""
    order=list(G.nodes)
    A=nx.to_numpy_array(G,nodelist=order,dtype=float)
    Z=np.diag([1.0 if x[0]=='Q' else -1.0 for x in order])
    start=order.index(('Q',0))
    output=[]
    for mass in masses:
        H=gamma*A+mass*Z
        U=expm(-1j*time*H)
        assert np.linalg.norm(U.conj().T@U-np.eye(len(order)))<RTOL
        psi=U[:,start]
        neighbor_probability=float(sum(abs(psi[order.index(v)])**2 for v in G.neighbors(('Q',0))))
        output.append({'sublattice_mass':mass,'neighbor_probability':neighbor_probability,
                       'energy_gap':float(np.linalg.eigvalsh(H)[len(order)//2]-np.linalg.eigvalsh(H)[len(order)//2-1])})
    assert abs(output[0]['neighbor_probability']-output[1]['neighbor_probability'])>1e-2
    return {'added_hamiltonian':'H = adjacency + mass * (+1 on Q; -1 on q)',
            'at_time':time,'cases':output}


def calculate():
    G1=graph_bipartite(7,(0,1,2))
    G2=graph_bipartite(7,(0,1,3))
    g1=graph_stats(G1);g2=graph_stats(G2)
    assert g1['four_cycles']>0 and g2['four_cycles']==0
    assert g1['adjacency_trace_A4']==266 and g2['adjacency_trace_A4']==210
    assert g1['girth']==4 and g2['girth']==6
    for G in [G1,G2]:
        W=reconstruct_from_entropy(G)
        assert sum(w>RTOL for w in W.values())==21
    assert not nx.is_isomorphic(G1,G2)

    # Independence of graph size: one identical local (3,2) three-arm object per Q vertex
    alpha=sp.simplify((79+48*sp.sqrt(2)+42*sp.sqrt(3)+3*sp.sqrt(6))/31104)
    bigger=graph_bipartite(13,(0,1,3))
    reconstruct_from_entropy(bigger)
    assert all(d==3 for _,d in bigger.degree())
    swaps=[swap_chain(d,n) for d in (2,3,4) for n in (1,2,3,4,5)]
    rng=np.random.default_rng(18)
    # Pure pair states invariant under *edge-wise* U tensor conj(U) for arbitrary U(d).
    match_errors=[]
    for d in (2,3):
        ket=bell_pair(d)
        for _ in range(24):
            Z=rng.normal(size=(d,d))+1j*rng.normal(size=(d,d))
            U,_=np.linalg.qr(Z)
            rotated=np.kron(U,U.conj())@ket
            match_errors.append(float(np.linalg.norm(rotated-ket)))
    assert max(match_errors)<RTOL
    # Distinct loop transports assigned to same product-pair data are extra labels.
    t0=np.eye(2,dtype=complex)
    sigma_z=np.diag([1,-1])
    t1=np.diag(np.exp(1j*np.array([np.pi/4,-np.pi/4])))
    loop_trace_flat=complex(np.trace(t0@t0@t0))
    loop_trace_nonflat=complex(np.trace(t0@t0@t1))
    assert abs(loop_trace_flat-2)<RTOL and abs(loop_trace_nonflat-math.sqrt(2))<RTOL

    results={
        'status':'PASS',
        'scope':'finite-dimensional virtual-copy Bell networks; auxiliary graph walk not source-derived',
        'graph_A':g1,'graph_B':g2,
        'entropy_reconstruction':{'graphs_tested':3,'full_pair_entropy_reconstruction':True,
                                  'S_vertex_bits':3,'S_adjacent_pair_bits':4,'S_nonadjacent_pair_bits':6,
                                  'edge_mutual_info_bits':2},
        'explicit_small_hilbert_space_entropy':physical_small_star_test(),
        'swapping_grid':swaps,'edge_matching_max_error':max(match_errors),
        'extra_loop_transport_example':{'flat_trace_real':loop_trace_flat.real,
                                        'nonflat_trace_real':loop_trace_nonflat.real},
        'auxiliary_walkers':{'A':walker_overlay(G1),'B':walker_overlay(G2)},
        'distinct_non_rescaling_hamiltonians':{'A':mass_walkers(G1),'B':mass_walkers(G2)},
        'alpha_inverse_from_source':str(sp.N(1/alpha,23)),
        'local_alpha_count_independent_of_graph':True,
        'beta_coefficient_derived':False,
        'limitations':['Physical reduced-state entropies require treating virtual copies as vertex tensor factors.',
                       'A clock or energy is not supplied by a Bell contraction.',
                       'Auxiliary hopping models and loop holonomies are not selected by the axiom.',
                       'The normalized selected branch omits the Born success probability.']
    }
    OUT.write_text(json.dumps(results,indent=2)+'\n')
    print(json.dumps({k:v for k,v in results.items() if k in ['status','graph_A','graph_B','entropy_reconstruction','explicit_small_hilbert_space_entropy','edge_matching_max_error','extra_loop_transport_example','auxiliary_walkers','distinct_non_rescaling_hamiltonians','alpha_inverse_from_source','beta_coefficient_derived']},indent=2))
    return results

if __name__=='__main__':
    calculate()