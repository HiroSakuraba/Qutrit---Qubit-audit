# Chain Theory → QED beta program, v18
## Entanglement-reconstructed graph geometry, path contraction, and absence of a selected physical scale

**Date:** 25 October 2026  
**Primary source:** A. D. Wissner-Gross, *A Parameter-Free Candidate for the Fine-Structure Constant from Entangled Qudits* (23 October 2026), `ChainTheoryAlpha.pdf`.  
**Prior work:** v2–v17; especially the v17 gauge-intertwiner classification.  
**Research status:** Mathematical theorems and numerical tests of *virtual-copy Bell-pair networks*, plus explicit **added** dynamical countermodels. No emergence of a physical clock, relativistic dispersion, or the QED beta-function coefficient from the source axiom has been demonstrated.

## Abstract

The original axiom specifies finitely many maximally entangled qudit links, shared endpoints modeled with one virtual copy per link, and an act of distillation implemented by a selected Bell projection. This does not determine an external time parameter, spectral dimension, or continuous momentum scale. It **does**, in the virtual-copy tensor-product interpretation, determine an exact cut-entropy function on any *given* graph. The entropy of a vertex subset equals the sum of `log(m_e)` over its boundary edges, and singleton and pair entropies reconstruct every weighted edge. This gives a precise positive route from entanglement to *combinatorial network geometry*, contingent on taking virtual-copy reduced states seriously.

That result is not the same as deriving an operational spacetime geometry. The paper expressly prohibits estimation from multiple copies and has not shown how an observer obtains entropy values using its stipulated single-copy reads. Moreover, equally admissible three-regular qutrit–qubit link networks have different diameters, girths, and Laplacian spectra; arbitrary extra nearest-neighbor Hamiltonians produce different propagation on the same network. Bell swapping across any path has a dimension-and-join-count amplitude, but its normalized endpoint pair does not retain that count. Therefore none of the tests fixes a unique physical distance-to-energy conversion or QED running coefficient.

## 1. Primary-source inventory (what is and is not assumed)

The source defines a link `(d1,d2)` as the common-support maximally entangled pair `|ω_m⟩`, `m=min(d1,d2)` (Letter Eq. 1). The shared qudit has **one copy per pair**, the source's PEPS-like convention (Letter page 1; Supplement S1). The only designated act is projection onto `P_m=|ω_m⟩⟨ω_m|`, with no failure branch in the ontology (Supplement S1, Corollary S2). Readable types are `(2,2)` and `(2,3)`; the unique charged two-read, one-unread *local object* is the three-arm qutrit hub (Letter Theorem 1). Its normalized readout (Letter Eq. 6) is

\[
\alpha_{\rm cand}=\frac{79+48\sqrt2+42\sqrt3+3\sqrt6}{31104}.
\]

The paper's Supplement S8 expressly states it does **not** construct a function of momentum transfer or a propagator; the identification with `α(0)` is its Photon Conjecture. It also states that statistical expectation-value estimation and fresh repetitions of one named link are unavailable (Supplement S1, Corollaries S4–S5). All external dynamics below are consequently called out as hypotheses, not alleged deductions.

## 2. The cut-entropy theorem (positive result)

**Hypotheses.** Let `G=(V,E)` be a finite *simple* graph. Each edge `e=(i,j)` carries the pure maximally entangled state

\[
|\omega_{m_e}\rangle_{e_i,e_j}=\frac1{\sqrt{m_e}}
\sum_{a=1}^{m_e}|a\rangle_{e_i}|a\rangle_{e_j},\qquad m_e\ge2.
\]

Use the paper's independent virtual copy on each incident edge and form the product virtual-link state

\[
|\Omega_G\rangle=\bigotimes_{e\in E}|\omega_{m_e}\rangle.
\]

For a vertex subset `A⊂V`, place all virtual legs incident to vertices in `A` on one side of a bipartition, and the remaining virtual legs on the other. Write `∂A` for edges connecting `A` to its complement.

**Theorem 18.1 (exact edge-cut entropy).**

\[
\boxed{S(A)=\sum_{e\in\partial A}\log m_e.}
\]

**Proof.** Each pair contained wholly on one side contributes a pure tensor factor and entropy zero. Each pair cut by the bipartition has reduced density matrix `I_{m_e}/m_e` on its Schmidt support, hence entropy `log m_e`. Tensor-product entropy is additive. There are no cross-link correlations because the virtual-link state is a tensor product. QED.

For distinct vertices `i,j`, define the edge weight `w_ij=log m_ij` if an edge exists and `w_ij=0` otherwise. The cut function gives

\[
S(i)=\sum_{k\ne i}w_{ik},\quad
S(j)=\sum_{k\ne j}w_{jk},\quad
S(\{i,j\})=S(i)+S(j)-2w_{ij}.
\]

Therefore:

\[
\boxed{w_{ij}=\frac12\big[S(i)+S(j)-S(\{i,j\})\big]
=\frac12 I(i:j).}
\]

Here `I(i:j)` is the quantum mutual information between all virtual copies at the two named vertices. Thus singleton and pair entropies reconstruct every edge and its Schmidt-rank weight exactly. In particular the adjacency relation `i∼j` is equivalent to `I(i:j)>0`. One may define a shortest-path *combinatorial* distance using any positive edge-length function of these weights, but the axioms do not select the conversion into physical lengths. The edge entropy is *not* a physical length by dimensional analysis.

**Qualified status.** The theorem is exact for the virtual-copy tensor product. The source also calls all copies at one shared end “one qudit,” and does not provide an operational tensor-product subsystem access rule or a many-copy entropy estimator. It is therefore not presently a theorem that an observer inside Chain Theory can measure `S(A)` via permissible distillations. Any quotient/vertex projection identifying the virtual copies can also modify the tensor-product entropy formula and requires a fresh calculation.

### 2.1 Explicit finite-Hilbert-space check

Two qutrit–qubit links with identical rank-two support are represented on `(C³⊗C²)⊗(C³⊗C²)` by `|ω₂⟩⊗|ω₂⟩`. Direct singular-value decomposition of the associated tensor gives entropies in bits:

- `S(Q₁)=S(Q₂)=1`;
- `S(Q₁Q₂)=2`, since those vertex virtual legs share no direct edge;
- `S(Q₁q₁)=0`, since that entire Bell pair is pure and disconnected from its complement;
- `S(q₁q₂)=2`.

All agree with Theorem 18.1 within `<7×10⁻¹⁶` numerical entropy error.

## 3. Global geometry is not selected by identical local charged objects

Construct the 3-regular bipartite graph `B(n,S)` with qutrits `Q_i` and qubits `q_j`, `i,j∈Z/nZ`, and edges `Q_i∼q_{i+s}` for each `s∈S`. Every edge is a `(3,2)` link with Schmidt rank two. Each qutrit has three links and therefore supports the source's *local* three-arm object (upon choosing two arms as reads). Take `n=7`:

- **A:** `S_A={0,1,2}`;
- **B:** `S_B={0,1,3}`.

The graphs have equal numbers of qudits and edges, same local valence and link type, but different global graph invariants:

| invariant | A | B |
|:--|--:|--:|
| vertices | 14 | 14 |
| edges | 21 | 21 |
| degree at each vertex | 3 | 3 |
| shortest cycle | 4 | 6 |
| number of 4-cycles | 7 | 0 |
| diameter | 4 | 3 |
| `Tr(adjacency⁴)` | 266 | 210 |
| first nonzero combinatorial-Laplacian eigenvalue | 0.7530203963 | 1.5857864376 |

The numerical gap is not a physical mass or energy; it is a graph spectrum of a graph Laplacian **chosen as a diagnostic**. Different gaps also produce different heat-kernel return probabilities and apparent finite-graph spectral-dimension diagnostics. This establishes an *underselection* of graph geometry by the local object, not that no geometry can be recovered once a global network is supplied. The entropy reconstruction above succeeds on *both* graphs and on the larger `B(13,{0,1,3})`.

## 4. Bell swapping determines conditional amplitudes, not a time-evolution operator

For `d`-dimensional common-support Bell pairs, write `B_d=I_d/√d` for the coefficient matrix. Contracting two adjacent Bell pairs over a Bell bra gives

\[
B_d B_d^* B_d=\frac1d B_d.
\]

For a path of `n` links joined by `n−1` such contractions,

\[
\boxed{\langle\omega_d|^{\otimes(n-1)}
\big(|\omega_d\rangle^{\otimes n}\big)
=d^{-(n-1)}|\omega_d\rangle_{\rm endpoints}.}
\]

A particular sequence of these projection outcomes would have ordinary Born probability `d^{-2(n-1)}`. However, under the source's always-selected, normalized postselection ontology, **the final state is simply `|ω_d⟩` independent of path length**, and the failure probabilities are not available as ordinary outcome frequencies. Consequently

\[
\mathcal L_{\rm amplitude}(n)=(n-1)\log d
\]

is a *mathematical* bookkeeping variable derived from unnormalized amplitudes, not yet an operational metric, clock, or inverse momentum scale. The script checks the contraction exactly up to floating point for `d=2,3,4`, `n=1,…,5`. Changing the order of disjoint virtual-copy joins does not establish physical time order; this is an associative tensor contraction, as already analyzed in v5.

## 5. Exact pair-state matching invariance does not determine propagating dynamics

For arbitrary `U∈U(d)`, the source's matching action obeys

\[
(U\otimes\overline U)|\omega_d\rangle=|\omega_d\rangle.
\]

This is a basis-relabeling symmetry. In the script, 48 Haar-style random unitary trials for `d=2,3` have maximum state-vector residual `<6×10⁻¹⁶`. Since the pair state is unchanged, the matching action cannot by itself distinguish elapsed time or create a new wave packet at a distant vertex.

To demonstrate what *adding* dynamics would entail, take an **auxiliary one-particle Hilbert space** `C^{|V|}` with a basis indexed by vertices and the graph's adjacency matrix `A_G`. Define the Hermitian local Hamiltonian

\[
H_{\gamma,\Delta}=\gamma A_G+\Delta Z,\quad
Z|Q_i\rangle=+|Q_i\rangle,\quad Z|q_i\rangle=-|q_i\rangle.
\]

It obeys nearest-neighbor hopping and respects translations/permutations preserving the two vertex types. Neither `γ`, `Δ`, nor an evolution parameter is present in the Chain Theory axiom. The `ΔZ` term is a different dynamical law, not just a relabeling of time. With `γ=1`, an added time `t=0.6`, and a walker initially at `Q₀`, the total probability on its three neighbors is

| graph | `Δ=0` | `Δ=1.5` |
|:--|--:|--:|
| A | 0.4818022227 | 0.3453224025 |
| B | 0.5868649943 | 0.4272752440 |

Both evolutions are unitary to numerical precision. They simply illustrate that network adjacency does not itself select *dynamics*. They are **not** known to preserve the original “nothing exists but maximally entangled pairs” ontology during evolution and are not offered as completions satisfying every axiom.

An independent auxiliary example assigns transport matrices to an abstract loop of three virtual-copy frames. The same underlying product-pair state permits `Tr H=2` for identity transports or `Tr H=√2` for a nontrivial `SU(2)` loop; these transports are extra variables, not derived loop data. This reprises the limited holonomy observation of v10.

## 6. Consequence for the QED one-loop coefficient

The source's local record formula is independent of the number of identical three-arm qutrit hubs elsewhere in a network. Thus *if* one postulates a “scale” `μ(G)` from graph size while continuing to assign each local object the same source value, then

\[
\alpha_{\rm local}(\mu)=\alpha_{\rm cand},\qquad
\frac{d\alpha_{\rm local}}{d\log\mu}=0.
\]

This conditional zero is not a prediction of QED running, nor a no-go against every extension. It shows that a simple graph-size relabeling of the unmodified local record formula cannot recover the nonzero one-Dirac-fermion target

\[
\beta(\alpha)=\frac{2}{3\pi}\alpha^2+O(\alpha^3).
\]

To obtain the target independently, a future model must select (i) an operational energy or momentum variable, (ii) propagating charged matter with spin/statistics and the correct light-species count, (iii) a physical two-point function and a legitimate scale transformation, and (iv) the logarithmic response's coefficient rather than inserting known QED integrals. The v2–v17 results establish that simply choosing conventional Dirac dynamics and its heat-kernel normalization is insufficient for an *independent* derivation.

In particular, a dimensionless rescaling `ℓ→cℓ` of a graph-resolution parameter changes the inferred derivative `dα/dℓ` by `1/c`. Without specifying how `ℓ` is related to `log(physical energy)`, its beta-function coefficient cannot be compared to the physical QED value. Unlike a harmless change of physical energy *units* (`μ→cμ`, which leaves `d/dlogμ` unchanged), this is a change in the **mapping of graph resolution to physical scale**.

## 7. Next falsifiable construction

The next useful mathematical experiment is to treat source-allowed reads as a category of operations on finite network configurations and ask whether an intrinsically defined, positive, compositional *transfer map* exists on the observable algebra. Its necessary requirements are:

1. Functoriality/covariance under the paper's independent pair matchings and its shared-virtual-copy identification.
2. No-signaling when interpreted as physical local operations, without silently replacing the paper's selected postselection by a different instrument.
3. A nontrivial one-parameter semigroup (or a discrete-to-continuum limit) acting on physical observables, not just frame labels.
4. A source-determined mapping from graph/entropy resolution to physical momentum, to remove arbitrary time/energy calibration.
5. Charged matter and a Ward identity for the resulting two-point function.

Only after such a transfer map is derived can the QED beta function be a genuine target of independent computation.

## 8. Computational verification

- Script: `chain_beta_extension_v18.py` (NumPy, SciPy, NetworkX, SymPy).
- Machine results: `chain_beta_extension_v18_results.json`.
- Console summary: `chain_beta_extension_v18_stdout.txt`.
- Exact and integer checks: graph edge, degree, girth, cycle, and adjacency counts; symbolic evaluation of the source candidate; analytic tensor identities. The v18 script does not independently re-enumerate the source’s 172 record directions (those checks are retained in earlier suites).
- Numerical checks: SVD entropies on the explicit `(3×2)⊗(3×2)` state, full pair-entropy graph reconstruction on three graphs, 15 swap path/dimension cases, 48 unitary matching cases, unitarity and different transition probabilities of additional Hamiltonians.
- The script's `PASS` certifies **only the assertions it implements**. It neither constructs physical QED from the original axiom nor validates the photon identification conjecture.

### Related external research (context only)

G. M. D'Ariano and P. Perinotti, *Derivation of the Dirac equation from principles of information processing*, Physical Review A **90**, 062106 (2014), DOI `10.1103/PhysRevA.90.062106`, studies quantum cellular automata with **added** unitarity, locality, homogeneity, and isotropy assumptions. It demonstrates a plausible methodological route but does not establish those dynamical postulates as deductions of Wissner-Gross's single Chain Theory axiom.