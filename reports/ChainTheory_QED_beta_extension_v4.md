# Chain Theory and the first QED beta-function coefficient — v4

**24 October 2026 · Operational geometry, entanglement swapping, and causal dynamics**  
**Research status:** Exact finite-dimensional calculations and counterexamples completed. **No independent derivation of the QED beta function from the original qudit axiom has been obtained.**  
**Reproduction:** `python chain_beta_extension_v4.py > chain_beta_extension_v4_results.txt` (Python 3.10+, NumPy, NetworkX).

## Abstract

This note advances the open question from v2–v3: can the shared-qudit and distillation definitions in A. D. Wissner-Gross, *A Parameter-Free Candidate for the Fine-Structure Constant from Entangled Qudits* (23 October 2026), generate a preferred notion of locality and causal propagation without separately postulating a quantum cellular automaton? We prove six limited but useful results.

1. An on-link projection leaves its maximally entangled link unchanged. A Bell projection joining *different* links reorganizes entanglement but has a nonunit normalized projection step.
2. Joining two `(3,2)` links at the qutrit copies has unnormalized output `|ω₂⟩/√6`, and the usual laboratory probability for that specific outcome is `1/6`. This independently reproduces the normalization explicitly used in the source's Supplement §S7. For a chain of `n` equal-dimensional Bell links, the unnormalized endpoint amplitude is `d^{-(n−1)}|ω_d⟩`.
3. On ordinary overlapping tensor factors, successive Bell projections have a nonzero commutator and satisfy a simple exact sandwich identity. This creates an order-sensitive operation pair, but does not supply a time arrow or identify the source's virtual shared copies.
4. Two non-isomorphic 14-vertex, 21-edge, degree-three bipartite networks can meet the identical qutrit–qubit local template, while having girth 4 and 6. The original local link selection does not single out even a global combinatorial geometry.
5. With additional assumptions of one-dimensional lattice space and discrete unitary time, a continuous family of distinct quantum walks shares exact finite-speed causality, translation invariance, and parity symmetry. These constraints do not select the kinetic parameter.
6. Even with a chosen ring geometry, the two simplest disjoint nearest-neighbor blocking choices are exchanged by a one-site translation. No individual such blocking is translation-invariant.

These do **not** rule out a *future strengthened theory*. They separate what is fixed by the source from what would have to be derived or explicitly added.

## 1. Source versus added structure

### 1.1 Statements supported by the uploaded source

**[SOURCE]** The model's only axiom is that everything consists of finitely many maximally entangled pairs of qudits, some sharing named qudits, known only by distillation. A link of type `(d₁,d₂)` uses a maximally entangled state on the common support of dimension `m=min(d₁,d₂)` (main Letter pp. 1–2; Supplement §S1). The source treats each incident link of a shared qudit as having its own copy and defines how these copies compose; this is a network/PEPS-style bookkeeping rule, *not* the usual assertion that one physical `d`-dimensional factor is independently maximally entangled with several separate neighbors.

**[SOURCE]** Distillation is a rank-one projection onto a pair state. The model excludes its complementary outcome from the ontology, rather than representing the full deterministic Bell measurement and its corrections. Supplement §S1, Eq. (S1), explicitly gives entanglement-swapping amplitude `1/d` for equal dimensions; Supplement §S7, Eq. (S2), uses an overall `√6` normalization for a qutrit–qubit join and says the ordinary laboratory success probability is `1/6`.

**[SOURCE]** The source states that the two insertions of its two-point record are not time ordered. Supplement §S8 concedes that the model constructs no momentum-dependent propagator, charged loop, loop momentum, or vertex. Its electromagnetic identification remains a Photon Conjecture. It therefore fixes a static scalar `α_cand(0)`, not a physical time evolution or a renormalization-group flow.

**[ADDED]** When we inspect ordinary Born probabilities, lattice time, spatial graphs, the unitary coin, overlapping projectors on one physical middle factor, or a deterministic normalized postselected channel, those are **external** physical concepts or mathematical probes. None is claimed to be forced by the source axiom.

### 1.2 Specific original QED target

For one relativistic four-component Dirac field of electric charge one in 3+1 spacetime dimensions,

\[
\beta(e)=\mu\frac{de}{d\mu}=\frac{e^3}{12\pi^2}+O(e^5),\qquad
\beta(\alpha)=\frac{2}{3\pi}\alpha^2+O(\alpha^3).
\]

The source's representation `(3,2)_{1/6}` gives `Tr(Q_em²)=5/3` **after** adding `Q_em=T₃+Y`; this is charge bookkeeping, not a spin-statistics or loop normalization theorem. The original axiom still has no `μ` derivative.

## 2. Exact projections and a candidate operational path weight

### Proposition 1: A read of its own pair does not propagate information

Let

\[
|\omega_d\rangle=\frac1{\sqrt d}\sum_{j=0}^{d-1}|j\rangle|j\rangle,
\quad P_d=|\omega_d\rangle\langle\omega_d|.
\]

Then

\[
P_d|\omega_d\rangle=|\omega_d\rangle.
\]

**Proof.** The vector has norm one, so `P_d|ω_d⟩=|ω_d⟩⟨ω_d|ω_d⟩=|ω_d⟩`. A projection onto an already existing pair therefore leaves its state invariant. □

The source's nontrivial network act is **joining different incident link copies**, not repeatedly projecting an already existing link in a way that automatically generates a kinetic term.

### Proposition 2: General unequal-support entanglement-swapping amplitude

Let `D ≥ m ≥ 2`, and take two `(D,m)` links with D-dimensional ends `Q₁,Q₂` and m-dimensional ends `q₁,q₂`:

\[
|\omega_m\rangle_{Q_1q_1}\otimes|\omega_m\rangle_{Q_2q_2}
=\frac1m\sum_{i,j=0}^{m-1}|i\rangle_{Q_1}|i\rangle_{q_1}|j\rangle_{Q_2}|j\rangle_{q_2}.
\]

Join the D-dimensional copies by taking the Bell bra `⟨ω_D|` on `Q₁Q₂`:

\[
\begin{aligned}
{}_{Q_1Q_2}\langle\omega_D|\Bigl(|\omega_m\rangle_{Q_1q_1}|\omega_m\rangle_{Q_2q_2}\Bigr)
&=\frac1{m\sqrt D}\sum_{i=0}^{m-1}|i\rangle_{q_1}|i\rangle_{q_2}\\
&=\boxed{\frac1{\sqrt{mD}}|\omega_m\rangle_{q_1q_2}}.
\end{aligned}
\]

The ordinary Born success probability of **that outcome** is consequently

\[
\boxed{p=\frac1{mD}}.
\]

At `D=3`, `m=2`, the output is precisely `|ω₂⟩/√6`, `p=1/6`; restoring unit norm multiplies by `√6`, as in the source's §S7. At `D=m=d` it reduces to the source's Eq. (S1): amplitude `1/d`, probability `1/d²`.

**[TESTED]** Direct four-index tensor contractions verify the formula for `(2,2)`, `(3,2)`, `(3,3)`, `(4,2)`, and `(5,3)`, with residuals below `5 × 10⁻¹⁶`.

### Proposition 3: Chain contraction tracks the number of joins, but not a clock

For `n` equal-dimensional Bell links in a chain, project each adjacent pair of internal link ends with the normalized Bell bra. Iterating Proposition 2 yields

\[
\boxed{\mathcal A_n=d^{-(n-1)}|\omega_d\rangle},\qquad
\boxed{p_n=d^{-2(n-1)}}.
\]

**Proof by induction.** At `n=1`, `𝒜₁=|ω_d⟩`. Appending an independent Bell link and contracting the connecting ends multiplies the remaining Bell amplitude by `1/d`, so `𝒜ₙ₊₁=(1/d)𝒜ₙ`. □

**[TESTED]** Exact contraction formulas evaluated numerically for `d=2,3,4` and `1≤n≤5`, with all norm and probability residuals `<2×10⁻¹²`.

A potential *operational path-weight functional* is

\[
L_d(n):=-\ln\|\mathcal A_n\|=(n-1)\ln d.
\]

But this is **not** a derived physical distance or duration. The source makes the selected branch an accepted act and normalizes its output. After normalization the endpoint state is always the same `|ω_d⟩` for this simple chain; without keeping track of the pre-normalization amplitudes, `n` is not recoverable from the endpoint state. Even if the source readout retains amplitude information, assigning one join to one tick and `\ln d` to a metric unit are extra interpretive steps. Moreover, a set of graph paths is not a coherent Feynman sum until relative phases and a composition law are supplied.

## 2.1. Noncommuting Bell projections: a candidate event-conflict relation

**[ADDED]** Represent three qudits on the *ordinary* tensor product `H₁⊗H₂⊗H₃` and put

\[
A=P_{12}\otimes I_3,\qquad B=I_1\otimes P_{23},
\]

where each `P` projects onto an equal-dimensional maximally entangled pair. These two projections overlap on the middle factor. Their composition satisfies the exact identity

\[
\boxed{ABA=\frac1{d^2}A}.
\]

**Proof.** Inserting the explicit Bell-pair expansions, the two contractions through the common index produce one factor `1/d` each; the remaining outer projection is `A`. This is an operator equality on `H₁⊗H₂⊗H₃`, not just an equality of matrix elements in a chosen state. □

Since both `A` and `B` are Hermitian projections, the Hilbert–Schmidt norm of their commutator is

\[
\begin{aligned}
\|[A,B]\|_{\rm HS}^2
&=\operatorname{tr}[(AB-BA)^\dagger(AB-BA)]\\
&=2\operatorname{tr}(AB)-2\operatorname{tr}(ABAB).
\end{aligned}
\]

The sandwich identity and `rank A=d` give

\[
\operatorname{tr}(AB)=\operatorname{tr}(ABA)=\frac{\operatorname{tr}A}{d^2}=\frac1d,
\qquad
\operatorname{tr}(ABAB)=\frac1{d^3}.
\]

Therefore

\[
\boxed{\|[A,B]\|_{\rm HS}^2=\frac{2(d^2-1)}{d^3}>0,\quad d\ge2.}
\]

**[TESTED]** At `d=2,3,4`, the respective squared commutator norms are `3/4`, `16/27`, `15/32` and the `ABA=A/d²` residuals are below `2×10⁻¹⁶`.

**Interpretation and essential caveat.** A pair of overlapping Bell projections is order-sensitive, so this offers an *event-conflict* relation to study. It does not produce an oriented arrow of time: `A` preceding `B` and `B` preceding `A` are both algebraically definable. More importantly, the user's source treats incidents on a shared named qudit as distinct virtual copies, whereas this lemma treats the same middle tensor factor as physically shared. Equating them would be an **additional operational axiom** and must not be silently assumed. The lemma is useful as a design test for a causal completion, not as a theorem that Chain Theory already predicts time.

## 3. Local link type does not fix global topology

### Proposition 4: Two graphs, same local structure, different cycle spectrum

Define two bipartite graphs `G_A,G_B` with qutrit nodes `Q_i`, qubit nodes `q_i`, and indices `i∈ℤ₇`. Give each edge type `(3,2)` and the maximally entangled rank-two pair state on that edge.

- `G_A`: `Q_i` attaches to `q_{i+s}` for `s∈{0,1,2}`.
- `G_B`: `Q_i` attaches to `q_{i+s}` for `s∈{0,1,3}`.

Both graphs have 14 nodes, 21 edges, are connected, are bipartite, and have degree 3 at **every** node. In particular, every qutrit sees a three-link star of exactly the source's locally selected `(3,2)` type; two of its three links may be marked as reads when studying the source's two-point object.

**[TESTED]** Their invariants differ:

| Invariant | `G_A` | `G_B` |
|---|---:|---:|
| Vertices | 14 | 14 |
| Edges | 21 | 21 |
| Degree | 3 | 3 |
| Four-cycles | 7 | 0 |
| Girth | 4 | 6 |
| `tr(A⁴)` | 266 | 210 |

**Algebraic explanation.** Two qutrit vertices support a 4-cycle exactly when they have at least two common qubit neighbors. For a circulant attachment set `S⊂ℤ₇`, common-neighbor multiplicities are governed by differences of elements of `S`. In `S_B={0,1,3}`, the six ordered nonzero differences are all distinct modulo 7; hence no two distinct qutrit vertices share two neighbors and there are no 4-cycles. For `S_A={0,1,2}`, one difference repeats, producing seven 4-cycles. Each four-cycle adds eight oriented closed length-four walks, accounting for `266−210=8×7`.

**Scope.** These are **combinatorial** counterexamples for networks with the same local allowed link template. They respect the source's per-link copy language as a formal pair network. This does not establish that both graphs satisfy every possible future emergent spacetime law, nor does it resolve the physical identification of several virtual copies as one qudit. It establishes that such an additional selection law has not been provided by the existing local state axiom.

## 4. Even explicit causality and symmetry do not choose a unique walk

### Proposition 5: A continuously adjustable family of equally causal local updates

**[ADDED]** Take a spatial ring with `N=25` sites, each carrying a two-component amplitude `ψ_x∈ℂ²`, and a discrete clock. Let

\[
C_m=e^{-im\sigma_x},\qquad
(S\psi)_{x,+}=\psi_{x-1,+},\quad(S\psi)_{x,-}=\psi_{x+1,-},
\]

and `U_m = S C_m`. The unitary acts exactly one lattice edge per step after an onsite mixing operation; its Fourier symbol is

\[
\boxed{U_m(k)=e^{-ik\sigma_z}e^{-im\sigma_x}}.
\]

The following are exact for **every real** `m`:

1. `U_m†U_m=I`, because both shift and onsite coin are unitary.
2. Finite-speed causality: an initial state at position zero has zero support beyond graph distance `t` after `t` steps (before wraparound).
3. Translation invariance: `U_m` commutes with the one-site translation `T`.
4. Parity invariance: `U_m` commutes with `P`, where `(Pψ)_{x,+}=ψ_{-x,-}` and vice versa.

These do **not** select `m`. One-step probabilities for initial `|x=0,+⟩` are `Pr(x=+1)=cos²m`, `Pr(x=−1)=sin²m`:

| `m` | `Pr(+1)` | `Pr(−1)` |
|---:|---:|---:|
| 0 | 1.000000 | 0 |
| 0.21 | 0.956544 | 0.043456 |
| 0.61 | 0.671823 | 0.328177 |

**[TESTED]** Full `50×50` matrices give unitarity, translation, and parity residuals all below `3×10⁻¹⁷`; at steps `1` through `6` the probability strictly outside the nearest-neighbor causal cone is zero to machine precision. The three dynamics are demonstrably distinct.

Near `k=m=0`, this update recovers the free `1+1` Dirac dispersion, as the v3 experiment did for the reversed ordering `C_m S`. But parity-symmetric local evolution still leaves at least the coin parameter undetermined. In 3+1 dimensions there are additional choices of graph, spinor content, and lattice species.

The literature has achieved stronger conditional uniqueness results by **adding** unitarity, homogeneity, locality, discrete isotropy, minimal complexity, and a specified class of graphs: D'Ariano and Perinotti, *Physical Review A* **90**, 062106 (2014); Raynal, *Physical Review A* **95**, 062344 (2017). Those theorems do not follow solely from this Chain Theory source's pair-and-read axiom.

## 5. A small renormalization obstruction

For a ring of `N` even sites, a disjoint nearest-neighbor pairing can be

\[
\mathcal B_0=\{(0,1),(2,3),\ldots,(N-2,N-1)\},
\]

or instead

\[
\mathcal B_1=\{(1,2),(3,4),\ldots,(N-1,0)\}.
\]

Translation by one site exchanges these pairings. Therefore no individual choice among these two pairwise blockings is invariant under the ring's full one-site translation symmetry, although their orbit as a set is invariant. **[TESTED]** For `N=6,8,10`, direct set enumeration confirms this exactly.

This is not a no-go theorem for arbitrary real-space renormalization: one can average over blockings, choose a different covariant blocking, or impose a staggered reference frame. The point is that a **canonical** RG operation is not supplied merely by the existence of finite pairs. RG is indispensable for a beta function.

## 6. The postselection caveat: what is and is not ruled out

The original source already acknowledges that its no-failure-branch postselection is not the complete ordinary laboratory protocol, and notes that fixed normalization holds on its specified restricted configurations. We should not falsely turn an ordinary postselection objection into a contradiction of a claim the paper never made.

**[ADDED warning test].** If one promotes a normalized Bell projection to an **unconditionally deterministic map on arbitrary quantum inputs**, it permits signaling. Take `AB=|ω₂⟩` and independently prepare `C=|c⟩` for `c=0` or `1`. Then

\[
{}_{BC}\langle\omega_2|\left(|\omega_2\rangle_{AB}|c\rangle_C\right)=\frac12|c\rangle_A.
\]

After branch normalization, the state of remote `A` becomes `|c⟩`, whereas its initial marginal was `I/2`. Normal Bell postselection cannot be used for signaling because it succeeds only probabilistically and identifying that branch requires classical communication. If the same output could be imposed deterministically for arbitrary preparations, the remote marginal would change under the sender's choice.

**[TESTED]** The tensor contraction and remote output are exact numerically for both `c` values. **Limit:** a freely prepared unpaired `C` lies outside the source's restricted 'nothing but maximally entangled pairs' configurations. Thus this is a *consistency condition on extensions with free matter/initial-state preparation*, not a proof of superluminal signaling inside the original toy ontology.

For a future physical theory admitting charged excitations and arbitrary local preparation, it will be necessary to specify an operationally sound causal update, perhaps using a complete measurement-and-correction channel instead of literal deterministic postselection.

## 7. State of the QED coefficient and the next exact gate

| Required ingredient | Supplied by original qudit paper? | v4 result |
|---|---|---|
| Finite Hilbert spaces and pair-state algebra | Yes | Checked source projection normalization |
| Acting hypercharge-like generator | Yes | Retained `(3,2)₁/₆` link |
| Graph incidence and local link type | Yes | Different global graphs remain possible |
| Unique physical spacetime geometry | No | Explicit global counterexample for local templates |
| Physical causal evolution/clock | No | Many admissible *added* local walks |
| Canonical RG coarse-graining scale | No | Pair blocking has a symmetry ambiguity |
| Fermionic chirality and unique light spectrum | No | v3 counterexamples remain |
| Momentum-dependent electromagnetic polarization | No | No original-axiom formula |
| `β(e)=e³/(12π²)` independently | **No** | **Not derived** |

### Most valuable next problem: formalize a causal operational completion

A constructive program would:

1. **Specify the state category and composition:** define the tensor factors behind 'shared copies are one qudit', the allowed pair contractions, and the status of postselected complementary branches. Prove closure and probability consistency on the admissible configurations.
2. **Build a derived adjacency criterion:** characterize which pairs of read events can affect one another from the axiom, without importing a target lattice. Test whether the two competing 14-vertex circulant graphs are excluded by any source-derived invariant.
3. **Extract a clock or causal partial order:** investigate whether noncommuting read operations canonically order events. The current source explicitly gives its two read insertions no time order. The overlap-projector commutator in §2.1 is a possible conflict test but its physical tensor-factor meaning and orientation remain additional choices. A sequential QCA requires this new structure.
4. **Derive a natural unitary/circuit generator, if possible:** show uniqueness under the exact source matching group, not solely under familiar external symmetries. A minimality or locality assumption, if added, must be openly named and justified.
5. **Only afterward build charge dynamics and renormalization:** derive a conserved current and Ward identity at finite momentum, and obtain the logarithmic polarization slope without inserting the Dirac determinant or the QED beta coefficient by hand.

**Go/no-go test:** If two physically inequivalent propagating, gauge-covariant, causal completions remain compatible with all original read and exchange constraints, the original axiom does not select one beta function. A revised axiom would then be required; further precision on the zero-momentum record count would not resolve that underdetermination.

## 8. Reproducibility and limits

`chain_beta_extension_v4.py` tests ordinary state-vector and tensor contractions, graph invariants, translation/parity/unitarity of finite matrices, finite causal cones, and elementary RG partition combinatorics. It depends on NumPy and NetworkX, uses no measured fine-structure constant, and does not require an external field-theory package.

These are finite algebraic checks, **not** a simulation of vacuum-polarization loops or of an interacting 3+1D QFT continuum limit. All numerical identities have closed-form proofs above. Floating-point residuals cannot establish an original-axiom beta-function theorem, and the v4 results should not be presented as such.

### Cited literature

- A. D. Wissner-Gross, *A Parameter-Free Candidate for the Fine-Structure Constant from Entangled Qudits*, 23 October 2026, user-provided 24-page manuscript, main paper and Supplement §§S1, S3, S7, S8.
- G. M. D'Ariano and P. Perinotti, *Derivation of the Dirac equation from principles of information processing*, **Phys. Rev. A 90**, 062106 (2014), DOI: 10.1103/PhysRevA.90.062106. <https://journals.aps.org/pra/abstract/10.1103/PhysRevA.90.062106>
- P. Raynal, *Simple derivation of the Weyl and Dirac quantum cellular automata*, **Phys. Rev. A 95**, 062344 (2017), DOI: 10.1103/PhysRevA.95.062344. <https://journals.aps.org/pra/abstract/10.1103/PhysRevA.95.062344>
- P. Arrighi, V. Nesme, and R. Werner, *Unitarity plus causality implies localizability*, **J. Comput. Syst. Sci. 77** (2011), 372–378, DOI: 10.1016/j.jcss.2010.05.004. <https://www.sciencedirect.com/science/article/pii/S0022000010000589>

### Files

- `chain_beta_extension_v4.py` — tests and closed-form example definitions.
- `chain_beta_extension_v4_results.txt` — captured reproducible output.
- `ChainTheory_QED_beta_extension_v4.md` — this technical note.