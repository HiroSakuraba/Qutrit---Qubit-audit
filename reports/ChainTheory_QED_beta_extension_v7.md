# Chain Theory and the QED beta function — v7
## Fixed qutrit-join impossibility, covariant support transport, and the residual phase ambiguity

**Date:** 24 October 2026  
**Research status:** Exact finite-dimensional no-go statement with an explicit covariant operational *completion*. No independent QED beta-function derivation.  
**Primary text:** A. D. Wissner-Gross, *A Parameter-Free Candidate for the Fine-Structure Constant from Entangled Qudits* (23 October 2026), `ChainTheoryAlpha.pdf`, especially pp. 1–5 and Supplemental §§S1, S3, S7–S8.  
**Reproduction:** `chain_beta_extension_v7.py` and `chain_beta_extension_v7_results.json` (Python 3, NumPy, SymPy).

## 0. Executive findings

1. **General rank-one fixed-join obstruction.** On the class of **independently oriented** (3,2) maximally entangled pair states, *no nonzero fixed rank-one bra on the two qutrit factors* maps every nonzero selected outcome to a maximally entangled remote qubit pair. This holds for **every** choice of fixed bra, not only the paper's Bell projector. The proof is by elementary rank-nullity and is given below. Importantly, this is a **conditional theorem**: the source does not specify the physical admissibility of independent support rotations versus joint relabelling of the entire join.
2. **A literal finite-twirl identity.** Twenty-four explicit unitary basis changes (cyclic shifts and diagonal signs) on **one** of the two qutrits map the rank-one Bell projector under averaging to the full-rank operator `I9/9`. This is an exact identity. It proves that discarding an unspecified *independent* relative frame loses a rank-one Bell selection. It does **not** apply to the correlated `U ⊗ conjugate(U)` symmetry which fixes the original Bell state.
3. **A successful covariant construction exists.** If one supplies oriented embeddings `U1,U2 in U(3)` extending the two link supports and transforms the join with them, the join has success probability `1/6` and produces a maximally entangled remote qubit pair for **all** such supports. Transporting the v6 nine-Kraus, measurement-and-repreparation instrument also preserves no-signalling and gives exclusively pair states on each realized branch of these transformed source inputs. This is a family of input-frame-adapted instruments, **not** a unique fixed physical instrument selected by the original axiom.
4. **The SU(3) volume form partly resolves the missing third direction.** Given a **phased** 3x2 isometry `V`, the normalized conjugate cross product of its two columns gives a canonical third column, yielding an SU(3) matrix. It is exactly covariant under SU(3). However, the input *physical ray* does not select the phase of `V`. Under `V -> exp(i phi) V` the two occupied columns gain `exp(i phi)` while the unused complement gains `exp(-2i phi)`. At `phi=pi/3` the corresponding full-qutrit join projectors have trace distance `2 sqrt(2)/3` despite being generated from identical input pair rays. The **original selected remote amplitude is still identical**. Whether these differing projectors are gauge-equivalent in the full theory requires an operational definition of admissible join and outcome observables.
5. **The projective stabilizer gives a second impossibility result.** A join defined solely from the two physical link rays cannot be simultaneously `SU(3)`-equivariant, rank one, and maximally entangled on the full `3x3` space: the stabilizer of each two-dimensional occupied plane forbids rank-three Schmidt support. This result holds under the explicit assumption that the independent stabilizer is gauge.
6. **The 172-record value has not been disproved or established as an unconditional coupling.** If all record projectors, exchange maps and readout operators are simultaneously conjugated, their spectra and normalized traces are trivially invariant. The original source does not define that full frame transport as a physical procedure, nor does it define readouts for other instrument outcomes. Therefore the physical running of alpha remains underdetermined.

## 1. Source axioms and what they do *not* settle

The source starts with maximally entangled pair states on common support, a shared qudit represented by one virtual copy per link, reads implemented through projections, and a selected join of two qutrit ends. In its configuration, two links each have local support dimension two inside a three-dimensional qutrit, and the qutrit join is `P3 = |omega3><omega3|`. The source excludes branches outside its pair-only ontology and states that the relative basis matchings of different pairs are not lined up (main text p. 2, paragraph on matchings). The source states in Supplemental S8 that it constructs no momentum-dependent propagator.

Two distinct mathematical uses of the word *matching* must not be conflated:

- **Correlated Bell-state matching:** `(U ⊗ conjugate(U)) |omega3> = |omega3>`. This preserves the Bell projector exactly.
- **Independent input-support orientation:** act on the qutrit end of one `(3,2)` input by `G1` and on the other by an independent `G2`, but hold the join bra fixed. The original projector generally does not preserve Bell-pair outputs in this experiment.

The source has not given a complete intervention calculus resolving which description defines physically distinct operations. The theorems below explicitly condition on the latter interpretation rather than claiming a contradiction with the former.

## 2. General fixed rank-one join no-go

Let `V1,V2` be complex `3x2` isometries, so `Vi† Vi = I2`. Define the input link states

\[
|\psi(V_i)\rangle=\frac1{\sqrt2}\sum_{r=0}^1V_i|r\rangle_{A_i}\otimes |r\rangle_{B_i}.
\]

Every such state is maximally entangled between its two-dimensional support in `A_i` and the qubit `B_i`.

Let a **fixed** nonzero rank-one selected outcome be `P_B = |J_B><J_B|`, where `|J_B> = sum_{a,b=0}^2 B_{ab}|a,b>` for a fixed nonzero `3x3` matrix `B`. Applying its bra on `A1 A2` leaves an unnormalized pure two-qubit state whose `2x2` coefficient matrix is

\[
\boxed{C(V_1,V_2;B)=\tfrac12 V_1^T\overline{B}\,V_2.}
\]

The selected probability is `p = tr(C† C)`. Conditional on `p>0`, the remote pair is maximally entangled if and only if `C C† = (p/2) I2`. In particular it must have rank two.

**Theorem.** For every nonzero `B`, there exist isometries `V1,V2` such that `p>0` and `rank(C)=1`.

**Proof.** Since `B ≠ 0`, choose a two-plane `V1` for which `M=V1^T conjugate(B)` is nonzero. Because `M` has shape `2x3`, there exists a unit vector `w0 in ker M`. Since `M ≠ 0`, choose a unit vector `w1 orthogonal to w0` with `Mw1 ≠ 0` (otherwise `M` would vanish on the whole orthogonal complement and also on `w0`, contradicting `M≠0`). Let `V2=[w0,w1]`, a `3x2` isometry. Then

\[
C=\frac12[Mw_0,Mw_1]=\frac12[0,Mw_1]
\]

is nonzero but has rank exactly one. QED.

**Numerical stress:** 128 seeded complex constructions at each of ranks 1, 2 and 3 for `B`, with the kernel constructed using singular-value decomposition. Across 384 constructed counterexamples, the maximum computed remote concurrence was `1.03e-15`; all branches had strictly positive probability. This is a check of the constructive proof, not a statistical proof by random sampling.

**Exact instance for the paper's `P3`:** Choose first support `span{e0,e2}`, second `span{e0,e1}`. Then `C = diag(1,0)/(2 sqrt(3))`; `p=1/12` and the selected remote state is `|00>`, not a Bell pair.

## 3. Exact 24-element matching twirl

Let `X` cyclically permute `e0→e1→e2→e0`, and let `D_s = diag(s0,s1,s2)` with `s_i ∈ {+1,-1}`. Consider 24 unitaries

\[
\mathcal G=\{X^k D_s:\ k=0,1,2,\;s_i=\pm1\}.
\]

The 24-element average on the *first* qutrit of `P3` is

\[
\boxed{\frac1{24}\sum_{U\in\mathcal G}(U\otimes I_3)P_3(U^\dagger\otimes I_3)=\frac{I_9}{9}.}
\]

Proof: the eight sign choices dephase the first-qutrit matrix index; the three cyclic shifts make the three diagonal blocks equal. In formula, the induced map on any `3x3` matrix is `A -> (tr A) I3/3`. Since `tr_{A1} P3 = I3/3`, the output is `I3/3 ⊗ I3/3`.

This uses **independent** rotation of one qutrit relative to the join. Under the correlated `U ⊗ conjugate(U)` action `P3` is fixed, so it would be incorrect to use this identity as evidence that the *correlated isotropic twirl of the paper* destroys the selected state. See [QETLAB's definitions of isotropic and Werner twirls](https://qetlab.com/Twirl) for the distinction.

For the original aligned two-link input, the effective averaged POVM element `I9/9` gives success probability `1/9`; its conditional Bob state is `I4/4`. These values follow exactly from the full-rank averaged effect and the original maximally mixed Bob marginal. They are **not** the source's `p=1/6` selected process and are not presented as the paper's prescription.

## 4. Covariant support-following rank-three join

Choose unitary completions `U1,U2 ∈ U(3)` of the isometries: `Vi=Ui[:,0:2]`. Define

\[
|J(U_1,U_2)\rangle=(U_1\otimes U_2)|\omega_3\rangle,\qquad
P(U_1,U_2)=|J(U_1,U_2)\rangle\langle J(U_1,U_2)|.
\]

Then for the input `|psi(V1)> ⊗ |psi(V2)>`, direct contraction gives

\[
(\langle J(U_1,U_2)|\otimes I_B)|\psi(V_1)\psi(V_2)\rangle
=\frac1{\sqrt6}|\omega_2\rangle_B.
\]

Hence `p=1/6`, with remote concurrence exactly one. Under independent physical transformations `G1,G2`, if the join is co-transformed as `P -> (G1⊗G2)P(G1⊗G2)†`, the selected amplitude is unchanged. This is genuine **covariance** of a family of joins, not invariance of one fixed projector.

The 128 Haar-unitary numerical trials gave:

- maximum selected-branch probability/concurrence error `7.8e-16`;
- maximum transported-nine-Kraus outcome probability error `5.0e-16`;
- maximum remote-state no-signalling residual `1.14e-15`;
- maximum remote Bell-concurrence error across nonzero branches `8.9e-16`.

The transported CPTP instrument is constructed from the v6 nine-dimensional orthonormal input basis `{v_j}` by setting `T=U1⊗U2`, `|v'_j>=T|v_j>`, and

\[
K'_j=|J(U_1,U_2)\rangle\langle v'_j|=T K_j T^\dagger.
\]

Every transformed instrument is CPTP: `sum_j K'_j†K'_j = I9`. The selected Kraus operator is precisely its transformed Bell projector `K'_0=P(U1,U2)`. For the correspondingly transformed **source input**, five branches occur with probabilities `(1/6,1/12,1/4,1/4,1/4)`, and every realized branch consists of a qutrit Bell pair and a qubit Bell pair. The uncommunicated remote state remains `I4/4`.

**Important qualifications:** (i) to choose the correct `T` for each physical input requires that the relative frames or a connection be available as operational data; (ii) the nonselected Kraus operators perform replacement/repreparation, explicitly excluded by the original source's Corollary S2; (iii) the theorem is about the family of transformed source inputs, **not** arbitrary entangled many-body configurations. Simply promoting this family to the original theory would be an additional dynamical/operational assumption.

## 5. Does the SU(3) epsilon tensor determine the missing third column?

This is a nontrivial improvement over simply asserting the need for an arbitrary third-dimensional frame. Given the **phased** isometry `V=[v0,v1]`, define

\[
w(V)=\overline{v_0\times v_1},\qquad
\widetilde U(V)=[v_0,v_1,w(V)].
\]

Here the cross product uses the ordinary complex-bilinear Levi-Civita tensor in a chosen `SU(3)`-oriented basis. Direct algebra shows `Ũ(V)† Ũ(V)=I3` and `det Ũ(V)=1`. Moreover for `G ∈ SU(3)`, `w(GV)=G w(V)`, so `Ũ(GV)=G Ũ(V)`.

This means the **SU(3) volume form itself provides a natural complement** to any specified phased two-column isometry. The numerical SU(3) equivariance residual across 128 random trials was at most `7e-16`.

However, the pair **ray** is unchanged by `V -> exp(i phi) V`. The canonical complement then changes as

\[
\widetilde U(e^{i\phi}V)
=\widetilde U(V)\,\mathrm{diag}(e^{i\phi},e^{i\phi},e^{-2i\phi}).
\]

At `phi=pi/3`, the first two columns gain `e^{i pi/3}` and the third gains `e^{-2i pi/3}=-e^{i pi/3}`. In the full rank-three Bell join this changes the sign of the `|22>` term relative to the occupied `|00>,|11>` terms. The two normalized join vectors have squared overlap `1/9`, yielding rank-one-projector trace distance

\[
\boxed{D(P_{\phi=0},P_{\phi=\pi/3})=\sqrt{1-1/9}=\frac{2\sqrt2}{3}.}
\]

The *physical input pair rays* are identical, and the *selected remote Bell amplitude* is unchanged. Distinct full join projectors can therefore arise from one ray-level input, although they may be gauge-related when all other measurement operators are transported. The paper's phase/basis freedom needs to determine whether that gauge equivalence is part of the formal operational definition.

The SU(3) epsilon construction is worth retaining as a candidate: it narrows the freedom substantially and may exploit the paper's identification group more economically than introducing an unconstrained U(3) connection. It does **not**, on its own, supply a well-defined rank-three join as a function of two **projective** input link states.


### 5A. A stronger stabilizer obstruction on the projective pair-state space

There is a direct way to express the phase issue as an equivariance theorem. Take the first input link in the standard plane `span{e0,e1}`. For every real `phi`, the SU(3) transformation

\[
G_\phi=\mathrm{diag}(e^{i\phi},e^{i\phi},e^{-2i\phi})
\]

acts on the input `(3,2)` pair ket simply as multiplication by `e^{i phi}`; hence it fixes the physical **ray**. The second input ray is unaffected. Therefore `(G_phi,I)` belongs to the stabilizer of this pair of *ray-level input states*.

Suppose a **rank-one, full-Schmidt-rank** qutrit–qutrit join projector `P(ray1,ray2)` were an equivariant function of those input rays alone. At the standard input it would have to satisfy

\[
(G_\phi\otimes I)P(G_\phi^\dagger\otimes I)=P
\qquad\text{for every }\phi.
\]

For a rank-one projector this means its spanning vector must be an eigenvector of `G_phi ⊗ I` for every `phi`. The first-qutrit eigenspaces are the two-dimensional plane with eigenvalue `e^{i phi}` and the one-dimensional complement with eigenvalue `e^{-2i phi}`. Because these eigenvalues are distinct for generic `phi`, any such eigenvector lies entirely in one of the two sectors. Its Schmidt rank across the two qutrits is at most two, contradicting the required rank three of a qutrit Bell state.

**Conclusion:** There is **no equivariant map from the two *projective input pair states alone* to a rank-one, maximally entangled (3,3) join projector**, under an action that treats the independent `SU(3)` stabilizer as a physical relabelling symmetry. Such a map may be defined on a larger space of *phased frames, reference identifications, or gauge connections*; that is exactly the extra structure that must be supplied or derived. This is a stronger mathematical statement than the fixed-bra theorem in §2, but its applicability to the source still depends on the actual gauge quotient and operational meaning of shared copies.

The special value `phi=pi/3` was checked symbolically in the v7 script: the source input changes only by a global phase, but the original full Bell projector changes by trace distance `2 sqrt(2)/3`.

## 6. Where the 172 records stand

The paper's fixed conditional readout remains

\[
N=24\cdot36\cdot36=31104,\quad
(n_1,n_{\sqrt2},n_{\sqrt3},n_{\sqrt6})=(79,48,42,3),
\]

\[
\alpha_{\rm cand}=(79+48\sqrt2+42\sqrt3+3\sqrt6)/31104.
\]

There are two different invariance questions:

- **Passive covariance of the complete record algebra:** If a full unitary map `W_R` transports *all* record operators, exchange actions and readout insertion as `M -> W_R M W_R†`, their traces and sector multiplicities remain unchanged. This is a tautological property of conjugation. The source asserts invariance of its counts under appropriate matchings, but does not define a complete operational channel with outcome averaging.
- **Actual physical readout after a no-signalling completion:** One needs a specified observable on every Kraus branch and a rule for summing uncommunicated outcomes. This is not supplied. The original `alpha_candidate` is exact **on its selected branch** and is not established as the unconditional photon coupling.

The v7 evidence does not give a new alpha value and makes no claim that the v6 conditional result is experimentally wrong.

## 7. Implications for the QED beta-function research program

The requested one-Dirac-species target remains

\[
\beta(e)=\frac{e^3}{12\pi^2}+O(e^5),\qquad
\beta(\alpha)=\frac{2}{3\pi}\alpha^2+O(\alpha^3).
\]

The original paper neither defines physical momentum nor supplies charged matter dynamics. The present v7 results add a sharper, earlier selection requirement: distinguish (a) symmetry under relabelling, (b) a relational join between independently specified supports, and (c) a *controllable* instrument with ordinary no-signalling.

A disciplined next research sequence is:

1. **Formalize the admissible symmetry group and objects:** Are individual link phases gauge, are qutrit support planes independent physical data, and what precisely acts on each virtual shared copy?
2. **Determine whether the SU(3) epsilon complement descends to the gauge quotient:** Test whether the source's central `U(1)` charge and its exchange definitions eliminate, compensate, or preserve the relative third-direction phase in every readout and join. A precise coequalizer/equivariance statement is preferable to another numerical value match.
3. **Specify an input-independent quantum instrument or prove none exists under the source rules:** A family selected using knowledge of each input is not an ordinary fixed channel unless that knowledge is obtained through allowed operations. Any new intervention must account for outcome records and nonselective dynamics.
4. **Only then introduce dynamics, scale, and fermion spectrum**, from which a beta function can be calculated. Do not insert one-loop factors by hand and call them a derivation from qudits.

## 8. Reproducibility and checks

Execute `python chain_beta_extension_v7.py > chain_beta_extension_v7_results.json` with NumPy and SymPy. The script checks exact SymPy results for the fixed-support counterexample, the 24-element twirl and complement-phase identity; it constructs 384 fixed-bra rank-one outcomes by nullspace methods; it tests 128 independent complex support orientations with transported Kraus instruments; and it checks SU(3)-volume-form equivariance on a further 128 random samples. Random seeds are fixed.

All implemented assertions pass. Machine floating-point residuals are about `1e-15` or smaller. These numerical checks corroborate the algebraic results. They do not prove the origin of spacetime, no-signalling from the source axioms, or the QED beta function.

## External background

- [QETLAB: Twirl](https://qetlab.com/Twirl) — distinguishes isotropic `U⊗conjugate(U)` twirling from Werner `U⊗U` twirling.
- Felix Leditzky, [Representation-theoretic methods in quantum information theory](https://www.felixleditzky.info/teaching/FT25/math595-repth-qit.pdf) — Haar averaging and representation-theoretic twirling. These are background sources; the no-go proofs above are independently written finite-dimensional algebra.