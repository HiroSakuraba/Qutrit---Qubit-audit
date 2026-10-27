# Chain Theory and the QED beta function — v5 operational-completion audit

**Date:** 24 October 2026  
**Status:** Exact conditional consistency result; **not** a derivation of the QED beta function.  
**Primary source:** A. D. Wissner-Gross, *A Parameter-Free Candidate for the Fine-Structure Constant from Entangled Qudits* (23 October 2026), supplied `ChainTheoryAlpha.pdf`, main text and Supplemental §§S1, S3, S7, S8.  
**Reproducibility:** `chain_beta_extension_v5.py`, `chain_beta_extension_v5_results.json`.

## Executive conclusion

The previous work (v2–v4) established that ordinary local qudit lattice models can reproduce conditional QED calculations after importing dynamical assumptions, but the Chain Theory axioms do not select an emergent geometry, clock, propagating fermion, or running charge. This update audits the **operational status of the source's no-failure-branch distillation** before attempting to obtain any such dynamics.

The strongest new result is a conditional obstruction with **only source-admissible pair configurations**: two maximally entangled `(3,2)` links have a maximally mixed joint remote-qubit state, while a normalized local qutrit–qutrit Bell projection makes their remote qubits a pure Bell pair. Those two states have trace distance **3/4**. If (and only if) Alice may freely select this operation while Bob has independently accessible measurements on his remote two-qubit system, deterministic selected-branch dynamics violates the standard no-signalling condition. Crucially, the original source **does not** assert those added spacetime, free-choice, and measurement conditions. Thus this is not a contradiction internal to its ontology; it is an explicit constraint on an operational extension.

A second independent calculation contracts up to five equal-dimensional Bell links in **every possible order of their permitted disjoint virtual-copy joins**. The final unnormalized state is exactly the same in all orders. The maximally entangled pair is invariant under all `U ⊗ conjugate(U)` matchings. These facts show that the pair calculus supplies **composable contractions**, but not an intrinsic time orientation.

## 1. Facts taken directly from the source

1. A type `(d_1,d_2)` link is a pure maximally entangled pair on common support `m=min(d_1,d_2)`; a shared qudit has one *virtual copy per incident link* (main text pp. 1–2, Supplement S1).
2. Distillation is the projection `P_m = |ω_m⟩⟨ω_m|` on a pair. The model retains the projected configuration and explicitly excludes the complementary outcome as a state of its ontology (main text pp. 1–2; Supplement S1, especially postselection discussion).
3. Its two qutrit–qubit links contracted through the qutrit ends have **ordinary** selected-outcome norm squared `1/6`; the model renormalizes the selected vector by `√6` (Supplement S7, its joint state, and Eq. S2).
4. The pair is invariant under its matching group `U ⊗ conjugate(U)` (main text p. 2, Supplement Lemma S1).
5. The author says that the two insertions have **no time order**, and that the construction supplies no momentum-dependent propagator or physical momentum parameter (Supplement S7–S8).

These facts alone must not be augmented with classical causal structure, an arbitrary mixed-state ontology, accessible complementary measurement outcomes, or freely choosable interventions without marking them as additional assumptions.

## 2. Exact calculation on original pair-type inputs

Let `A_1` and `A_2` each be a **qutrit virtual copy**, sharing a common site under the source's language, and let `B_1,B_2` each be a qubit. Define two admissible links:

\[
|\Psi\rangle=|\omega_2\rangle_{A_1 B_1}|\omega_2\rangle_{A_2 B_2},\qquad
|\omega_2\rangle_{AB}=\frac{|00\rangle+|11\rangle}{\sqrt2}.
\]

The qutrit has basis states `0,1,2` and each `(3,2)` link occupies common support `0,1`. Before any projection, Bob's two-qubit marginal is

\[
\rho_B=\operatorname{tr}_{A_1A_2}|\Psi\rangle\langle\Psi|=I_4/4.
\]

Apply the source's qutrit-end Bell bra,

\[
\langle \omega_3|_{A_1A_2}=\frac{\langle00|+\langle11|+\langle22|}{\sqrt3}.
\]

The unnormalized remote qubit state is, by direct contraction,

\[
\begin{aligned}
\langle\omega_3|_{A_1A_2}\Psi\rangle
&=\frac1{2\sqrt3}(|00\rangle+|11\rangle)_{B_1B_2}\\
&=\frac1{\sqrt6}|\omega_2\rangle_{B_1B_2}.
\end{aligned}
\]

The ordinary Born success probability is `1/6`. After conditioning and normalization, the entire post-projection configuration comprises **two maximally entangled pairs**:

\[
|\omega_3\rangle_{A_1A_2}\otimes|\omega_2\rangle_{B_1B_2}.
\]

Thus Bob's **conditional** state is `P_B=|ω_2⟩⟨ω_2|`, while the initial marginal was `I_4/4`. Their eigenvalue difference is `(+3/4,-1/4,-1/4,-1/4)`, giving

\[
\boxed{D_{\rm tr}(P_B,I_4/4)=\frac12(3/4+3/4)=3/4.}
\]

One explicit distinguishing test is Bob's Bell projector. Its probability is `1` if the selected branch is known to occur and `1/4` with no swap. Under standard equal-prior state discrimination, these density matrices have the Helstrom optimal success probability `(1 + 3/4)/2=7/8`; this is **not** a measurement protocol supplied by Chain Theory.

### A complete instrument restores no-signalling

Compare the source's always-selected projection with a standard **nonselective** two-outcome measurement on Alice,

\[
\mathcal E_A(\rho)=P_A\rho P_A+(I-P_A)\rho(I-P_A),
\quad P_A=|\omega_3\rangle\langle\omega_3|.
\]

Tracing out Alice yields precisely the original Bob marginal, `I_4/4`. The unnormalized joint remote states of the two branches are

\[
\sigma_\mathrm{yes}=\frac16 P_B,
\qquad
\sigma_\mathrm{no}=\frac14 I_4-\frac16 P_B.
\]

They have traces `1/6` and `5/6`. Therefore the normalized failed-branch remote state is

\[
\rho_{B|\mathrm{no}}=\frac{3}{10}I_4-\frac15 P_B.
\]

Its Bell-projector probability is `1/10`. Averaging over outcomes yields

\[
\Pr(P_B)=\frac16\cdot1+\frac56\cdot\frac1{10}=\frac14,
\]

exactly as required by no-signalling.

**Conditional consistency proposition.** If one simultaneously stipulates **(a)** independently accessible Alice/Bob subsystems, **(b)** Alice's freedom to choose between an inert operation and the normalized qutrit join, **(c)** guaranteed success of that join with no heralding requirement, and **(d)** Bob's ability to measure the two-qubit Bell projector, then these operations cannot satisfy standard no-signalling. Proof: Bob's predicted statistics differ by `3/4` depending on Alice's choice. Conditions (a), (b), and (d), especially their spacetime and operational content, **are additions**; the source does not establish them. No claim that physical superluminal signalling is already demonstrated follows.

This result is stronger and more source-faithful than the v4 exploratory example, which used an independently prepared unpaired state and therefore departed from the original axiom.

## 3. Join-order test on virtual copies

For `n` Bell links of dimension `d`, written as tensors with `2n` separately labelled ends, join each adjacent internal endpoint pair with the Bell bra. Each joining bra acts on a disjoint pair of **virtual tensor factors**; hence these contractions commute as linear tensor operations, regardless of which is evaluated first. Their final unnormalized output is

\[
\boxed{
\bigl(\langle\omega_d|\bigr)^{\otimes(n-1)}
\bigl(|\omega_d\rangle\bigr)^{\otimes n}
=d^{-(n-1)}|\omega_d\rangle_{\rm endpoints}.
}
\]

Proof by matrix multiplication: each Bell tensor is `I_d/√d`. Contraction across two neighbouring links contributes `I_d/(d√d)` on the endpoints, which is `1/d` times the normalized Bell tensor. Adding one more link and one joining bra multiplies the coefficient by another `1/d`. Tensor contractions over disjoint labels commute, so evaluation order is immaterial.

The test suite explicitly enumerates every join ordering for `d=2,3,4` and `n` up to 5. The largest case has `4!=24` orderings. All agree, with maximum Frobenius errors below `8×10^-17`.

**Scope:** This statement is about contraction ordering for independently labelled virtual ends, not arbitrary overlapping projections acting on the same physical factor (the latter need not commute, as v4 already showed). The source itself does not define a temporal orientation from either case. An event poset, a free intervention, and spacelike commutativity remain unassigned.

## 4. Matching invariance is not a propagator

The normalized Bell state obeys

\[
(U\otimes\overline U)|\omega_d\rangle=|\omega_d\rangle
\quad\text{for every }U\in U(d).
\]

This follows from `\sum_i U_{ai}\overline U_{bi}=δ_{ab}`. Tests with random complex Haar-type unitary matrices in dimensions 2, 3, and 4 give errors smaller than `10^-15`.

These matchings are changes of basis under which the pair state is unchanged. They cannot be treated as nontrivial propagation operators just because they form a unitary group: their entire orbit on `|ω_d⟩` is a single state.

There are also many possible **added** Hamiltonians that preserve this state. On a Bell pair, for example, both

\[
H_0=0,\qquad H_1=I-P_d
\]

leave `|ω_d⟩` stationary (the second has zero eigenvalue on the Bell state). Their evolutions `e^{-itH_0}` and `e^{-itH_1}` disagree off that state, but off-state disturbances and a time coordinate are not source-generated. The experiment illustrates underdetermination, not new dynamics inside the source ontology.

## 5. Why this blocks a QED beta derivation

The original value `α_cand(0)` is a finite normalized record sum. It carries no scale. A nonzero beta function requires at least a class of momentum-dependent charged excitations, a causal/operational current two-point function, and a scale-changing transformation. This v5 audit reveals an additional gate: **a consistent intervention calculus** specifying what 'an act' means, especially for the no-failure branch of a local join.

Reproducing the familiar one-loop factor `2/(3π)` by appending a Dirac determinant remains possible, as v2 demonstrated, but is not a derivation from these axioms. The qudit representation-theoretic `5/3` charge-squared trace and the eventual beta function should not be conflated.

### Required next decision — explicit alternatives, not assumed equivalence

- **A. Acausal/global-constraint interpretation.** Treat distillations as correlations or boundary conditions, not freely controlled interventions. This may avoid operational signalling; it does not yet produce a clock, locality, or time-directed propagation.
- **B. Causal instrument completion.** Introduce full trace-preserving measurement outcomes and heralded classical communication/corrections. Standard no-signalling is restored. This **modifies the literal no-failure branch** and may affect the interpretation of the record sum; its preservation must be checked rather than assumed.
- **C. Restricted intervention theory.** Permit selected-branch operations only under a new rule that guarantees no operational signalling even where separate subsystems are meaningful. This requires a precise new selection principle and should be treated as a second axiom unless derived from the first.

A rigorous future uniqueness theorem must state a class of allowed configurations, a set of events/instruments, a locality relation, and a composition law. It must then show there is **one** equivalence class of causal propagating charged theories in that class. At present neither the existence nor uniqueness of such a completion has been proved.

## 6. Verification and limitations

`chain_beta_extension_v5.py` uses only NumPy and SymPy; it checks exact rational two-qubit density-matrix identities, enumerates contraction orders, verifies matching invariance, and constructs inequivalent off-subspace unitary completions. `chain_beta_extension_v5_results.json` is the captured machine-readable output. The earlier v2–v4 files remain unmodified.

No simulations of physical propagation, four-dimensional renormalization, or charged particle loops are reported here. All references to 'signalling' are explicitly conditional on standard operational assumptions absent from the source. The Chain Theory formula for `α(0)` is neither confirmed nor falsified by this test.

## References

1. Wissner-Gross, *A Parameter-Free Candidate for the Fine-Structure Constant from Entangled Qudits* (23 Oct 2026), `ChainTheoryAlpha.pdf`, main and Supplement §§S1, S7, S8 (user-provided primary source).
2. Żukowski, Zeilinger, Horne, and Ekert, *Event-Ready-Detectors Bell Experiment via Entanglement Swapping*, Phys. Rev. Lett. 71, 4287 (1993).
3. Aaronson, *Quantum computing, postselection, and probabilistic polynomial-time*, Proc. R. Soc. A 461, 3473 (2005), https://doi.org/10.1098/rspa.2005.1546.
4. IETF RFC 9340, *Architectural Principles for a Quantum Internet*, §4.4.2, https://www.rfc-editor.org/rfc/rfc9340.html (standard Bell-measurement outcomes and the necessity of communicating outcome labels).