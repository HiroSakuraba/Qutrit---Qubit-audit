# Chain Theory → QED beta program, v13: a physical quantum instrument and the frame-invariance obstruction

**Date:** 24 October 2026  
**Primary source:** A. D. Wissner-Gross, *A Parameter-Free Candidate for the Fine-Structure Constant from Entangled Qudits* (`ChainTheoryAlpha.pdf`, 23 October 2026).  
**Related work:** v2–v12 in the companion archive.  
**Status:** Exact constructions and conditional obstruction proved; reproducible SymPy/NumPy checks passed. The first QED beta-function coefficient has **not** been derived from the source axioms.

## Abstract

We construct trace-preserving, completely positive (CPTP) instruments acting directly on the **physical** two-qutrit system in the paper's three-link, 216-dimensional quantum state. The selected Kraus operator is exactly the paper's qutrit–qutrit Bell projector `P3`; its probability is `1/6` and its normalized output is the three-pair joined state. An unconditioned local instrument provably preserves the remote reduced state. The aggregate **failure** state on the remote qubits is uniquely fixed, even though individual failure outcomes and their updates are highly nonunique. A nine-Kraus measure-and-reprepare completion realizes maximally entangled pairs on every nonzero outcome for the **aligned reference input**, but is forbidden by the paper's distillation-only rule and fails to preserve maximal entanglement for general independently oriented input supports.

A precise representation-theoretic conditional no-go follows. If *all independent* `SU(3) × SU(3)` changes of qutrit frames are treated as gauge equivalences of a **fixed** instrument and each finite outcome is invariant, its effects must be scalar multiples of the identity; the rank-one selected Bell effect cannot be one of them. A correlated `U ⊗ conj(U)` subgroup, by contrast, leaves that Bell projector invariant. A finite nine-element qutrit Weyl twirl verifies the difference exactly. The original source states matching invariance for the pair states; it does not establish that the entire independent product group must act as gauge on the already joined operation. Accordingly this is an explicit **conditional test**, not a falsification of the paper.

Finally, outcome probability `1/6` and the original 172-record scalar `alpha_cand≈0.00729735256366` are different quantities. A physical POVM can be calibrated to report the latter by **inserting that value** as a classical post-processing probability, but such calibration cannot count as a derivation. No physical construction here identifies the record-weight sum with the photon coupling or determines its running with energy.

## 1. Source assumptions and physical Hilbert space

The source's three links are qutrit–qubit maximally entangled pairs, with independent virtual copies `Q1,Q2,Q3` of a shared qutrit, remote qubits `q1,q2`, and unread qubit `f`. The paper chooses the qutrit Bell projector `P3=|omega3><omega3|` as its distillation join, normalizes the resulting postselected state by `sqrt(6)`, and excludes the complementary outcome from its ontology (Supplemental §§S1 and S7). The paper's record algebra is instead a **31,104-dimensional operator space** and has 172 distinguished weighted directions (Supplemental §§S6–S8). The two spaces must not be identified automatically.

Regroup the ordinary three-link physical Hilbert space as

\[
\mathcal H_{\rm phys}=\mathcal H_A\otimes\mathcal H_B\otimes\mathcal H_C,
\quad (\dim \mathcal H_A,\dim \mathcal H_B,\dim \mathcal H_C)=(9,4,6),
\quad \dim\mathcal H_{\rm phys}=216,
\]

where `A=Q1Q2`, `B=q1q2`, and `C=Q3 f`. Let

\[
|\omega_d\rangle=d^{-1/2}\sum_{j=0}^{d-1}|jj\rangle.
\]

For the aligned support convention `|0>,|1>` at each `(3,2)` link, the unjoined input is

\[
|\Psi_{\rm in}\rangle=
|\omega_2\rangle_{Q_1q_1}
|\omega_2\rangle_{Q_2q_2}
|\omega_2\rangle_{Q_3 f}.
\]

Apply the **literal** physical projector `P3` to the qutrit pair `A`:

\[
\boxed{
(P_3\otimes I_{BC})|\Psi_{\rm in}\rangle
=\frac1{\sqrt6}
|\omega_3\rangle_{Q_1Q_2}
|\omega_2\rangle_{q_1q_2}
|\omega_2\rangle_{Q_3 f}
}.
\]

The selected outcome has probability `p0=1/6`, exactly as in the source's laboratory postselection comparison. This identity was checked as a **216-component vector**, not only as a reduced state.

## 2. A minimal physical CPTP completion

Define two Kraus operators on `A`:

\[
K_0=P_3,\qquad K_f=I_9-P_3.
\]

Both are Hermitian projectors and

\[
K_0^\dagger K_0+K_f^\dagger K_f=P_3+(I_9-P_3)=I_9.
\]

The instrument maps on the actual physical density matrix are

\[
\mathcal I_i(\rho)=(K_i\otimes I_{BC})\rho(K_i^\dagger\otimes I_{BC}),
\qquad i\in\{0,f\}.
\]

Each map is completely positive, their sum is trace-preserving, and the maps have a standard finite-dimensional Stinespring dilation after adjoining an outcome register. The physical interpretation of that outcome register is an **additional assumption**, not a consequence of the source's no-failure ontology.

**No signaling, exact.** For any joint density operator `rho_AB`, and any local CPTP map `Phi_A`,

\[
\operatorname{Tr}_A[(\Phi_A\otimes\mathrm{id}_B)(\rho_{AB})]
=\operatorname{Tr}_A\rho_{AB}.
\]

For the aligned input, Bob's unconditioned qubits remain

\[
\boxed{\rho_B=I_4/4}.
\]

This does not hold branch by branch; outcome-conditional steering is fully compatible with no-signaling when outcome information is unavailable remotely.

### A universal constraint on the failure branch

Let `P_B=|omega2><omega2|` on the remote qubits. The successful unnormalized remote state is `(1/6)P_B`; the unconditional state is `I4/4`. Therefore the **coarse-grained** failure branch is constrained for *every* local trace-preserving completion retaining exactly the successful Kraus `K0=P3`:

\[
\begin{aligned}
\rho_{B,f}^{\rm unnormalized}&=\frac14I_4-\frac16P_B,\\
\operatorname{tr}\rho_{B,f}^{\rm unnormalized}&=\frac56,\\
\boxed{\rho_{B|f}}&=\boxed{\frac3{10}I_4-\frac15P_B}.
\end{aligned}
\]

Its eigenvalue on `|omega2>` is `1/10`; the other three are `3/10`. Thus its Bell fidelity is `1/10`. The aggregate failure state is fixed **even though individual failure outcomes and Alice's postmeasurement states are not**. This is a genuine nontrivial consequence of the input state plus local CPTP preservation of the exact selected effect. It is not a result of the source axiom by itself, because that axiom disallows the failure branch.

## 3. An explicit nine-Kraus completion that is pair-preserving on the aligned input

Define the occupied two-qutrit Bell ket

\[
|b_0\rangle=(|00\rangle+|11\rangle)/\sqrt2,
\qquad |J\rangle=|\omega_3\rangle,
\qquad |\eta\rangle=\sqrt3|b_0\rangle-\sqrt2|J\rangle.
\]

Complete the qutrit-pair orthonormal basis with the other three two-level Bell kets and the four basis vectors `|02>,|12>,|20>,|21>`. This gives nine `|v_j>` with `|v_0>=|J>`. Define

\[
K_j=|J\rangle\langle v_j|,\qquad j=0,\ldots,8.
\]

Then `K0=|J><J|=P3` exactly and `sum_j K_j†K_j=I9`. On the aligned input, the outcome probabilities are exactly

\[
\boxed{(1/6,\;1/12,\;1/4,\;1/4,\;1/4,\;0,0,0,0)}.
\]

Each of the five nonzero branches prepares `|omega3>` on `A`, a pure maximally entangled two-qubit state on `B`, and leaves the unread link `C` unchanged. The remote reduced state averaged over all outcomes remains `I4/4`.

**Scope limitations:** This completion uses a replacement/repreparation operation excluded by the source's Corollary S2. The preservation of entangled pairs holds for the particular aligned input; it is not a universal statement for all support orientations. The existence of this completion therefore proves compatibility with *one added operational theory*, not uniqueness or derivation from the original axiom.

## 4. Gauge-frame test: a conditional Schur-lemma obstruction

A fixed POVM outcome effect `E` on `A=C3⊗C3` is invariant under an independent local frame group if

\[
(g_1\otimes g_2)E(g_1^\dagger\otimes g_2^\dagger)=E
\quad \forall(g_1,g_2)\in SU(3)\times SU(3).
\]

The representation of the product group on `C3⊗C3` is irreducible. Schur's lemma implies

\[
\boxed{E=cI_9}.
\]

Thus a **fixed**, individually invariant, finite-outcome measurement for that entire independent product group can have no rank-one Bell effect. In particular, `P3` fails this requirement. By contrast,

\[
(U\otimes\overline U)|\omega_3\rangle=|\omega_3\rangle
\quad\Rightarrow\quad
(U\otimes\overline U)P_3(U^\dagger\otimes U^{\mathsf T})=P_3.
\]

This correlated transformation is compatible with the original Bell join. A **covariant family** of effects `E_g=g E g†`, transforming with its frame labels, is also valid, but specifying these labels/identifications is the additional structural question.

**Interpretive caution:** The source explicitly includes conjugate matchings on *each maximally entangled pair*. It does not establish that independent `SU(3)` transformations on the two qutrit copies must be treated as gauge identifications of the *fixed already-joined projector*. The Schur result is therefore conditional on the gauge-action hypothesis and is not an unconditional inconsistency theorem about the source.

### Exact finite Weyl test

Let `X|j>=|j+1 mod 3>`, `Z|j>=omega^j|j>` and `omega=e^(2πi/3)`. The nine states

\[
|J_{ab}\rangle=(X^aZ^b\otimes I_3)|\omega_3\rangle,\qquad a,b\in\{0,1,2\},
\]

are orthonormal. Hence

\[
\boxed{\frac19\sum_{a,b}|J_{ab}\rangle\langle J_{ab}|=\frac19I_9}.
\]

For the aligned physical input, the frame-averaged selected effect has probability `1/9` and leaves Bob with conditional state `I4/4`. The original frame-fixed selected effect has probability `1/6` and leaves Bob's conditional state `P_B`. The **difference is an operational distinction** between a known relational frame and an average over an unknown frame; it is not a prediction of variable `alpha`.

## 5. Support-orientation and frame-transport stress tests

For 64 seeded independent Haar-random qutrit support pairs, the physical two-outcome instrument retained remote no-signaling within `1.34e-16` Frobenius norm. Under simultaneous transport of physical inputs and the local instrument, the covariance residual was below `3.80e-16`. For the **fixed**, untransported `P3`:

- selected branch probability varied across the sampled orientations from `0.0836223` to `0.1587641`;
- remote pure-state concurrence varied from approximately `0.11737` to `0.99876`;
- the exact support choice `span(e0,e2)` versus `span(e0,e1)` had success probability `1/12` and a **product**, concurrence-zero remote state.

These computations are consistent with the v7 fixed-join no-go. Frame covariance is restored when the *instrument is transported together with the state*, but that transport has to be specified.

## 6. The 172-record constant does not emerge as a physical probability

The original source calculates

\[
\alpha_{\rm cand}=\frac{79+48\sqrt2+42\sqrt3+3\sqrt6}{31104}
=0.00729735256366226\ldots,
\]

as a **normalized trace in the operator readout space**. Our physical selected Bell outcome has probability `1/6`, not `alpha_cand`. Neither the two-outcome Lüders instrument nor the nine-Kraus completion singles out the record-weight sum as a probability or as a photon pole residue.

For clarity, it is straightforward but circular to reproduce the numerical value as a physical event probability. Conditional on Bell success, toss an independent classical coin with acceptance probability

\[
q=6\alpha_{\rm cand}\approx 0.04378411538.
\]

A three-outcome physical instrument can implement this via

\[
K_{\rm accept}=\sqrt q\,P_3,\qquad
K_{\rm reject-click}=\sqrt{1-q}\,P_3,\qquad
K_{\rm other}=I_9-P_3.
\]

Its acceptance probability is `q*(1/6)=alpha_cand`. However, `q` was **chosen using the desired answer**. Indeed any `t∈[0,1/6]` can replace `alpha_cand` by setting `q=6t`. Such matching has zero explanatory force.

One could instead attach numerical scores to the measurement branches and tune their weighted expectation. Again, a canonical score assignment connecting the source's 172 operator-space records to ordinary physical measurement outcomes is **not provided**. The source's trace is a well-defined algebraic number, not yet an experimentally derived Born probability.

## 7. Implications for the QED beta function

The first coefficient for a single unit-charge Dirac fermion,

\[
\beta(\alpha)=\frac{2}{3\pi}\alpha^2+O(\alpha^3),
\]

requires a physical electromagnetic two-point response, fermionic degrees of freedom, and scale dependence. Neither the source nor our v13 CPTP completion supplies these. Deriving the correct `alpha(0)` from a record count and deriving its running are logically separate problems.

**Pass/fail milestones for v14:**

1. Define which of the source's matching operations are physical gauge equivalences of the *complete joined state and instrument*, including all three shared-qutrit copies.
2. Construct an operational map from the 31,104-dimensional readout operator space to physical measurement effects/score distributions on `H_phys` that respects complete positivity and gauge transport.
3. Require that this map fix the 172 weights and bare normalization without inputting `alpha_cand` as a calibration parameter.
4. Only after these requirements, introduce dynamical time, charged fermions and an independently selected scale transformation needed for the beta function.

## 8. Reproducibility and assertions

Run `python chain_beta_extension_v13.py`; dependencies: Python 3, NumPy, SymPy. The script performs:

- exact selected join identity on **all 216 physical amplitudes**;
- exact Kraus completeness, branch probabilities and coarse failure density matrix;
- exact Bell concurrence checks for the aligned nine-Kraus completion;
- exact cyclotomic (`Q(sqrt(-3))`) verification that nine qutrit Weyl rotations form an orthonormal Bell basis, plus exact correlated matching invariance;
- 64 seeded random independent-support tests, remote no-signaling and physical-frame-transport covariance;
- an exact support-misalignment counterexample evaluated numerically;
- parameter-calibration nonuniqueness.

An independent regression run of `chain_beta_extension_v2.py` through `chain_beta_extension_v12.py` reported **11/11 successful process exits**. That verifies code execution, **not** every physical interpretation or theorem claimed in prior notes. The script emits `chain_beta_extension_v13_results.json`; regression summaries are in `chain_beta_extension_v13_regressions.json`.

## Conclusion

A local physical CPTP completion exists with the source Bell projector preserved exactly as a selective outcome. That completion restores no-signaling and fixes one nontrivial quantity—the **coarse-grained remote failure state**—but it is not unique, is not supplied by the source's no-failure ontology, and does not reproduce the 172-record number without an additional observable and score definition. Full independence of qutrit-frame changes is incompatible with a fixed rank-one Bell effect; correlated frames or a transported instrument can evade this conditional obstruction. The independent derivation of the QED beta-function coefficient remains **open**.