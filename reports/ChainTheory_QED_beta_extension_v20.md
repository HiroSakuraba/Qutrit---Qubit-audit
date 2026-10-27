# Chain Theory → QED beta program, v20
## Same Bell distillation, inequivalent covariant transfer laws: a conditional dynamical-selection obstruction

**Date:** 25 October 2026  
**Primary source:** A. D. Wissner-Gross, *A Parameter-Free Candidate for the Fine-Structure Constant from Entangled Qudits* (23 October 2026), `ChainTheoryAlpha.pdf`.  
**Previous work:** v2–v19, particularly v13 (CPTP completion) and v18–v19 (network geometry, overlap transfer, continuum ambiguity).  
**Status:** Exact finite-dimensional theorem and reproducible tests. **Not** an emergent Hamiltonian, a proof that no extension can ever work, or a derivation of the QED beta coefficient.

## Abstract

The source paper's selected join is the rank-one projector `P=|ω₃><ω₃|` on two qutrit read copies. It assigns the normalized state resulting from this single branch and explicitly excludes the orthogonal complement as an admissible configuration. Section S7 additionally says that the model has no time order. These features make the step from a distillation rule to propagation a genuine open question.

We prove an exact **conditional non-uniqueness theorem**: there is a continuum of completely positive, trace-preserving (CPTP), local, Bell-branch-preserving, correlated-`SU(3)`-covariant channels on the joining pair. Their selected Kraus operation is **identically** `ρ ↦ PρP`, but their action on the orthogonal eight-dimensional sector differs. Their nontrivial discrete relaxation eigenvalue runs over `[0,1]`. In a simple optional Poisson-clock construction, the relaxation rate can be any `γs`, with `γ>0` introduced by hand. Thus the selected branch, locality, full trace preservation, and this Bell-stabilizer covariance do **not** choose physical propagation, decay, or an intrinsic time unit.

The counterexample is about **ordinary quantum-mechanical operational completions**, not about different allowed failure configurations in the paper: the source's ontology forbids such configurations. Its original 172-record normalized trace remains unchanged. No electron, photon propagator, or one-loop QED beta coefficient emerges.

## 1. What the primary source specifies—and what it does not

The source's Supplement S1 states that ordinary entanglement swapping has a Bell-measurement outcome with a Born probability, whereas its model keeps only the selected projection, with **no failure branch**. Supplement S7 gives the joined state

\[
|\Omega\rangle=\sqrt6\,(P_3)_{Q_1Q_2}
|\omega_2\rangle_{Q_1q_1}
|\omega_2\rangle_{Q_2q_2}
|\omega_2\rangle_{Q_3f}
\]

and states that this joint readout **has no time order**. In conventional quantum mechanics the selected `P₃` branch has probability `1/6` on this aligned input. The third, unread qutrit–qubit link is spectator throughout the following calculation.

The original rules do not provide a CPTP failure completion, iteration frequency, time coordinate, routing among network junctions, physical energy, fermion action, or electromagnetic two-point function. In particular, we cannot take a convenient conventional completion and silently call it a new consequence of the original axioms.

## 2. Exact theorem: an entire family of symmetry-respecting physical completions

Let

\[
\mathcal H_A=\mathbb C^3\otimes\mathbb C^3\cong\mathbb C^9,
\qquad |\omega_3\rangle=\frac1{\sqrt3}\sum_{a=0}^2|a,a\rangle,
\qquad P=|\omega_3\rangle\langle\omega_3|,
\qquad F=I_9-P.
\]

The complement `F` has rank eight. For each `s∈[0,1]`, define the linear map

\[
\boxed{\mathcal E_s(X)=PXP+(1-s)FXF+
\frac{s}{8}\operatorname{Tr}(FXF)F.}\tag{20.1}
\]

**Theorem 20.1 (operational non-uniqueness).** For every `s∈[0,1]`, Eq. (20.1) is a CPTP map on the joining qutrit pair with the following properties:

1. It admits an instrument having **exactly** `K_{\rm Bell}=P` as a distinguished Kraus operator, so the probability and state of this branch are unmodified.
2. It acts only on the joining qutrits, hence cannot signal to remote qubits when outcomes are not conditioned on.
3. It is covariant under every correlated matching-frame transformation `W(U)=U⊗\overline U`, `U∈SU(3)`, even with the physical measurement frame held fixed: `\mathcal E_s(WXW†)=W\mathcal E_s(X)W†`.
4. Its 81-dimensional superoperator spectrum is

\[
\boxed{\operatorname{spec}(\mathcal E_s)=
\{1^{\times2},\;(1-s)^{\times63},\;0^{\times16}\},}
\tag{20.2}
\]

where degeneracies coincide at the endpoints. Different `s` give inequivalent maps.

**Proof (CPTP and preservation of the source branch).** Choose any orthonormal basis `{|f_i>}_{i=1}^8` for the range of `F`. Use the Kraus operators

\[
\boxed{
 K_{\rm Bell}=P,\qquad
 K_{\rm hold}=\sqrt{1-s}\,F,\qquad
 K_{ij}=\sqrt{s/8}\,|f_i\rangle\langle f_j|,
 \quad 1\le i,j\le8.
}\tag{20.3}
\]

The 64 `K_ij` satisfy

\[
\sum_{ij}K_{ij}^{\dagger}K_{ij}=sF,
\]

because each input basis vector is counted eight times. Hence

\[
\sum_a K_a^{\dagger}K_a=P+(1-s)F+sF=I_9.
\]

This proves complete positivity and trace preservation. Summing the corresponding operations gives (20.1). The selected Kraus operator is identically `P`, independent of `s`. `K_{\rm hold}` and the `K_{ij}` are additional operational failure branches; the original theory does not admit them as physical configurations. QED.

**Proof (covariance).** The maximally entangled ket obeys

\[
(U\otimes\overline U)|\omega_3\rangle=|\omega_3\rangle.
\]

Thus `W P W†=P`, `W F W†=F`. The first two terms of (20.1) transform by conjugation; the scalar `Tr(FXF)` is invariant. This gives exact covariance for the **unconditioned channel and its coarse-grained Bell/failure instrument**. An individual labelled `K_ij` need not be invariant: the failure Kraus decomposition itself rotates. QED.

**Proof (spectrum).** Split `End(C⁹)` into four blocks relative to `C⁹ = im P ⊕ im F`:

* Bell population `PXP`, dimension 1, eigenvalue 1;
* scalar population `Tr(FXF)F/8`, dimension 1, eigenvalue 1;
* traceless operators within the 8-dimensional complement, dimension `8²−1=63`, eigenvalue `1−s`;
* the off-diagonal Bell/complement coherences `PXF+FXP`, complex dimension `2×8=16`, eigenvalue 0.

The dimensions sum to `1+1+63+16=81`. QED.

**Scope of the covariance claim.** The theorem concerns the correlated `SU(3)` stabilizer of the original Bell join, not the independent full `SU(3)×SU(3)` action and not the unspecified full matching representation of all 31,104 readout directions. It proves that even this nontrivial physically relevant set of constraints leaves operational freedom; it is not a complete classification of all source-theoretic dynamics.

## 3. Operationally inequivalent outputs: two exact witnesses

At `s=0`, Eq. (20.1) is the Lüders-type nonselective Bell measurement

\[
\mathcal E_0(X)=PXP+FXF,
\]

which retains all information inside `im F`. At `s=1`, it is

\[
\mathcal E_1(X)=PXP+\operatorname{Tr}(FXF)F/8,
\]

which erases all fine-grained information inside `im F`.

For any unit vector `|f>∈im F`, the input `ρ_f=|f><f|` gives

\[
\mathcal E_0(\rho_f)=\rho_f,
\qquad
\mathcal E_1(\rho_f)=F/8.
\]

Their trace distance is **exactly**

\[
\boxed{\frac12\|\rho_f-F/8\|_1=\frac78.}\tag{20.4}
\]

This follows from eigenvalues `7/8` and `−1/8` (sevenfold) on `im F`.

The difference is also observable on the source's own unjoined input. Tracing out the read qubits gives the local density matrix

\[
\rho_A=\operatorname{diag}(1/4,1/4,0,1/4,1/4,0,0,0,0),
\]

in the computational basis `|00>,|01>,...,|22>`. Its two output channels differ by a Hermitian matrix whose **exact** eigenvalues are

\[
\{7/48\; (3\text{ times}),-1/48\;(1\text{ time}),-5/48\;(4\text{ times}),0\;(1\text{ time})\}.
\]

Consequently,

\[
\boxed{\frac12\|\mathcal E_0(\rho_A)-\mathcal E_1(\rho_A)\|_1=\frac7{16}.}\tag{20.5}
\]

Unlike an arbitrary test state in the failure complement, (20.5) concerns the local marginal of the source's actual two read links. This was independently checked with exact SymPy rationals.

## 4. Original Bell outcome and no-signaling: verified without changing the source input

Reorder the original read-link pure state as local two-qutrit coordinates versus remote two-qubit coordinates. Its coefficient matrix is

\[
C_{(3i+j),(2i+j)}=\tfrac12\quad(i,j\in\{0,1\}),
\]

with all other entries zero. The local projector produces an unnormalized remote density matrix `C†PC`, with

\[
\operatorname{Tr}(C^{\dagger}PC)=\frac16,
\qquad
\frac{C^{\dagger}PC}{1/6}
=|\omega_2\rangle\langle\omega_2|.
\tag{20.6}
\]

This result is **independent of `s`** because the distinguished selected Kraus operator does not change. Since every `E_s` is trace preserving on the joining pair, Bob's unconditional state stays `I₄/4`. The normalized remote state conditional on the **coarse-grained failure event** is also universal:

\[
\boxed{\rho_{B|\mathrm{fail}}=
\frac3{10}I_4-\frac15|\omega_2\rangle\langle\omega_2|.}
\tag{20.7}
\]

But the **local failure output** varies: Eq. (20.5) proves it. Therefore, fixing all the selected-branch and coarse-grained remote statistics still does not determine the full process.

The unread third qutrit–qubit pair is a tensor factor and is preserved by all these local instruments. Restoring it gives the original 216-dimensional three-pair physical system without altering these identities.

## 5. Optional iteration produces arbitrary relaxation

Iterating the unconditioned `E_s` for a fixed number `n≥1` sends an initial pure failure-sector state to

\[
\boxed{
\mathcal E_s^n(\rho_f)
=(1-s)^n\rho_f+\big[1-(1-s)^n\big]\,F/8.
}\tag{20.8}
\]

The distance to its maximally mixed fixed point is

\[
\boxed{D_n=\frac78(1-s)^n.}\tag{20.9}
\]

**Nothing in the original Bell outcome determines `s`**, so nothing determines this relaxation spectrum per arbitrarily labelled step.

To display the still larger clock ambiguity, *add* a Poisson jump clock with rate `γ>0` and define

\[
\boxed{
\Lambda_t^{(s,\gamma)}
=\exp\!\left[\gamma t(\mathcal E_s-\mathrm{id})\right]
=e^{-\gamma t}\sum_{n=0}^{\infty}
\frac{(\gamma t)^n}{n!}\mathcal E_s^n.
}\tag{20.10}
\]

Every term is CPTP and correlated-`SU(3)`-covariant, and the Poisson weights are a normalized nonnegative distribution. Thus (20.10) is a CPTP, covariant one-parameter quantum Markov semigroup. Its GKSL generator is explicitly

\[
\mathcal L_{s,\gamma}(\rho)
=\gamma\sum_aK_a\rho K_a^\dagger-\gamma\rho
=\gamma\sum_a\left(K_a\rho K_a^\dagger
-\tfrac12\{K_a^\dagger K_a,\rho\}\right),
\tag{20.11}
\]

where Eq. (20.3) supplies its jumps. For traceless `F`-sector components, the decay rate is **`γs`**; for Bell–failure coherence it is **`γ`**. Both `γ` and `s` are external choices. A change in `γ` rescales any imagined physical time, and `s` even changes the dimensionless ratio of the two relaxation rates.

This shows more than trivial lack of seconds: a dimensionless dynamical ratio remains undetermined. No preferred continuous-time limit or Hamiltonian has been derived. The semigroup is an illustrative completion, **not** a time evolution asserted to exist inside the source's no-time-order ontology.

## 6. Relationship to v19 and gauge-field emergence

Version 19 showed that support-overlap paths can encode nontrivial holonomy but have nonunique continuum transfer even under uniformly converging support geometry. Version 20 attacks the problem from the complementary direction: even at a **fixed, aligned, maximally entangled Bell join**, a continuum of distinct operational dynamics retains the same selected event and matching-frame covariance. The two ambiguities are logically independent:

1. **Geometric regularity ambiguity (v19):** relative support orientations alone do not select a unique continuum transfer without a fluctuation rule.
2. **Instrument/dynamics ambiguity (v20):** the selected Bell branch and a substantial symmetry requirement do not select a unique CPTP step or intrinsic time rate.

One cannot resolve the second simply by choosing a smoother support path, nor resolve the first simply by imposing trace preservation.

In the original source the numerical candidate

\[
\alpha_{\rm cand}=
\frac{79+48\sqrt2+42\sqrt3+3\sqrt6}{31104}
\]

is a **dimensionless normalized record trace**. It is not a function of a physical energy variable. Neither (20.1) nor (20.10) provides a photon, charged-fermion dispersion law, momentum-dependent polarization function, or a renormalization prescription. Hence they do not yield the first QED coefficient

\[
\beta(\alpha)=\frac{2}{3\pi}\alpha^2+O(\alpha^3)
\]

from the original source axioms. Changing `γ` or `s` cannot legitimately be fitted to that known coefficient and then called an independent derivation.

## 7. Reproducibility and exact/numerical checks

Run `python chain_beta_extension_v20.py`; it writes `chain_beta_extension_v20_results.json`. The code uses NumPy and SymPy, with a fixed random seed. Checks include:

- **Symbolic exact:** the Bell and failure projectors, the original-input trace-distance spectrum `7/16`, and the independent pure-failure distance `7/8`.
- **Explicit Kraus sums:** 65 operators on `C⁹`, including 64 matrix-unit jumps. Exact completeness identities and numerical residuals were checked at five values of `s`.
- **Physical source input:** the `1/6` Bell success, exact remote Bell state, universal remote coarse-grained failure density matrix, and non-signaling for three completions.
- **Covariance:** 160 sampled correlated `SU(3)` basis rotations, maximum Frobenius residual about `6.18×10⁻¹⁶`.
- **Full map spectrum:** numerical diagonalization of the `81×81` superoperator at `s=0,0.35,1`, matching the exact analytic multiplicities.
- **Iteration:** eight-step decay at four distinct `s` values.
- **Continuous-time toy completion:** composition and Poisson/Kraus expansions for four `(s,γ,t)` combinations, with semigroup composition residual under `8×10⁻¹⁷`.

The attached regression manifest records the current runs of v2–v20. A previously generated output JSON in this package should not be treated as an additional independent proof: it records the implementation's checks. The central claims have explicit algebraic proofs above.

## 8. Recommendation after v20

Continuing to add arbitrary quantum channels will not fix the missing physics. The next investigation should seek a **source-derived compositional principle** strong enough to select both *which links can act next* and *how their actions are weighted*. Candidate mathematical formulations include a causal partial order, an information-theoretic least-action rule, or a canonical quantum stochastic process. Each would be **a new postulate unless it can be shown to follow from the original pair/distillation definitions**.

A falsifiable milestone is to derive a finite-network transition operator with a uniquely fixed spectrum, then establish a controlled continuum scaling law with **no free dimensionless relaxation parameters**. Only after that would a charged matter/current sector and a QED running-coupling calculation become productive targets.

**Conclusion.** The original selected Bell projection is mathematically precise and survives every v20 test. Its use as an elementary dynamical law remains underdetermined: complete positivity, locality, no-signaling, and correlated-matching covariance do not fix propagation, attenuation, or time. This is a conditional no-uniqueness theorem about its physical completions, not a contradiction of the original source's no-failure ontology and not a derivation of electromagnetism.