# Chain Theory → QED beta program — v16
## Complete classification of invariant two-qutrit outcome effects, Choi-certified CP extensions, and a Bell-first score no-go

**Date:** October 25, 2026  
**Primary source:** A. D. Wissner-Gross, *A Parameter-Free Candidate for the Fine-Structure Constant from Entangled Qudits* (October 23, 2026), `ChainTheoryAlpha.pdf`, especially Letter Eqs. (2)–(6) and Supplement S7–S8.  
**Prior numerical construction:** v11 record embedding, v14 unital CP frame bridge, v15 localized swap-algebra effects and exact Bell-first instrument.  
**Status:** Fully proved exact statements within an *explicitly supplemented operational model*. No derivation of the physical coupling or QED beta function from the original axiom.

## Abstract

Let the original record weights be indexed by `w=0,1,√2,√3,√6`, with multiplicities `(30932,79,48,42,3)` in the paper's `N=31104`-dimensional operator readout space. Fix a nine-dimensional physical two-qutrit Hilbert space `H_A=C³⊗C³`, Bell projector `P=|ω₃><ω₃|`, and `F=I-P`. Impose **fixed output-side invariance** of physical effects under *all correlated* `W(U)=U⊗conj(U)` for `U∈SU(3)`; impose the conventional POVM positivity and normalization; and impose an additional **spectral-trace calibration** `Tr(E_w)/9=n_w/N`. The complete classification of these effects is `E_w=a_w P+b_w F`, `a_w,b_w≥0`, `Σa=Σb=1`, and `a_w+8b_w=9n_w/N`. In particular, none of these conditions uniquely determines a measurement.

An explicit unital completely positive extension from the full record matrix algebra to physical effects exists for each member of this family. It is described by a positive, factored Choi matrix and does not require allocating a `279936×279936` array. This establishes existence but emphatically not uniqueness of a *full* gauge-equivariant bridge under the original, incompletely specified action on the 31,104 record directions.

For the source's original unjoined input, `p=Tr(ρ_A P)=1/6`. The attainable expectation value of the five record labels ranges over `[15α_cand/16,3α_cand/2]`, including the candidate by choice rather than necessity. A **stronger conditional no-go** holds if one requires the distinct Bell outcome `K_Bell=P`, scores the record labels only on the complementary branches, and assigns Bell outcome score zero: the maximum weighted expectation is `15α_cand/16`, strictly smaller than `α_cand`. Achieving the latter by adding a Bell-outcome score needs `s_B≥3α_cand/8`, an extra calibration.

## 1. What comes from the paper and what is added

From the source paper, retain only the exact record weights, counts, the proposed normalized trace, and the Bell projector used in its selected join:

\[
N=31104,\quad(n_0,n_1,n_{\sqrt2},n_{\sqrt3},n_{\sqrt6})=(30932,79,48,42,3),
\]
\[
\alpha_{\mathrm{cand}}=\frac{79+48\sqrt2+42\sqrt3+3\sqrt6}{31104}
=0.00729735256366226025\ldots.
\]

The source **does not** require a five-outcome POVM on two physical qutrits, the trace-calibration condition, an outcome-scoring rule, or trace preservation of the selected Bell projection. These are additional choices explored here to test whether they could make the proposed number physical. Nor does it define the full action of independent link matchings on each 31,104-dimensional record projector as an operational measurement.

The physical two-qutrit space has dimension 9, and we write
\[
P=|\omega_3\rangle\langle\omega_3|,\qquad F=I_9-P,
\]
with ranks 1 and 8. The source's original, unjoined qutrit state, after tracing out the paired qubits, is
\[
\rho_A=\left(\tfrac12\operatorname{diag}(1,1,0)\right)^{\otimes2},
\quad \operatorname{tr}(\rho_A P)=\tfrac16.
\]

## 2. Classification theorem (fixed output-side invariance)

**Theorem 1.** Let `{E_w}` be a finite POVM on `H_A` invariant under every `W(U)=U⊗\overline U`, `U∈SU(3)`. Then
\[
E_w=a_wP+b_wF,\quad a_w,b_w\in\mathbb R_{\ge0},
\quad\sum_w a_w=\sum_w b_w=1.
\]
Conversely, every such family is a positive, normalized, fixed-invariant POVM.

**Proof.** The representation of `SU(3)` on `C³⊗\overline{C³}` is the direct sum of its one-dimensional singlet and its eight-dimensional adjoint irreducible representation. Its commutant is `span(P,F)`. Any invariant Hermitian effect is therefore `aP+bF` by Schur's lemma. Because these are mutually orthogonal projectors, positivity is equivalent to `a,b≥0`. Summing effects to identity enforces `Σa=Σb=1`. Sufficiency is immediate. □

**Extra calibration assumption.** Require the normalized physical trace of each effect to equal the source's record-subspace trace fraction:
\[
\frac{\operatorname{tr}E_w}{9}=\frac{n_w}{N}.
\]
Then
\[
\boxed{a_w+8b_w=\mu_w:=\frac{9n_w}{N}.}
\]
For `w≠0`, positivity allows `0≤a_w≤μ_w`; these parameters can be chosen independently because `Σ_{w≠0} μ_w=9\times172/31104=43/864<1`. Set
\[
a_0=1-\sum_{w\ne0}a_w,\quad
b_w=\frac{\mu_w-a_w}{8}.
\]
All coefficients are nonnegative, and `Σb=1` follows from `Σμ=9` and `Σa=1`. Thus the calibrated family has **four independent continuous parameters**, corresponding to the four nonzero record labels. This is an exact characterization of the spectral-outcome effects under the stated constraints, not a classification of all CP superoperators on the full source algebra.

## 3. A factored Choi certificate for explicit unital CP bridges

Write `P_w^R` for the mutually orthogonal eigenspace projectors of the 31,104-dimensional source record operator `M`, with ranks `n_w`. Define the **Heisenberg** map
\[
\boxed{\Phi(X)=\sum_w \frac{\operatorname{Tr}(P_w^R X)}{n_w}\,E_w,}
\quad\Phi:B(\mathbb C^{31104})\longrightarrow B(\mathbb C^9).
\]
This is a concrete extension from the full input matrix algebra. Because each normalized trace `Tr(P_w^R X)/n_w` is a positive linear functional and every `E_w` is positive, `Φ` is completely positive. It is unital since `Φ(I_R)=ΣE_w=I_9`. It maps `P_w^R` precisely to `E_w`, while discarding all input coherences between source-weight sectors.

Its Choi matrix in the convention `J_Φ=Σ_ij |i><j|⊗Φ(|i><j|)` is exactly
\[
\boxed{J_\Phi=\sum_w \frac{(P_w^R)^{\mathsf T}}{n_w}\otimes E_w\succeq0.}
\]
In an eigenbasis of the input record operator, each of the `n_w` blocks has one eigenvalue `a_w/n_w` and eight eigenvalues `b_w/n_w`. They are all nonnegative. Moreover,
\[
\operatorname{Tr}_{\mathrm{input}}J_\Phi=I_9,
\]
which checks unitality in the Choi representation. There is no reason to allocate its square matrix of dimension 279,936 to establish this factorization.

This is a **Choi-certified family**, not a uniqueness result. For example, the `a_w=0` extreme has Choi rank 279,764, whereas the `a_w=μ_w` extreme has rank 278,560; the uniform scalar map has full Choi rank 279,936. All three satisfy the same output symmetry, trace calibration, locality of the output effects, and positivity.

**Gauge-action qualification.** The map is fixed invariant on the *output* under correlated qutrit transformations. To claim actual gauge equivariance for the **untransported** input record algebra, one needs the precise source representation `R(g)` and must verify `Φ(R(g)XR(g)†)=W(g)Φ(X)W(g)†`. The source has not established this for all record sectors, and neither have we. Defining a transported family `Φ_g` gives covariance by construction, not uniqueness of a fixed physical observable.

## 4. Exact score interval on the source's original input

For a physical density matrix `ρ` define `p=Tr(ρP)`. A weighted score measured using the invariant POVM is
\[
\mathbb E_\rho[w]=\sum_{w\ne0} w\big[p a_w+(1-p)b_w\big].
\]
Using `b_w=(μ_w-a_w)/8` yields
\[
\mathbb E_\rho[w]=\frac{1-p}{8}\sum_{w\ne0}w\mu_w
+\frac{9p-1}{8}\sum_{w\ne0}w a_w.
\]
For the source input `p=1/6` and `Σ w μ_w=9α_cand`, this becomes
\[
\boxed{\mathbb E_{\rho_A}[w]
=\frac{15}{16}\alpha_{\mathrm{cand}}
+\frac1{16}\sum_{w\ne0}w a_w.}
\]
Every term is nonnegative and `0≤a_w≤μ_w`; hence the **sharp**, attained range is
\[
\boxed{\frac{15}{16}\alpha_{\mathrm{cand}}
\le\mathbb E_{\rho_A}[w]
\le\frac32\alpha_{\mathrm{cand}}.}
\]
The endpoints are obtained by setting every positive-label `a_w=0` or `a_w=μ_w`, respectively. There is no selection of the candidate within that interval.

**Three explicit examples:**

1. **Outside-Bell extreme:** `a_w=0` for each positive label; weighted score on the original input is `0.006841268028433369`.
2. **State-independent scalar POVM:** `E_w=(n_w/N)I_9` for every `w`; weighted score is identically `α_cand` for every physical state. This reproduces the target solely through the imposed trace calibration and therefore supplies no dynamical electromagnetic explanation.
3. **Inside-Bell extreme:** `a_w=μ_w` for each positive label; weighted score on the original input is `0.010946028845493390`.

All three are positive, normalized, and invariant under the full correlated `SU(3)` matching action on the two joining qutrits. They are mathematically different measurements.

## 5. Exact Bell-first score no-go under trace calibration

In v15 we required a **distinct**, literal Bell outcome `K_Bell=P`. For an invariant family of effects we may refine its complement with Kraus operators
\[
K_w=\sqrt{E_w}\,F=\sqrt{b_w}\,F, \qquad
\sum_wK_w^\dagger K_w=F,\qquad K_{\mathrm{Bell}}^\dagger K_{\mathrm{Bell}}=P.
\]
The resulting local instrument is completely positive and trace preserving. It satisfies no-signaling when outcomes are not communicated. On the source input its Bell probability remains `p=1/6`.

If the Bell outcome has score **zero**, and the four positive record weights are scored only on the failure outcomes, the total expectation is
\[
\mathbb E_{\rm sequential}[w]
=(1-p)\sum_{w\ne0}w b_w
=\frac56\sum_{w\ne0}w\frac{\mu_w-a_w}{8}.
\]
Since `a_w≥0`,
\[
\boxed{0\le\mathbb E_{\rm sequential}[w]
\le\frac56\frac98\alpha_{\mathrm{cand}}
=\frac{15}{16}\alpha_{\mathrm{cand}}
<\alpha_{\mathrm{cand}}.}
\]
This upper bound is sharp and is achieved by the outside-Bell extreme.

If the source's candidate is to be recovered by assigning a nonzero *new* numeric score `s_B` to the selected Bell outcome, the most favorable possible failure score still implies
\[
\frac16 s_B+\frac{15}{16}\alpha_{\mathrm{cand}}
\ge\alpha_{\mathrm{cand}},
\quad\boxed{s_B\ge\frac38\alpha_{\mathrm{cand}}
=0.002736507211373348\ldots.}
\]
This is an **additional scoring convention**. The source axiom does not select it.

The impossibility statement is deliberately narrow: it uses the specified state, distinct selected Bell branch, zero score for that branch, and individual multiplicity-based trace calibration. Dropping any of those requirements may evade it. It does not prove that the original paper's normalized trace is mathematically inconsistent; it proves only that this particular physical interpretation cannot reproduce it.

## 6. Exact matched-frame averaging of the v15 effects

The v15 local effects were `c_w I+d_w S`, where `S` swaps the two qutrits. They commute with `P` but fail fixed invariance under all `W(U)=U⊗\bar U`.

The exact `SU(3)` Haar twirl on the output transforms the swap as
\[
\boxed{\mathcal T(S)=P+\frac14 F.}
\]
**Derivation:** `S|ω₃>=|ω₃>` and `Tr S=3`. The invariant operator must have coefficient 1 on `P` and coefficient `(3−1)/8=1/4` on `F`.

Therefore the twirled v15 effects are
\[
\mathcal T(E_w)=(c_w+d_w)P+(c_w+\tfrac14 d_w)F.
\]
They remain positive and trace calibrated, now have fixed correlated-`SU(3)` invariance, and remain non-unique relative to the scalar POVM and other invariant families.

On the source input their stand-alone weighted score is `0.007307649381702107`, whereas the scalar POVM gives `0.007297352563662260`; their Bell-first zero-score values are `0.006063965772985473` and `0.006081127136385217`, respectively. In particular, **twirling to restore fixed-frame invariance does not force the original target value on the source state**.

## 7. Implementation and verification

Script: `chain_beta_extension_v16.py`, writing `chain_beta_extension_v16_results.json`.

* Exact SymPy derivation of `P,F`, source input `p=1/6`, all five effects and traces, complete classification constraints, sharp weighted-score endpoints, and the Bell-first no-go.
* Explicit completely positive, unital extension with a factored positive Choi matrix. Choi dimensions and ranks certified algebraically, without attempting memory-prohibitive dense allocation.
* Exact `SU(3)`-twirled v15 swap effects and their weighted expectations.
* Eighty seeded random `SU(3)` matched-frame tests. Worst fixed-effect covariance/invariance residual `3.42×10⁻¹⁵` and worst `P3` residual `1.09×10⁻¹⁵`.
* Clear accounting for assumptions not present in the paper.

No beta-function calculation is attempted here, because the paper still specifies no intrinsic momentum/energy scale, causal charged-matter dynamics, or polarization spectral measure. Importing the one-loop QED coefficient would not answer the original research question.

## 8. What would settle the next question?

The **next strongest test** should explicitly construct the physical input–output matching action on *every* sector of the original `31104`-dimensional record algebra, not merely the labels `P_w`, and solve the intertwiner/Choi constraints under that full action. It must also identify a physically privileged scoring operator from the distillation axiom. If the only forced physical bridges remain depolarizing/trace-evaluation maps, matching `α_cand` will be a mathematical identity without predictive dynamics. If a distinguished nontrivial bridge survives, calculate its expectation on the paper's initial state and assess whether it is compatible with the `P3` branch without modifying record multiplicities.

A scale-dependent propagating theory remains indispensable for obtaining the first QED beta-function coefficient independently. The v16 result concerns measurement uniqueness, not energy running.