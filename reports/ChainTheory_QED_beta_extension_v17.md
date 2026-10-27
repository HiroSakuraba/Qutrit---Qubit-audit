# Chain Theory → QED beta research, v17
## Full factorized gauge-intertwiner classification and a non-scalar completely positive counterexample

**Date:** 25 October 2026  
**Primary source:** A. D. Wissner-Gross, *A Parameter-Free Candidate for the Fine-Structure Constant from Entangled Qudits* (23 October 2026), `ChainTheoryAlpha.pdf`, Letter Eqs. (1)–(6), Supplement S3, S6, S7, S8, S11.  
**Prior research:** v9 (172-of-175 read-arm model), v11 (31,104-dimensional embedding), v14 (physical operator-frame bridge), v15 (local effects), v16 (correlated-SU(3) invariant effect classification).  
**Status:** Exact representation-theoretic classification for one specified *factorized local group action*, and a constructive non-uniqueness theorem for covariant unital completely positive maps. **Not** a derivation of the electromagnetic coupling or the one-loop QED beta coefficient from the original axiom.

## Abstract and result

The paper supplies a 31,104-dimensional operator readout space with 172 selected directions, weights and multiplicities `(79,48,42,3)`, and the scalar candidate

\[
\alpha_{\rm cand}=\frac{79+48\sqrt2+42\sqrt3+3\sqrt6}{31104}.
\]

We equip the **complete** readout space with the natural factorized local unitary-conjugation representation induced by its three qutrit–qubit links. Under the group `G=(SU(3)×SU(2))^3`, the source Hilbert space splits into 16 inequivalent irreducible representations. A Heisenberg completely positive map from `B(R)` to local two-qutrit observables `B(C³⊗C³)` is equivariant precisely when its Choi matrix belongs to the commutant of the combined gauge action. The resulting Choi representation has **36 irreducible block types** and a **144-dimensional Hermitian commutant**. Normalized trace preservation imposes 16 independent linear conditions. Because the depolarizing channel has strictly positive Choi matrix, the admissible equivariant UCP maps contain a **128-dimensional affine interior**.

Two explicit kinds of counterexample demonstrate that the freedom is physical and not just a parameter-count artifact.

1. **Gauge-invariant scalar maps:** two distinct full-matrix-algebra UCP maps preserve all five record-sector trace fractions, including the zero sector, while differing on other full-algebra inputs.
2. **Non-scalar maps:** an explicit gauge-equivariant, unital, normalized-trace-preserving CP map preserves all five sector trace fractions and produces a different record-weighted expectation on the paper's original unjoined two-qutrit state. Its strictly positive Choi eigenvalues are calculated from a 27-dimensional Casimir factor, not guessed or inferred from numerics.

The second example is particularly decisive: under these augmented axioms the original weighted trace is **not uniquely mapped to the physical two-qutrit expectation**, even after full-space equivariance is imposed.

### Scope of the gauge-action assumption

This is a **source-motivated operational completion**. The source describes a basis matching on each maximally entangled pair and the operator readout factors, but **does not define an autonomous physical channel from the entire matrix algebra `B(R)` to measurements on the two joining qutrits**. We therefore take the standard local adjoint action on each slot's scalar/traceless operator Hilbert space. It contains independent local basis changes that are *more general* than pair-state-preserving matchings. Our explicit intertwiners commute with this enlarged action and consequently with its matching subgroups, but the calculation does **not** prove that this is the unique operational gauge action intended in the paper. The normalization of record eigenvalues as physical scores remains an extra hypothesis.

## 1. Complete source factorization and 16 representation sectors

In the v11 embedding the source operator Hilbert space is

\[
\begin{aligned}
\mathcal R&=\mathcal C_{Q_3 f}\otimes\mathcal L_{Q_1}\otimes\mathcal L_{q_1}
             \otimes\mathcal L_{Q_2}\otimes\mathcal L_{q_2},\\
\dim\mathcal C_{Q_3 f}&=(3^2-1)(2^2-1)=24,\\
\mathcal L_{Q_i}&={\rm End}(\mathbb C^3)_{\rm HS}=\mathbf1\oplus\mathbf8,\\
\mathcal L_{q_i}&={\rm End}(\mathbb C^2)_{\rm HS}=\mathbf1\oplus\mathbf3,\\
\dim\mathcal R&=24(1+8)(1+3)(1+8)(1+3)=31104.
\end{aligned}
\]

`C` is the connected unread-link sector, with the representation `(8,3)` under the unread `SU(3)×SU(2)`. Every read slot transforms in its identity or adjoint sector under an **independent** local factor of `G`. Because the group factors act independently, each choice of the four identity/adjoint labels gives an irreducible tensor-product representation inequivalent to all others.

Write `b=(b_Q1,b_q1,b_Q2,b_q2)∈{0,1}^4`, with `0` for scalar and `1` for traceless. Its invariant projector `Π_b` has rank

\[
\boxed{d_b=24\,8^{b_{Q_1}}3^{b_{q_1}}8^{b_{Q_2}}3^{b_{q_2}}.}
\]

There are 16 such irreducibles, and the code verifies `Σ_b d_b=31104` exactly. In particular, the all-scalar read-arm sector has rank 24, while the four-traceless sector has rank 13,824.

The v11 selected record space uses only one-slot qubit directions and two-slot traceless product directions; it has no component in either of those extreme sectors. This is checked directly against every canonical record column, with **exact rational overlaps** `Tr(Π_b P_w)` for all 16 patterns and all five weight sectors `w=0,1,√2,√3,√6`. The per-weight sums reproduce `(30932,79,48,42,3)`.

## 2. Full Choi block classification

The output Hilbert space of the two joining qutrits is

\[
\mathcal H_A=\mathbb C^3_{Q_1}\otimes\mathbb C^3_{Q_2},\qquad\dim\mathcal H_A=9,
\]

transforming under independent fundamental representations `U_{Q1}⊗U_{Q2}`. For a Heisenberg map

\[
\Phi:B(\mathcal R)\to B(\mathcal H_A),
\]

choose the Choi convention

\[
J_\Phi=\sum_{ij}|i\rangle\langle j|\otimes\Phi(|i\rangle\langle j|)
\quad\in B(\overline{\mathcal R}\otimes\mathcal H_A).
\]

Full gauge equivariance is equivalent to `[J_Φ,\overline{R(g)}⊗W(g)]=0` for every `g∈G`. The relevant local `SU(3)` decomposition is

\[
\boxed{(\mathbf1\oplus\mathbf8)\otimes\mathbf3
       =2\,\mathbf3\oplus\overline{\mathbf6}\oplus\mathbf{15}.}
\]

The unread representation stays `(8,3)`; each qubit read factor stays `1⊕3`. Thus the full Choi representation has `3×3×2×2=36` **inequivalent irrep types**. For each joining qutrit the multiplicities of its types `(3,bar6,15)` are `(2,1,1)`; the read-qubit and unread factors have multiplicity one. Hence the complex dimension of the commutant equals

\[
\boxed{\dim_\mathbb C\operatorname{Comm}(\overline R\otimes W)
=4(2^2+1^2+1^2)^2=144.}
\]

The real vector space of Hermitian invariant Choi matrices also has dimension 144. In a group-adapted basis **every** equivariant completely positive map has Choi form

\[
\boxed{J_\Phi=\bigoplus_{\lambda=1}^{36}
I_{d_\lambda}\otimes H_\lambda,\qquad H_\lambda\succeq0,}
\]

where each PSD multiplicity block has size 1, 2, or 4. This is the complete representation-theoretic **cone classification for the chosen group**, up to the ordinary unital/trace conditions; it is not an enumeration of 279,936-square matrices.

Normalized trace preservation means

\[
\frac{1}{9}\operatorname{Tr}_{\mathcal H_A}\Phi(X)
=\frac{1}{31104}\operatorname{Tr}_{\mathcal R}X
\qquad\forall X,
\]

or equivalently

\[
\boxed{\operatorname{Tr}_{\mathcal H_A}J_\Phi=
\frac9{31104}I_{\mathcal R}.}
\]

Because the input representation has 16 inequivalent multiplicity-one sectors, the invariant left side has 16 independent scalar components. These 16 linear constraints are independent: invariant Choi operators `Π_b⊗I_9` have partial output traces `9Π_b`. Thus the normalized-trace affine space has dimension

\[
\boxed{144-16=128.}
\]

Unitality follows from equivariance and normalized trace preservation: `Φ(I_R)` commutes with the independent fundamental `SU(3)×SU(3)` output action, so is scalar `cI_9`; its trace equals 9, forcing `c=1`. The completely depolarizing map's Choi matrix `I/(31104)` is positive definite. Therefore the 128-dimensional affine constraint space intersects the CP cone in a relative open set, not only on a degenerate boundary.

**Interpretation:** this settles the algebraic dimension of the full gauge-equivariant Choi family **conditional on the stated factorized action**. The original source does not establish a unique physical choice among its members.

## 3. Two different invariant scalar channels with exactly the same record statistics

Let `P_w^R` denote the paper-motivated full-R spectral projectors onto the five weight sectors. Define

\[
\sigma_0=\frac{I_R}{N},\qquad
\sigma_t=\frac{I_R}{N}+t\left(\frac{\Pi_{0000}}{24}
 -\frac{\Pi_{1111}}{13824}\right),\quad N=31104.
\]

For `t=1/4`, the spectrum of `σ_t` is positive: it equals `325/31104` on `Π0000`, `7/497664` on `Π1111`, and `1/31104` everywhere else. Its trace is one; it commutes with **all** independent local matching transformations.

Construct the Heisenberg channels

\[
\Phi_{\sigma}(X)=\operatorname{Tr}(\sigma X)I_9.
\]

Each is unital, completely positive, and fully equivariant (the output is a scalar operator). They are **distinct** on the full record matrix algebra, e.g. on `Π0000`. Yet every nonzero `P_w^R` is orthogonal to the two extreme source sectors, and the zero sector's overlap cancels. Hence exactly

\[
\boxed{\operatorname{Tr}(\sigma_tP_w^R)=\frac{n_w}{N}
\quad\text{for all five }w.}
\]

Therefore neither gauge covariance nor exact matching to the five source trace fractions fixes the full channel. Scalar channels alone do not, however, establish different physical *scores* because these two examples agree on the spectral-label subalgebra. For that, use the non-scalar counterexample below.

## 4. Constructing a non-scalar full-gauge intertwiner

Choose Hermitian traceless Gell-Mann generators `ℓ_a (a=1,...,8)` normalized by

\[
\operatorname{tr}(\ell_a\ell_b)=\delta_{ab}.
\]

Let `F_a` be their Hermitian generators in the 8-dimensional adjoint representation of `SU(3)` acting on traceless qutrit operators. Pad each `F_a` by zero on the scalar 1-dimensional summand to obtain `9×9` Hermitian matrices, also called `F_a`. The exact Lie-algebra normalization is

\[
\boxed{\operatorname{Tr}_8(F_aF_b)=6\delta_{ab}.}
\]

Reordering only tensor *factors*, write the source space as

\[
\mathcal R\cong(\mathbf1\oplus\mathbf8)_{Q_1}\otimes
(\mathbf1\oplus\mathbf8)_{Q_2}\otimes\mathbb C^{384}_{\rm rest},
\]

where `384=24×4×4`. Define

\[
A_{ab}=\frac1{384}\left(F_a\otimes F_b\otimes I_{384}\right),\qquad
B_{ab}=\ell_a\otimes\ell_b.
\]

For any real `ε`, set

\[
\boxed{\Phi_\varepsilon(X)=\frac{\operatorname{Tr}X}{31104}I_9
+\varepsilon\sum_{a,b=1}^8\operatorname{Tr}(A_{ab}X)B_{ab}.}
\tag{17.1}
\]

The adjoint covariance of the `F_a` and the fundamental conjugation covariance of the `ℓ_a` imply

\[
\Phi_\varepsilon(R(g)XR(g)^\dagger)
=W(g)\Phi_\varepsilon(X)W(g)^\dagger.
\]

The other local factors act trivially on the `I384` in the source probes and on the output, so the relation holds for all six independent local group factors. Because `Tr Aab=Tr Bab=0`, the map is unital **and** preserves normalized trace for every input, not only for the five record projectors.

### 4.1 Complete positivity proved from the Choi spectrum

The Choi matrix of (17.1) factors, up to a tensor-factor permutation, as

\[
J_{\Phi_\varepsilon}=\frac{I_{9N}}{N}
+\frac{\varepsilon}{384}\bigl(K\otimes K\otimes I_{384}\bigr),
\qquad K=\sum_aF_a^{\mathsf T}\otimes\ell_a.
\]

The `27×27` matrix `K` has spectrum

\[
\boxed{\operatorname{spec}(K)=
\{(-1)^{\times15},0^{\times3},1^{\times6},3^{\times3}\}.}
\]

One derivation uses `8⊗3 = 3⊕bar6⊕15` and quadratic Casimirs `C_2(3)=4/3`, `C_2(bar6)=10/3`, `C_2(15)=16/3`, `C_2(8)=3`, noting the transpose of the purely imaginary adjoint generators. An independent numerical diagonalization reproduces these eigenvalues and multiplicities. The paired product eigenvalues range from `−3` to `9`.

Since `N=81×384`, all Choi eigenvalues are

\[
\boxed{\lambda(J)=\frac{1+81\varepsilon\lambda_i\lambda_j}{N}.}
\]

For real `ε`, this gives the exact CP interval

\[
\boxed{-\frac1{729}\le\varepsilon\le\frac1{243}.}
\]

In particular `ε=1/1000` is strictly positive, with

\[
\boxed{\lambda_{\min}(J)=\frac{757}{31\,104\,000}>0,\qquad
\lambda_{\max}(J)=\frac{1729}{31\,104\,000}.}
\]

The calculation is exact once the `SU(3)` representation decomposition is established. No huge Choi matrix or speculative numerical positivity criterion is needed.

## 5. Record effects, source-state expectation, and non-uniqueness

The first record block contains 64 `Q1Q2` letters. Its symmetric (`36`-dimensional) and antisymmetric (`28`-dimensional) pieces receive the source weights `√3` and `1`. The adjoint generators have zero trace, so all other record blocks have vanishing contraction with `Aab`. On the `Q1Q2` symmetric/antisymmetric projectors the traces are respectively `+3δab/384` and `−3δab/384`. Therefore the five *physical* effects are exactly

\[
\begin{aligned}
E_0&=\frac{30932}{N}I_9,\\
E_1&=\frac{79}{N}I_9-\frac{3\varepsilon}{384}B,\\
E_{\sqrt2}&=\frac{48}{N}I_9,\\
E_{\sqrt3}&=\frac{42}{N}I_9+\frac{3\varepsilon}{384}B,\\
E_{\sqrt6}&=\frac3N I_9,\qquad
B=\sum_a\ell_a\otimes\ell_a.
\end{aligned}
\]

The `SU(3)` Fierz identity gives

\[
\boxed{B=S-\frac13I_9,}
\]

where `S|a,b>=|b,a>` is the qutrit swap. The effects are positive because they are images of positive projectors under the CP channel, sum to identity, and preserve every normalized physical trace `Tr(E_w)/9=n_w/N`.

Consequently the physical score observable is

\[
\boxed{\Phi_\varepsilon(M)=\alpha_{\rm cand}I_9+
\frac{\varepsilon(\sqrt3-1)}{128}
\left(S-\frac13I_9\right).}
\tag{17.2}
\]

For the source's original **unjoined** two-qutrit state,

\[
\rho_A=\left[\frac12\operatorname{diag}(1,1,0)\right]^{\otimes2},\quad
\operatorname{Tr}(\rho_AS)=\frac12,\quad
\operatorname{Tr}(\rho_AP_3)=\frac16,
\]

so

\[
\boxed{\langle w\rangle_{\rho_A,\varepsilon}
=\alpha_{\rm cand}+\frac{\varepsilon(\sqrt3-1)}{768}.}
\]

At `ε=0` this is the source's target merely because we chose uniform trace evaluation. At `ε=1/1000`, it is

\[
\begin{aligned}
\alpha_{\rm cand}&=0.007297352563662260\ldots,\\
\langle w\rangle_{\varepsilon=1/1000}&=0.007298305754817949\ldots,\\
\Delta&=\frac{\sqrt3-1}{768000}=9.5319115569\times10^{-7}.
\end{aligned}
\]

Thus two UCP **full** gauge intertwiners obey the same normalized-trace and source-record multiplicity constraints and have different physical scores. This is not an independent prediction that the real electromagnetic constant varies: `ε` labels alternative, unselected operational bridges. The paper still has no physical prescription identifying *any* of these weighted scores with Thomson-limit electromagnetism.

A Bell-first refinement is algebraically possible in the aligned basis, because all five effects commute with `P3`. Choosing `K_Bell=P3` and failure Kraus operations `K_w=√E_w (I−P3)` is trace-preserving and preserves the selected branch. Its *failure-only weighted score* is not Eq. (17.2); giving the Bell branch score zero invokes the v16 obstruction. More generally, transporting the instrument with the gauge frame is distinct from requiring each fixed outcome effect to be separately invariant under the independent group.

## 6. Numerical and regression audit

`chain_beta_extension_v17.py` performs the following reproducible checks:

- All 16 source representation-sector ranks and their exact dimension sum.
- The 36 Choi irrep types, dimension sum `9×31104=279936`, invariant commutant dimension 144, and normalized affine dimension 128.
- Exact rational overlaps of all 172 record columns with the scalar/traceless sectors, including the 30,932-dimensional complement.
- Exact positivity and record-sector calibration of two different invariant scalar density functionals.
- `SU(3)` adjoint generator normalization, the `K` Casimir spectrum and strictly positive CP Choi eigenvalue bounds for `ε=1/1000`.
- The `SU(3)` Fierz identity `Σℓ_a⊗ℓ_a=S−I/3`.
- Positivity, completeness, calibration, and Bell-commutation of the five physical effects.
- Twelve independently sampled `SU(3)×SU(3)` transformation pairs on a complete `81`-dimensional joining-qutrit-operator subfactor, using random Hermitian `81×81` test operators. The other `384`-dimensional factor is handled **analytically** by tensor identities, not by a random scan of the huge matrix algebra.
- Independent `64×64` qutrit adjoint-sector contractions for three generators.
- All 16 research scripts v2–v17 ran successfully (`chain_v17_regression.json`, log).

The maximum observed random transformation residual for the unscaled intertwiner on test matrices was `1.18×10^-13` in Hilbert–Schmidt norm; with the physical CP perturbation it was `7.41×10^-16`. These are **supporting checks** of an algebraic covariance result, not evidence of uniqueness.

## 7. What was established, and what was not

**Established in the specified operational completion**

1. Exact decomposition of the full source readout representation, with 16 irreducible sectors.
2. Exact 36-isotypic Choi classification, 144-dimensional Hermitian commutant, and a 128-dimensional affine family after normalized-trace conditions.
3. Multiple distinct, gauge-equivariant, unital, completely positive full-algebra maps respecting the record trace calibration.
4. An explicit non-scalar map changing the weighted score on the source's physical input, with a strict analytic CP certificate and checked covariance.
5. A Bell-first completion can preserve the source's exact selected projector, but does not resolve the physical scoring rule.

**Not established**

1. That the source's own axioms determine the factorized operational gauge action uniquely.
2. That an arbitrary 172-record score is a physical probability, scattering amplitude, or observable electromagnetic coupling.
3. That the source's single Bell-projection outcome defines a trace-preserving measurement without added operations.
4. That a microscopic spacetime, Dirac fermion, physical scale, or renormalization-group flow arises from distillation.
5. That `β(α)=(2/3π)α²+...` follows independently from the qudit axiom.

### Recommended v18 research direction

The gauge-intertwiner uniqueness route has now yielded a conditional, constructive **non-uniqueness theorem**. Another round of representation-specific POVM searches will not establish physical QED running. Future work should either (a) prove an operational selection principle deriving a **particular** CP map and outcome scoring rule from the axiom, or (b) switch to an explicit spacetime/scale construction and show how propagating charged states emerge without importing Dirac dynamics. The latter is the actual missing prerequisite for a genuinely independent beta-function coefficient.

### References

- A. D. Wissner-Gross, *A Parameter-Free Candidate for the Fine-Structure Constant from Entangled Qudits* (uploaded 23 October 2026 paper), especially Supplemental S6–S8 and S11.
- Standard `SU(3)` conventions for normalized Gell-Mann matrices and structure constants; see the CERN-hosted notes `https://s3.cern.ch/inspire-prod-files-e/e445fb8205e1f10e186e4407a6dfd4cf`.
- The standard `SU(3)` completeness/Fierz identity, e.g. *European Physical Journal C* article `https://link.springer.com/article/10.1140/epjc/s10052-026-15366-z`, Eq. (66); our identity uses `ℓ_a=λ_a/√2`.
- Companion files `chain_beta_extension_v17.py`, `chain_beta_extension_v17_results.json`, `chain_v17_regression.json`.