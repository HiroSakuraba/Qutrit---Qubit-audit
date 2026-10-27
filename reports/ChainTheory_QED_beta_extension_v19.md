# Chain Theory → QED beta program, v19
## Positive overlap transfers from Bell distillation, gauge-covariant holonomy, and nonunique continuum limits

**Date:** 25 October 2026  
**Primary source:** A. D. Wissner-Gross, *A Parameter-Free Candidate for the Fine-Structure Constant from Entangled Qudits* (23 October 2026), `ChainTheoryAlpha.pdf`.  
**Previous results:** v2–v18; especially v7 on misaligned supports, v10 on relational connection, and v18 on graph geometry and missing propagation.  
**Status:** Exact finite-dimensional conditional theorems and numerical checks, **not** a derivation of physical time, a preferred continuum theory, or the QED beta function.

### Abstract

For two two-dimensional supports in a common qutrit space, the overlap matrix `K=V_i† V_j` supplies a contractive, compositional operator and a completely positive **trace-nonincreasing** map `ρ↦KρK†`. With a matching conjugate-frame convention, the singular values of `K` determine the outcome of the source paper's qutrit Bell projection. They are exactly `(1,c)` with `c` the modulus of the inner product of the unit vectors normal to the two supports. The Bell probability is `(1+c²)/12`, and the selected remote state's concurrence is `2c/(1+c²)`. The result recovers the paper's `1/6` probability when the supports align, but for `c<1` the output ceases to be maximally entangled, so the candidate map is **not** universally a permissible state transition under the source's pair-only ontology.

A chosen sequence of overlaps has an exact, positive subchannel composition. Under local changes of the two-dimensional frames, transfer matrices transform covariantly. In an example of three real rank-two supports in `C³`, the polar part of a closed product has *nontrivial* and exactly calculable holonomy `tr W=16/sqrt(65)`, without assigning an independent gauge matrix to any edge. This is a conditional positive construction: its common ambient reference and chosen path are additional relational inputs, not furnished by the source axiom.

The decisive negative result concerns the continuum: smooth infinitesimal changes of support produce asymptotically unitary transfer, while a sequence of zigzagging supports, converging *uniformly to the same constant support*, yields `diag(1,exp(-σ²/2))`. Hence a continuum map is **not determined by the uniform limiting support configuration**. A quadratic-variation/regularity prescription is needed. No physical time, scale, charged particle spectrum, or `2/(3π)` coefficient follows.

## 1. Primary-source boundary

The Letter defines a maximally entangled `(d1,d2)` link as a pair ket on its `m=min(d1,d2)` common support (Eq. 1). In the source's unusual ontology, the same shared qudit is represented by distinct **virtual copies** on its several links, and distillation selects a Bell projection with no admissible failure branch (Supplement S1, Corollaries S1–S5). The specifically selected charged local object uses two read links, an unread third link, and a qutrit shared hub (Theorem 1). In Supplement S7 the Bell join on two aligned qutrit–qubit links has ordinary Born success probability `1/6`; the model conditions on this branch and normalizes by `sqrt(6)`.

These are the source ingredients used below. The source does **not** specify: (i) a common physical identification of two independently supplied qutrit support subspaces, (ii) an ordered path or foliation, (iii) a time-evolution generator, (iv) a trace-preserving physical instrument, (v) a Dirac fermion, or (vi) a momentum-dependent propagator (see Supplement S8). Each addition is therefore marked as a *conditional extension* rather than a logical consequence of the axiom.

## 2. Exact two-link contraction and principal-angle theorem

Let `V,W:C²→C³` satisfy `V†V=W†W=I₂`. We can represent a link state by

\[
 |\psi(V)\rangle=\frac1{\sqrt2}\sum_{a=0}^1 V|a\rangle\otimes|a\rangle.
\]

Take the fixed qutrit Bell bra

\[
 \langle\omega_3|=\frac1{\sqrt3}\sum_{r=0}^2\langle r,r|.
\]

For **real** supports `V,W`, direct contraction of the two source pair kets gives a remote two-qubit coefficient matrix

\[
 C=\frac1{2\sqrt3}V^{\mathsf T}W.
\tag{19.1}
\]

For general complex supports use the conjugate representation on the second qutrit leg, consistent with the Bell pairing's `U⊗Ū` invariance. The coefficient matrix is then the complex conjugate of `K=V†W`, multiplied by `1/(2sqrt(3))`, so all singular-value results remain unchanged. **The requirement to identify which second-leg frame is conjugate is part of the relational-comparison convention.** It cannot be inferred from two unlabelled input rays alone.

**Theorem 19.1 (two planes in three dimensions).** Let `n_V,n_W` be unit vectors spanning the respective one-dimensional orthogonal complements of `im V, im W`. The singular values of `K=V†W` are

\[
 \boxed{\operatorname{sv}(K)=\bigl(1,c\bigr),\qquad c=|\langle n_V,n_W\rangle|\in[0,1].}
\tag{19.2}
\]

**Proof.** Two rank-two subspaces of a 3-dimensional vector space intersect in dimension at least one. Their common unit vector gives singular value 1. For the other squared singular value, use `WW†=I−|n_W><n_W|` to write

\[
 KK^\dagger=V^\dagger WW^\dagger V=I_2−V^\dagger |n_W\rangle\langle n_W|V.
\]

The nontrivial eigenvalue is `1−||V†n_W||²=|<n_V,n_W>|²=c²`, as required. QED.

**Corollary 19.2 (Bell success and output entanglement).** With the source's qutrit Bell bra,

\[
 \boxed{p=\|C\|_{\mathrm F}^2=\frac{1+c^2}{12},\qquad
 \mathcal C_{\mathrm{remote}}=\frac{2c}{1+c^2}.}
\tag{19.3}
\]

The normalized two-qubit Schmidt amplitudes are `1/sqrt(1+c²)` and `c/sqrt(1+c²)`. Consequently `c=1` gives a Bell pair and `p=1/6`, whereas `c=0` gives a product state and `p=1/12`. The latter is not an allowed maximally entangled output of the source axiom. One must restrict relative supports or supplement the distillation rule; the raw fixed Bell projector does not resolve this internally.

The Python script checks (19.1) by a literal four-index tensor contraction, verifies (19.3) at four angles, and checks (19.2) on 128 independent complex random isometry pairs. Maximum singular-value error was less than `8×10⁻¹⁶`.

## 3. A completely positive compositional **subchannel**, not a physical clock

Assume the common ambient comparison has been specified. Set `K_ij=V_i†V_j` and

\[
 \boxed{\mathcal T_{ij}(\rho)=K_{ij}\rho K_{ij}^\dagger.}
\tag{19.4}
\]

For any positive `ρ`, the output is positive. Moreover `K_ij†K_ij≤I`, because `V_i V_i†` is an orthogonal projector. Thus the operation is completely positive and trace-nonincreasing (CPTNI). For an ordered sequence,

\[
 \mathcal T_{12}\circ\mathcal T_{23}\circ\cdots\circ\mathcal T_{n-1,n}(\rho)
 =K_{\mathrm{path}}\rho K_{\mathrm{path}}^\dagger,
 \quad K_{\mathrm{path}}=K_{12}\cdots K_{n-1,n}.
\tag{19.5}
\]

The composition order reflects the convention that `K_ij` maps the coordinate column at endpoint `j` into the coordinate column at endpoint `i`. This is an **ordered contraction**, not an elapsed physical time. Under reparametrizations `V_i→V_i g_i`, with `g_i∈U(2)`, one has

\[
 K_{ij}\to g_i^\dagger K_{ij}g_j,
 \qquad K_{\mathrm{path}}\to g_{\mathrm{start}}^\dagger K_{\mathrm{path}}g_{\mathrm{end}}.
\tag{19.6}
\]

For `c<1`, normalizing a selected outcome makes it a nonlinear transformation of density matrices:

\[
 \widehat{\mathcal T}(\rho)=\frac{K\rho K^\dagger}{\operatorname{tr}(K\rho K^\dagger)}.
\tag{19.7}
\]

For example take `K=diag(1,c)`, `0<c<1`. The normalized outputs of basis projectors `|0><0|` and `|1><1|` are the same projectors, but the output of their equal mixture is `diag(1,c²)/(1+c²)`, not `I/2`. This difference is exact and is not statistical noise.

A standard trace-preserving completion adds another Kraus operator

\[
 K_{\rm fail}=\sqrt{I-K^\dagger K},
 \qquad K^\dagger K+K_{\rm fail}^\dagger K_{\rm fail}=I.
\tag{19.8}
\]

This is a legitimate conventional quantum instrument, but its complementary branch is **not** among the source's permitted outcomes. It is a candidate operational completion, not part of Chain Theory's original axiom.

## 4. A nontrivial geometric holonomy, conditional on support comparisons

In the real ambient qutrit basis `(e₀,e₁,e₂)`, choose

\[
 V_0=(e_0,e_1),\quad
 V_1=(e_0,c e_1+s e_2),\quad
 V_2=(c e_0+s e_2,e_1),
 \quad c=\frac{\sqrt3}{2},\;s=\frac12.
\tag{19.9}
\]

These are three legitimate rank-two support isometries. Their overlap matrices satisfy

\[
 K_{01}=\begin{pmatrix}1&0\\0&c\end{pmatrix},\quad
 K_{12}=\begin{pmatrix}c&0\\s^2&c\end{pmatrix},\quad
 K_{20}=\begin{pmatrix}c&0\\0&1\end{pmatrix}.
\tag{19.10}
\]

The closed-loop product is

\[
 \boxed{H=K_{01}K_{12}K_{20}=
 \frac34\begin{pmatrix}1&0\\1/4&1\end{pmatrix}.}
\tag{19.11}
\]

Its unitary polar factor is the **nontrivial** rotation

\[
 \boxed{W=\frac1{\sqrt{65}}
 \begin{pmatrix}8&-1\\1&8\end{pmatrix},\qquad
 \operatorname{tr}W=\frac{16}{\sqrt{65}}\approx1.9845557534.}
\tag{19.12}
\]

The polar rotation angle is `atan(1/8)≈0.124355` rad. Under independent local changes of support frames, the loop matrix is conjugated at its basepoint, and the polar rotation's trace is invariant. The exact SymPy identities and 64 randomly chosen U(2) frame triples verify this claim.

**Interpretation.** Once a common qutrit ambient identification and a closed sequence of supports are provided, overlap-derived relative transports can carry nontrivial geometric holonomy without supplying independent gauge matrices on the edges. This is a mathematical construction related to parallel transport on the Grassmannian. It does *not* show the original pair-only axiom specifies the common ambient embedding, selects this loop, or turns its holonomy into electromagnetic field strength. For singular `H`, the polar unitary would also require a nonunique extension.

## 5. The central v19 no-go: uniform support geometry does not fix continuum transfer

Take smooth real support frames

\[
 V(\theta)=\bigl(e_0,\cos\theta\,e_1+\sin\theta\,e_2\bigr).
\tag{19.13}
\]

Then

\[
 K(\theta_i,\theta_{i+1})
 =\operatorname{diag}\bigl(1,\cos(\theta_{i+1}-\theta_i)\bigr).
\tag{19.14}
\]

**Smooth path.** Set `θ_j=jΘ/n`. The product is

\[
 K_n^{\rm smooth}=\operatorname{diag}\bigl(1,[\cos(\Theta/n)]^n\bigr)
 \underset{n\to\infty}{\longrightarrow}I_2,
\tag{19.15}
\]

because `n log cos(Θ/n)=−Θ²/(2n)+O(n⁻³)→0`. The ordered product transports the moving frame without attenuation in this particular real gauge. More generally, a smooth `V(x)` gives a path-ordered unitary continuum limit generated by the anti-Hermitian connection `A(x)=V(x)†V'(x)`. That limit depends on a chosen smooth support path, not a physical Hamiltonian.

**Uniformly convergent zigzag path.** Let `n` be even, fix `σ>0`, and set

\[
 \theta_j^{(n)}=\begin{cases}0,&j\text{ even},\\ \sigma/\sqrt n,&j\text{ odd} .\end{cases}
\tag{19.16}
\]

Every step has angle magnitude `σ/sqrt(n)`, and `sup_j|θ_j^{(n)}|=σ/sqrt(n)→0`. Therefore the *entire support path converges uniformly* to the constant support `V(0)`, with exactly the same endpoint support for every even `n`.

Nevertheless,

\[
 K_n^{\rm zigzag}
 =\operatorname{diag}\bigl(1,[\cos(\sigma/\sqrt n)]^n\bigr)
 \underset{n\to\infty}{\longrightarrow}
 \boxed{\operatorname{diag}\bigl(1,e^{-\sigma^2/2}\bigr)}.
\tag{19.17}
\]

This follows from `log cos x=−x²/2+O(x⁴)`; the microscopic quadratic variation satisfies `Σ_j(Δθ_j)²=σ²` for all `n`. Compare the constant approximants `θ_j=0`, which converge to the *same* limiting support path and have transfer `I₂` at every `n`.

**Theorem 19.3 (continuum non-uniqueness).** There is no map from the uniform limit of the support path to the overlap-transfer limit that is independent of the approximating sequence, even in this smooth ambient two-plane family. Two sequences of polygonal supports converge uniformly to the same constant path but their matrix transfer limits differ by `diag(0,1−e^{-σ²/2})`.

The example directly proves the claim. QED.

A nonzero dissipative limit thus requires extra regularity or a choice of microscopic quadratic-variation law. A stochastic law could supply such a limit, but that would be an independent new axiom or emergent statistical theorem, not something fixed by the Bell-pair definition.

For `σ=1.2`, the limiting second singular value is `exp(−0.72)=0.4867522559599717`, and the normalized endpoint concurrence is approximately `0.787`. The numerical products at `n=20,80,320,1280,5120` converge to the exact value. In contrast the smooth path with `Θ=0.9` approaches unit gain.

These limits are akin to projection-induced geometric/quantum-Zeno transport; the broader Zeno literature analyzes limits of repeated projections only under additional assumptions about the projection family and time parameter. Relevant background: Facchi and Pascazio, *Quantum Zeno dynamics: mathematical and physical aspects*, arXiv:0903.3297; S. Becker, N. Datta, and R. Salzmann, *Quantum Zeno Effect in Open Quantum Systems*, Annales Henri Poincaré **22**, 3795–3840 (2021), doi:10.1007/s00023-021-01075-8. The present result is an elementary exact counterexample, not an invocation of their theorems.

## 6. Why this does not derive the QED beta function

The source computes an exact normalized trace for a fixed 172-record counting problem,

\[
 \alpha_{\rm cand}=
 \frac{79+48\sqrt2+42\sqrt3+3\sqrt6}{31104}.
\tag{19.18}
\]

The new positive transfer construction has **no derived energy variable** `μ`, no continuum spacetime dimension, no interacting charged Dirac field, no gauge-field propagator, and no vacuum-polarization renormalization scheme. Neither the loop holonomy of §4 nor the microscopic attenuation of §5 identifies an electromagnetic coupling; assigning the observed beta coefficient to their parameters would be a fit or an assumption.

The first QED coefficient for one electrically charged unit-charge Dirac fermion remains

\[
 \beta(\alpha)=\frac{2}{3\pi}\alpha^2+O(\alpha^3),
\]

**not derived by this qudit construction**. To test a genuine derivation, one must specify a propagation/continuum rule from physical observables, show a charged fermion spectrum and a massless gauge sector, compute a gauge-transverse two-point response in that scheme, and extract its universal logarithmic coefficient *without introducing the coefficient as an input*.

## 7. Reproducibility and falsifiable next steps

Executable test file: `chain_beta_extension_v19.py`. Output: `chain_beta_extension_v19_results.json`. Dependencies: NumPy and SymPy. It performs 27 named assertions, including exact three-plane loop and polar-unitary calculations, generic complex principal-angle and covariance trials, CP trace-nonincreasing checks, Born-projection contractions, explicit Kraus completeness, nonlinearity, finite overlap products, and smooth versus zigzag continuum convergence. Random tests use the recorded deterministic seed `20261025`. `chain_beta_extension_v18.py` was also rerun successfully during v19 development.

The next research milestone should **not** be to postulate another lattice Hamiltonian. Instead:

1. **Intrinsic alignment:** Check whether allowed matching equivalences plus the source's actual shared-copy join determine the overlap `V_i†V_j` uniquely, or whether it is extra geometric data. Analyze gauge-frame stabilizers and unoccupied qutrit directions.
2. **Operational closure:** Require every selected branch to return allowed maximally entangled pairs. For the candidate map, the exact condition is `c=1`; test whether source-compatible adaptive projections evade this condition without an external relative frame.
3. **Continuum selection:** Identify a source-internal principle that fixes the regularity/quadratic variation of support paths. Theorem 19.3 is an explicit falsification test for proposals depending only on the uniform limiting geometry.
4. **Physical scaling:** Only after (1)–(3), ask whether transfer eigenvalues relate to spectral dimension, momentum, a charged propagator, and a renormalization-group flow.

**Bottom line:** Positive compositional transfer, nontrivial geometric holonomy, and exact continuum counterexamples are available in an explicit support-overlap completion. The original axiom does not select the completion or the continuum regularization. The QED beta function remains unestablished.