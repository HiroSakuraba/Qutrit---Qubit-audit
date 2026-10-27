# Chain Theory and the first QED beta-function coefficient

**Research note — 24 October 2026**  
**Status:** Explicit finite-qudit gauge construction and conditional one-loop calculation; **not** a derivation from the original Chain Theory axiom alone.

## Abstract

We test whether the qutrit–qubit construction in A. D. Wissner-Gross, *A Parameter-Free Candidate for the Fine-Structure Constant from Entangled Qudits* (23 October 2026), determines the first quantum-electrodynamic (QED) beta-function coefficient. It does not: the original model defines a fixed normalized two-point readout but has neither momentum dependence nor a scale-transformation law. We construct a separately identified extension using two-qubit fermionic modes, a finite cyclic gauge-link qudit, and a four-dimensional gauge-covariant lattice Dirac operator. The finite gauge constraints are tested exactly up to floating-point precision. With **additional** continuum Dirac dynamics, the spinorial Seeley–DeWitt coefficient is independently reproduced using gamma-matrix traces and a Landau-level spectral sum. This gives the usual one-Dirac-fermion result, 

\[
\beta(e)=\mu\frac{de}{d\mu}=\frac{q^2e^3}{12\pi^2}+O(e^5),\qquad
\beta(\alpha)=\frac{2q^2}{3\pi}\alpha^2+O(\alpha^3).
\]

We also establish a **conditional non-uniqueness counterexample**: the same finite-link qudit gauge symmetry supports both the naive and Wilson four-dimensional discretizations. One has 16 light Dirac tastes, and the other has one in the continuum limit. Therefore gauge covariance and the finite-dimensional link algebra alone do not select the desired coefficient. Finally we calculate the charge-square and Dynkin-index factors contributed by the paper's \((\mathbf3,\mathbf2)_{1/6}\) representation, and reconstruct the Standard Model one-loop gauge coefficients **only after supplying the remaining Standard Model field content and loop formulas**.

## 1. Precisely what is supplied by the source

The source paper establishes a readout based on 172 records:

\[
\alpha_{\rm cand}=\frac{79+48\sqrt2+42\sqrt3+3\sqrt6}{31104},
\quad
\alpha_{\rm cand}^{-1}=137.035999189566155\ldots .
\]

It identifies this number with the zero-transfer electromagnetic coupling through its **Photon Conjecture** (main text, Theorem 2 and Conjecture 1). It identifies a central charge on its qutrit–qubit link and the group \((SU(3)\times SU(2)\times U(1))/\mathbb Z_6\), including the representation \((\mathbf3,\mathbf2)_{1/6}\) (Supplement S3). Its normalized record-algebra trace gives \(B(Q,Q)=1\) (Supplement S6).

Supplement S8 explicitly acknowledges that the proposed construction defines **no momentum-dependent propagator**, charged loop, vertex factor, or loop momentum. It does **not** derive an electroweak symmetry-breaking mechanism, physical spacetime, a fermion action, or a coarse-graining flow.

Hence the original value alone underdetermines the function \(e(\mu)\). The beta function of the unextended model is **undefined**, not zero.

### Formal underdetermination at the original level

For \(a=\alpha_{\rm cand}\), every smooth \(c\) produces a possible extension

\[
\alpha_c(Q)^{-1}=a^{-1}-\tfrac c2\ln(1+Q^2/m^2),
\]

which agrees at \(Q=0\) but satisfies

\[
Q\frac{d\alpha_c}{dQ} = c\alpha_c^2\frac{Q^2}{Q^2+m^2}
\longrightarrow c\alpha_c^2\quad (Q\gg m).
\]

This is a family of *mathematical interpolations*, not constructed physical theories; it proves only that the boundary value cannot pick a unique coefficient.

## 2. Added finite-qudit gauge dynamics

### 2.1 Why a finite qudit needs a cyclic gauge algebra

A unitary \(U\) and Hermitian \(E\) acting on a finite-dimensional Hilbert space cannot obey \([E,U]=U\): multiplying by \(U^\dagger\) and taking the trace gives \(0=\operatorname{tr}I=d\).

Instead let \(\omega=e^{2\pi i/K}\) and take a \(K\)-level link with basis \(|n\rangle\), \(n\in\mathbb Z_K\). Define

\[
Z|n\rangle=\omega^n|n\rangle,\qquad
X|n\rangle=|n-1\bmod K\rangle.
\]

Both are unitary and obey the exact finite Weyl relation

\[
XZX^\dagger=\omega Z.
\]

This provides an exact finite \(\mathbb Z_K\) gauge link. A suitable \(K\to\infty\) limit would be needed to recover ordinary continuous \(U(1)\) gauge fields; such a limit is **additional structure**.

### 2.2 A gauge-invariant matter–link Hamiltonian

Encode two spinless fermionic modes \(c_0,c_1\) by two qubits using a Jordan–Wigner map (the code verifies their canonical anticommutation relations). Define

\[
H=c_0^\dagger Zc_1+c_1^\dagger Z^\dagger c_0.
\]

With number operators \(n_j=c_j^\dagger c_j\), define local finite gauge transformations

\[
G_0=e^{-2\pi i n_0/K}\otimes X,\qquad
G_1=e^{-2\pi i n_1/K}\otimes X^\dagger.
\]

Then

\[
G_0HG_0^\dagger=H=G_1HG_1^\dagger,\quad
G_0G_1=G_1G_0,\quad G_j^K=I.
\]

This is a functioning, finite, exactly gauge-invariant gauge–matter **extension**. Its \(4K\)-dimensional Hilbert space is explicit. It does not purport to prove that the source axiom generates this Hamiltonian or the Jordan–Wigner fermion encoding.

**Tests:** \(K=3,5,7,11\), largest local gauge residual \(1.02\times10^{-15}\), as implemented in `chain_beta_extension_v2.py`.

## 3. A four-dimensional Dirac operator from a two-qubit spinor space

A four-component Euclidean spinor is algebraically representable by two qubits. With Pauli matrices \(\sigma_i\), choose

\[
\gamma_j=\sigma_1\otimes\sigma_j\;(j=1,2,3),\qquad
\gamma_4=\sigma_2\otimes I_2.
\]

Direct matrix multiplication verifies

\[
\{\gamma_\mu,\gamma_\nu\}=2\delta_{\mu\nu}I_4.
\]

This supplies a representation of \(\mathrm{Cl}_4(\mathbb C)\) **once four dimensions and fermion spin dynamics are chosen**. Two qubits do not by themselves explain four-dimensional spacetime or the existence of Dirac fermions.

On a hypercubic lattice, introduce unitary link phases \(U_{x,\mu}\) from \(\mathbb Z_K\). The lattice Wilson operator with spacing \(a\), bare fermion mass \(m\), and Wilson parameter \(r\) is

\[
\begin{split}
D_r(x,y)={}&\left(m+\frac{4r}{a}\right)\delta_{xy}I_4\\
&-\frac1{2a}\sum_{\mu=1}^{4}
\left[(rI_4-\gamma_\mu)U_{x,\mu}\delta_{x+\hat\mu,y}
+(rI_4+\gamma_\mu)U^*_{x-\hat\mu,\mu}\delta_{x-\hat\mu,y}\right].
\end{split}
\]

Under \(\psi_x\mapsto g_x\psi_x\) and \(U_{x,\mu}\mapsto g_xU_{x,\mu}g_{x+\hat\mu}^{-1}\),

\[
D_r[U^g]=G D_r[U]G^\dagger.
\]

Therefore the determinant is gauge invariant. This was independently verified numerically on periodic \(2^4\) and \(3^4\) site lattices, using exactly \(K\)-th-root-of-unity background link phases. The relative covariance residuals on \(2^4\) were below \(3.2\times10^{-16}\) for \(K=5,7\), both \(r=0\) and \(r=1\). On \(3^4\), \(K=5\), residuals were below \(1.5\times10^{-16}\). Determinant phases and log magnitudes were invariant within roundoff.

### 3.1 Explicit non-uniqueness: species doubling

With all links set to unity, the Fourier symbol is

\[
D_r(p)=m+\frac{r}{a}\sum_\mu(1-\cos ap_\mu)
+\frac{i}{a}\sum_\mu\gamma_\mu\sin(ap_\mu).
\]

At \(m=0\):

* **Naive** fermions \((r=0)\): all 16 Brillouin-zone corners \(ap_\mu\in\{0,\pi\}\) are zeros, representing 16 light Dirac tastes.
* **Wilson** fermions \((r=1)\): only the all-zero corner remains light. A corner containing \(k\) components equal to \(\pi\) acquires Wilson mass \(2k/a\), which diverges in the continuum limit.

Both operators use the **same finite gauge-link qudits and exact gauge transformation rules**. Their one-loop coefficients, **conditional on the usual Dirac continuum limit**, are respectively

\[
\beta_{r=0}(e)=16\frac{q^2e^3}{12\pi^2}+\cdots,
\qquad
\beta_{r=1}(e)=\frac{q^2e^3}{12\pi^2}+\cdots.
\]

Thus the finite gauge algebra, charge normalization, and gauge invariance do not determine the coefficient. A physical principle selecting the kinetic operator and the light particle spectrum is indispensable. This does **not** establish that either lattice completion is derived from the original maximally entangled-pair axiom.

## 4. Derive the one-Dirac coefficient from added dynamics

This section uses the continuum limit **of the selected Wilson–Dirac kinetic operator**, together with the standard fermion determinant and heat-kernel asymptotics. It is a conditional calculation, not an inference from the 172 records.

To avoid mixing coupling normalizations, write the unrescaled Abelian gauge field with

\[
D_\mu=\partial_\mu+iq A_\mu,\qquad
S_{\rm gauge}=\frac1{4e^2}\int d^4x\,F_{\mu\nu}F_{\mu\nu}.
\]

The Euclidean Dirac operator is \(\not D=\gamma_\mu D_\mu\). Its positive Laplace-type square can be written

\[
P=-\not D^2=-(D_\mu D_\mu+E),\quad
E=\frac{iq}{4}[\gamma_\mu,\gamma_\nu]F_{\mu\nu},\quad
\Omega_{\mu\nu}=[D_\mu,D_\nu]=iqF_{\mu\nu}I_4.
\]

For flat space and a slowly varying field, the relevant local heat-kernel coefficient **before** the overall \((4\pi t)^{-2}\) factor is

\[
\operatorname{tr}a_4\big|_{F^2}=
\operatorname{tr}\left(\tfrac12E^2+\tfrac1{12}\Omega_{\mu\nu}\Omega_{\mu\nu}\right).
\]

The Clifford trace \(\operatorname{tr}(\gamma_{\mu\nu}\gamma_{\rho\sigma})=-4(\delta_{\mu\rho}\delta_{\nu\sigma}-\delta_{\mu\sigma}\delta_{\nu\rho})\) yields

\[
\begin{aligned}
\operatorname{tr}E^2 &=2q^2F_{\mu\nu}F_{\mu\nu},\\
\operatorname{tr}\Omega_{\mu\nu}\Omega_{\mu\nu}&=-4q^2F_{\mu\nu}F_{\mu\nu}.
\end{aligned}
\]

Therefore

\[
\boxed{\operatorname{tr}a_4\big|_{F^2}=
\left(1-\frac13\right)q^2F_{\mu\nu}F_{\mu\nu}=
\frac23q^2F_{\mu\nu}F_{\mu\nu}.}
\]

The script verifies this exact identity for a single magnetic component \(F_{12}=1\) and a general integer antisymmetric tensor with four independent components; its full Clifford calculations use exact SymPy arithmetic.

For a vector-like Dirac fermion the magnitude of its determinant contributes

\[
W[A]=-\tfrac12\operatorname{Tr}\ln(P+m^2),
\]

and hence its proper-time shell between UV cutoff \(\Lambda\) and renormalization scale \(\mu\) gives

\[
\begin{aligned}
\Delta W\big|_{F^2}
&=\frac12\int_{\Lambda^{-2}}^{\mu^{-2}}\frac{dt}{t}
\frac1{(4\pi t)^2}\int d^4x\;t^2\frac23q^2F^2\\
&=\frac{q^2}{24\pi^2}\ln\frac\Lambda\mu\int d^4x\,F^2.
\end{aligned}
\]

Comparing to \(\frac1{4e^2(\mu)}\int F^2\) determines

\[
\frac1{e^2(\mu)}-\frac1{e^2(\Lambda)}
=\frac{q^2}{6\pi^2}\ln\frac\Lambda\mu.
\]

Differentiating,

\[
-\frac{2}{e^3}\mu\frac{de}{d\mu}=-\frac{q^2}{6\pi^2},
\]

and therefore

\[
\boxed{\beta(e)=\frac{q^2e^3}{12\pi^2},\qquad
\beta(\alpha)=\frac{2q^2}{3\pi}\alpha^2.}
\]

**Interpretation:** No one-loop Feynman-diagram value was inserted as a fit. But the required continuum Dirac operator, dimensionality, determinant, and heat-kernel formalism were *assumed* as dynamical input. The result therefore cannot be labeled a parameter-free prediction of Chain Theory.

## 5. Independent Landau-spectrum check of the coefficient

For constant \(F_{12}=B\), the squared four-dimensional Dirac operator has transverse Landau levels \((2n+1)|qB|\) shifted by spin \(\pm|qB|\), with degeneracy \(|qB|/(2\pi)\) per transverse area. The two longitudinal directions contribute \((4\pi t)^{-1}\). Summing both spin branches gives

\[
\begin{aligned}
\frac{K_B(t)}V&=\frac{|qB|}{2\pi}\frac1{4\pi t}
2\sum_{n=0}^\infty
\left(e^{-2n|qB|t}+e^{-2(n+1)|qB|t}\right)\\
&=\frac4{(4\pi t)^2}\,x\coth x,\qquad x=|qB|t.
\end{aligned}
\]

Without referring to the local \(a_4\) formula, expand the spectrum:

\[
4x\coth x=4+\frac43 x^2-\frac4{45}x^4+\frac8{945}x^6+\cdots.
\]

The \(t^2\) coefficient is \(\frac43q^2B^2\), agreeing with \(\frac23q^2F^2\) since \(F^2=2B^2\). A direct sum over 1000 Landau levels was also checked for \(x\in\{0.4,0.2,0.1,0.05,0.025\}\). The recovered finite-\(x\) ratios approach \(1.333333\ldots\), with the last value \(1.3332777811\). This is an independent spectrum check *within the added continuum Dirac theory*, not within the original qudit postselection model.

## 6. What the (3,2) qudit link contributes to early-universe gauge running

The original paper's \((\mathbf3,\mathbf2)_{Y=1/6}\) representation has

\[
\operatorname{Tr}_{6}Y^2=6(1/6)^2=\boxed{1/6}.
\]

With the conventional electroweak relation \(Q_{\rm em}=T_3+Y\), **not derived by the paper**,

\[
\operatorname{Tr}_{6}Q_{\rm em}^2=3[(2/3)^2+(-1/3)^2]=\boxed{5/3}.
\]

Its Dynkin indices, taking \(T(\mathbf{fund})=1/2\), are

\[
T_{SU(2)}=3(1/2)=3/2,\qquad
T_{SU(3)}=2(1/2)=1.
\]

If **additionally** the representation is realized as one left-handed Weyl multiplet, it contributes the following terms to \(\mu\,dg_i/d\mu=(b_i/16\pi^2)g_i^3\):

| Coupling | Quark-doublet contribution to \(b_i\) |
|---|---:|
| \(g_Y\) (non-GUT normalized) | \((2/3)(1/6)=1/9\) |
| \(g_2\) | \((2/3)(3/2)=1\) |
| \(g_3\) | \((2/3)(1)=2/3\) |

The factor \(2/3\) per Weyl fermion comes from continuum quantum field theory; the representation traces are the only algebraic input connected directly to the paper.

For comparison, supplying the **entire** conventional Standard Model matter roster for one generation,

| Weyl multiplet | multiplicity \(n\) | hypercharge \(Y\) | \(nY^2\) |
|---|---:|---:|---:|
| \(Q_L\) | 6 | \(1/6\) | \(1/6\) |
| \(u_R\) | 3 | \(2/3\) | \(4/3\) |
| \(d_R\) | 3 | \(-1/3\) | \(1/3\) |
| \(L_L\) | 2 | \(-1/2\) | \(1/2\) |
| \(e_R\) | 1 | \(-1\) | \(1\) |
| **Total** | | | **\(10/3\)** |

and introducing three generations, one complex Higgs doublet, and the conventional non-Abelian gauge-field contribution gives

\[
\begin{aligned}
b_Y&=\frac23(3)(10/3)+\frac13(1/2)=\boxed{41/6},\\
b_2&=-\frac{11}{3}(2)+\frac23(3)(2)+\frac13(1/2)=\boxed{-19/6},\\
b_3&=-\frac{11}{3}(3)+\frac23(3)(2)=\boxed{-7}.
\end{aligned}
\]

These exact values verify the accounting and make the early-universe target explicit. **They are not an independent derivation of the Standard Model matter roster or loop coefficients from Chain Theory.** In addition, a chiral \(Q_L\) multiplet alone is anomalous; the vector-like Wilson–Dirac toy is not yet a faithful lattice formulation of the full chiral Standard Model.

## 7. Assumption ledger and decisive research bottleneck

| Element | Origin | Established by original axiom? |
|---|---|---|
| Maximally entangled qudit links, readout, qutrit–qubit gauge group, hypercharge-like generator | Source paper | Within paper's definitions |
| \(\alpha_{\rm cand}=\alpha(0)\) | Photon Conjecture | **No — conjectured** |
| Finite cyclic \(\mathbb Z_K\) gauge-link operator | This extension | No |
| Jordan–Wigner fermion modes and physical CAR dynamics | This extension | No |
| Four-dimensional lattice with spatial locality and physical scale \(a\) | This extension | No |
| Clifford module and Dirac/Wilson kinetic term with parameter \(r\) | This extension | No |
| Fermion determinant, proper-time shell, heat-kernel asymptotics | Continuum QFT inputs | No |
| Electric charge \(Q_{\rm em}=T_3+Y\) | Standard Model electroweak identification | No |
| Three generations, Higgs, and all remaining multiplets | Standard Model field content | No |

The mathematical bottleneck is not the Clifford trace: that is straightforward and exact. It is **selecting the physical Dirac dynamics and fermion multiplicity** from the entangled-qudit definitions. Species doubling provides a sharp counterexample to any assertion that finite-dimensional gauge covariance by itself suffices.

### Suggested falsifiable next test

Construct from the *original* allowed network operations a Euclidean transfer operator \(T\) and its associated long-wavelength inverse correlator. Without stipulating Wilson parameter \(r\), spacetime dimension, fermion anticommutation, or four-dimensional Dirac dynamics, test:

1. Does the transfer operator have gapless charged excitations with a linear, rotationally invariant dispersion?
2. Are those excitations fermionic, and is their low-energy multiplicity one, sixteen, or something else?
3. Is there an intrinsic network coarse-graining operation \(\mathcal C_b\) whose scale \(b\) is operationally definable?
4. Does the gauge-field two-point response exhibit a \(\ln b\) term with coefficient fixed independently of the chosen coarse-graining map?
5. Does that construction retain the source paper's zero-transfer readout \(\alpha_{\rm cand}\) and a single electromagnetic normalization?

**Decision criterion:** Without affirmative, independently computed answers, the result \(2/(3\pi)\) remains a consequence of supplemented continuum physics and does not validate the source paper's Photon Conjecture.

## 8. Reproducibility and sources

Files accompanying this note:

- `chain_beta_probe.py`: existing baseline probe; independently re-executed.
- `chain_beta_extension_v2.py`: explicit tests for exact two-qubit Clifford traces, a four-component field tensor, cyclic finite \(\mathbb Z_K\) gauge invariance, lattice Wilson–Dirac covariance and determinants, the 16-versus-1 doubling result, Landau-spectrum expansion and sums, and Standard Model representation accounting.
- `chain_beta_extension_v2_results.txt`: captured test output.

Execution: `python chain_beta_extension_v2.py` (NumPy and SymPy installed).

Source context:

- A. D. Wissner-Gross, *A Parameter-Free Candidate for the Fine-Structure Constant from Entangled Qudits*, 23 October 2026, uploaded `ChainTheoryAlpha.pdf`, especially main text §§Charge, Theorems 1–2, and Supplemental §§S3, S6, S8.
- *QED in the exact renormalization group*, PTEP (2021), β(e) = e³/(12π²): https://academic.oup.com/ptep/article/2021/12/123B06/6414015
- *Varieties and properties of central-branch Wilson fermions*, Phys. Rev. D 102, 034516 (2020), four-dimensional Wilson species structure: https://doi.org/10.1103/PhysRevD.102.034516
- *On the trace anomaly of a Weyl fermion in a gauge background*, Eur. Phys. J. C (2019), heat-kernel coefficients: https://link.springer.com/article/10.1140/epjc/s10052-019-6799-z
