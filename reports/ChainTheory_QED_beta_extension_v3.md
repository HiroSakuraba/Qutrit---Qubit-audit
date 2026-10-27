# Chain Theory, fermion dynamics, and uniqueness of the QED one-loop coefficient

**Research continuation, v3 — 24 October 2026**  
**Target:** Derive the first QED beta-function coefficient from the source qudit construction independently of ordinary QED.  
**Outcome:** An explicit qubit-walk *candidate* for emergent Dirac propagation was constructed, exact overlap-fermion covariance and Ginsparg–Wilson identities were tested, and two independent non-uniqueness obstructions were established. No derivation from the original Chain Theory axiom alone was obtained.

## 0. Epistemic legend

- **[SOURCE]** A result explicitly in A. D. Wissner-Gross, *A Parameter-Free Candidate for the Fine-Structure Constant from Entangled Qudits*, dated 23 October 2026 (the user-supplied 24-page PDF).
- **[DERIVED]** A deduction following algebraically from a clearly stated set of inputs.
- **[ADDED]** Hypothesis or structure not contained in the source qudit axiom.
- **[TESTED]** Exact symbolic or floating-point numerical check in `chain_beta_extension_v3.py`.
- **[CONDITIONAL]** Result relying on established continuum quantum-field-theoretic input.

The distinction is important: reproducing a one-loop factor after adding a Dirac determinant is **not** an independent explanation of why Chain Theory has that factor.

## 1. Precise starting point and goal

**[SOURCE]** The source postulates finitely many maximally entangled finite-dimensional qudit pairs, sharing and distillation. Among readable pairs, a qutrit–qubit link carries an acting central charge and a gauge *algebra* isomorphic to `su(3) ⊕ su(2) ⊕ u(1)`; the global group given is `S(U(3)×U(2))`. The model defines a normalized two-point readout of 172 weighted records, with `dim(R)=31104`, giving a parameter-free candidate for `alpha(0)`. The **Photon Conjecture** identifies that scalar with zero-momentum electromagnetic coupling.

**[SOURCE]** Supplement §S3 says the tensor-product link `C^3⊗C^2` transforms as `(3,2)_{1/6}`, with hypercharge `Y=1/6`. Supplement §S6 normalizes the conserved generator's trace; §S8 explicitly states that it does not construct a momentum-dependent propagator, charged loop, loop momentum, or vertex. Its `Z1=Z2` assertion is an *analogy* of roles; it does not establish a momentum-space Ward–Takahashi identity for physical correlation functions.

**Target** (one four-dimensional Dirac fermion of unit electric charge):

\[
\boxed{\mu\,\frac{de}{d\mu}=\frac{e^3}{12\pi^2}+O(e^5),\qquad \mu\frac{d\alpha}{d\mu}=\frac{2}{3\pi}\alpha^2+O(\alpha^3).}
\]

The original `alpha(0)` count fixes a boundary *number*, not the scale derivative.

## 2. Finite-dimensional Choi ambiguity: no preferred kinetic update

**[DERIVED, TESTED]** Let `|Ω_d⟩=1/√d Σ_j |j⟩|j⟩`. For any unitary `U` on `C^d`,

\[
|J(U)\rangle=(I\otimes U)|\Omega_d\rangle,
\quad
\rho_L=\rho_R=I_d/d.
\]

All these pure states have identical maximal bipartite entanglement (entropy `log₂ d`, purity `1/d`). Consequently, maximal entanglement by itself does not select a unitary update or a kinetic generator. Additional relational constraints elsewhere in a network *might* distinguish choices, but the lone pair and its entanglement invariant do not.

**[TESTED]** For `d=4`, choose `U_0=I_4` and `U_1=diag(1,1,1,-1)`. Both Choi states have reduced purity `1/4` and entanglement entropy 2 bits, despite mutual state-overlap magnitude `|⟨J(U_0)|J(U_1)⟩|=1/2`. Hence they are distinct states with the same entanglement statistics. This tests a *kinematic underdetermination*, not a model-equivalence theorem under all allowed full-network reads.

This is the first missing selection problem: neither a discrete clock nor a preferred propagation operator follows from equal Schmidt coefficients.

## 3. Constructive toy model: one qubit per site can exhibit emergent Dirac propagation

**[ADDED]** Introduce a one-dimensional translation-invariant lattice, discrete update time, conditional shift and onsite SU(2) coin. In momentum space, with qubit Pauli matrices,

\[
U(k,m)=e^{-im\sigma_x}e^{-ik\sigma_z}.
\]

Its unitarity is exact. One step of the conditional shift moves the two qubit components in opposite lattice directions; the coin rotates them on site. The exact quasienergy relation is

\[
\tfrac12\operatorname{tr}U=\cos m\cos k=\cos\omega(k,m).
\]

Expanding simultaneously at small `k` and `m`,

\[
U(k,m)=I-i(k\sigma_z+m\sigma_x)+O(k^2+m^2+|km|),
\quad
\omega(k,m)=\sqrt{k^2+m^2}+O((|k|+|m|)^3).
\]

Thus the effective long-wavelength Hamiltonian is the **1+1 dimensional Dirac Hamiltonian**

\[
H_{\rm eff}(k)=k\sigma_z+m\sigma_x.
\]

**[TESTED]** With `(m,k)=(0.7,0.9)ε` for `ε=0.2,0.1,0.05,0.025`, the absolute quasienergy errors relative to `sqrt(k²+m²)` were respectively `4.6623×10^-4`, `5.8082×10^-5`, `7.2542×10^-6`, `9.0659×10^-7`. Each halving reduces the error by ~8, consistent with cubic order.

**Limitations:** This is a one-dimensional *free* toy construction. It introduces the graph, translation, time step, coin, and locality. The source axiom does not select this walk, its parameters, or the corresponding 3+1-dimensional operator. This is evidence of a viable *engineering path*, not of unique dynamics or of the QED beta function.

## 4. Stronger non-uniqueness even with exact lattice chirality

**[ADDED]** Work in a four-dimensional Euclidean lattice, spacing `a=1`, with two-qubit representations of four gamma matrices. A Wilson operator has free momentum symbol

\[
D_W(p;m)=m+\sum_{\mu=1}^4(1-\cos p_\mu)
+i\sum_{\mu=1}^4\gamma_\mu\sin p_\mu.
\]

At the 16 Brillouin-zone corners, write `k` for the number of components equal to `π`. The scalar mass is `m+2k`; the kinetic sine terms vanish. Tuning `m=-2k` produces exactly `C(4,k)` massless Dirac corners:

| Wilson tuned mass (`a=1, r=1`) | `k` | Massless corners | Orientation parity `(-1)^k` |
|:--|--:|--:|--:|
| `0` | 0 | 1 | + |
| `-2` | 1 | 4 | − |
| `-4` | 2 | 6 | + |
| `-6` | 3 | 4 | − |
| `-8` | 4 | 1 | + |

**[TESTED]** These results are verified by exact enumeration of all 16 corners. They apply at the free critical points, not generically at arbitrary masses or gauge couplings.

To strengthen the symmetry requirements, define the massless Neuberger overlap operator

\[
H_W(m_0)=\gamma_5(D_W(0)-m_0 I),
\quad D_{\rm ov}(m_0)=I+\gamma_5\operatorname{sgn}H_W(m_0).
\]

When `H_W` has no zero eigenvalue, `V=γ₅ sign(H_W)` is unitary and satisfies `V†=γ₅Vγ₅`; as a direct algebraic consequence,

\[
\boxed{\gamma_5D_{\rm ov}+D_{\rm ov}\gamma_5=D_{\rm ov}\gamma_5D_{\rm ov}.}
\]

This is the Ginsparg–Wilson relation in normalization `a=1` for `D=1+V`. (Other literature conventions include an overall factor `m₀`, changing the right-hand normalization accordingly.)

At a free momentum corner with `k` components equal to `π`,

\[
D_{\rm ov}(p_{\rm corner})=1+\operatorname{sgn}(2k-m_0),
\]

which is zero precisely when `2k<m₀`, provided `m₀≠2k`. Thus the **exact free species count** is

\[
\boxed{N_{\rm light}(m_0)=\sum_{k=0}^4 \binom{4}{k}\mathbf1\{2k<m_0\}.}
\]

| Overlap mass interval | Free massless Dirac corners |
|:--|--:|
| `m₀<0` | 0 |
| `0<m₀<2` | 1 |
| `2<m₀<4` | 5 |
| `4<m₀<6` | 11 |
| `6<m₀<8` | 15 |
| `m₀>8` | 16 |

The threshold values `m₀=0,2,4,6,8` make `H_W` singular at certain free corners; the sign operator is not well defined there without extra prescription. **[TESTED]** All six open-interval counts were recovered by matrix evaluation at all corners. With generic momentum `(0.14,0.27,-0.32,0.55)`, the Ginsparg–Wilson defect is `2.95×10^-16`, and the overlap eigenvalues lie on `|λ−1|=1` to `2.22×10^-16`.

### 4.1 Gauge-background test, and its scope

**[ADDED, TESTED]** On a periodic `2^4=16`-site 4D lattice, each oriented link is assigned a classical background phase in `Z_K` for `K=5,7`. Test the Wilson operator's gauge covariance, form the full Hermitian `H_W`, evaluate its matrix sign by spectral decomposition, and check the overlap operator for covariance, `γ₅`-Hermiticity, the Ginsparg–Wilson identity, and its eigenvalue circle. For both K, the largest relative defect of covariance/GW/Hermiticity was below `1.5×10^-15`; the eigenvalue-circle deviation was below `1.0×10^-14`.

**Important distinction:** These 4D tests use *classical discrete background link phases*, not a quantum gauge-field register for every link. The separate two-site v2 test does quantize a finite `Z_K` link. Neither one tests the physical continuum limit as `K→∞`, volume `→∞`, and lattice spacing `a→0`; gauge coupling renormalization was not numerically extracted from these 16 sites.

**Consequence:** Requiring finite qudits, gauge covariance, and even this exact lattice chiral identity **still** does not independently select one light species. The common choice `0<m₀<2` yields one, but a source-derived rule is needed to choose that window and the overlap construction itself.

This complements the Nielsen–Ninomiya obstruction: naive local lattice fermions cannot maintain the entire desired set of chirality/locality/translation properties without complications, and Wilson or overlap formulations make additional choices.

## 5. The source quark-doublet representation is not by itself a consistent chiral gauge theory

**[SOURCE]** The qudit link has the gauge quantum numbers `(3,2)_{1/6}` as a *representation*, but the paper does not derive fermionic chirality or dynamical Weyl fields.

**[CONDITIONAL]** If one promotes exactly this representation to a *single left-handed Weyl quark doublet* `Q_L`, the SU(2) global-anomaly parity counts **three** fundamental left-handed weak doublets (one per color), which is odd. Witten's global SU(2) anomaly rules out such a chiral spectrum standing alone. The mixed anomaly is

\[
\mathcal A_{SU(2)^2Y}(Q_L)=3\left(\tfrac16\right)T(\mathbf2)
=3\left(\tfrac16\right)\left(\tfrac12\right)=\boxed{\tfrac14}.
\]

A lepton doublet `L_L=(1,2)_{-1/2}` contributes `-1/4`, yielding an even total number of doublets `3+1=4`, but **other perturbative anomaly constraints remain** until the rest of the fermion spectrum is supplied.

A conventional anomaly-free generation, written entirely using left-handed Weyl fields and **without** a sterile neutrino, is

\[
\boxed{Q_L:(3,2)_{1/6},\quad u^c:(\bar3,1)_{-2/3},\quad d^c:(\bar3,1)_{1/3},\quad L_L:(1,2)_{-1/2},\quad e^c:(1,1)_{1}.}
\]

**[TESTED]** For this *imported* matter roster, the gravitational–hypercharge sum, `U(1)_Y^3`, `SU(2)^2U(1)_Y`, `SU(3)^2U(1)_Y`, and `SU(3)^3` anomaly sums vanish exactly (SymPy rationals), and the SU(2) doublet count is even. This is not a prediction of the qudit model.

### 5.1 Even anomaly cancellation cannot uniquely fix beta-function coefficients

There is another rigorous non-uniqueness test conditional only on standard four-dimensional anomaly and beta-function bookkeeping. Add a vectorlike `U(1)`-charged Dirac pair. In the left-Weyl convention the two fields carry charges `+q` and `−q`. Their cubic and mixed gravitational anomalies cancel:

\[
q^3+(-q)^3=0,\qquad q+(-q)=0,
\]

while their positive charge-squared trace is

\[
q^2+(-q)^2=2q^2.
\]

For a unit-charge vectorlike pair, the first-loop Abelian beta coefficient changes by

\[
\Delta b_{U(1)}=\tfrac23\times2=\boxed{\tfrac43},\qquad
\Delta\beta(e)=\frac{e^3}{16\pi^2}\frac43=\boxed{\frac{e^3}{12\pi^2}}.
\]

This contribution applies at scales above the pair's mass (or in a mass-independent scheme before threshold matching); below a heavy pair's mass it decouples from the physical low-energy running. **[TESTED]** All anomaly sums and the positive charge-squared shift were checked with exact rational arithmetic.

Thus not even gauge invariance **plus anomaly cancellation** determines the magnitude of the beta function without a rule selecting the physical particle content and masses. This is independent of the overlap-parameter ambiguity.

## 6. What has and has not been derived

| Statement | Status |
|---|---|
| Source qutrit–qubit readout and candidate `alpha(0)` | **[SOURCE]**; source computation |
| Link representation `(3,2)_{1/6}` | **[SOURCE]** |
| Electric charge `T₃+Y` and the charge-square trace `5/3` | **[CONDITIONAL]** electroweak identification; trace derived |
| Finite-qudit Bell/Choi states allow multiple unitary evolutions | **[DERIVED, TESTED]** no unique update from entanglement |
| A local 1D qubit walk approaches a Dirac Hamiltonian | **[ADDED, TESTED]** requires translation/time/coin |
| Four-dimensional overlap GW identity and gauge covariance | **[ADDED, TESTED]** |
| Unique light fermion species under those requirements | **Disproved for this extension** by explicit examples |
| Viability of lone chiral `(3,2)_{1/6}` | **[CONDITIONAL]** anomalous without additional fields |
| Unique beta coefficient from anomaly cancellation | **Not implied**; vectorlike additions change it |
| Universal one-loop `4/3` Dirac factor from an assumed 4D Dirac determinant | **[CONDITIONAL]** previous v2 result |
| Universal one-loop factor from *only* original Chain Theory axiom | **Not established** |

## 7. Proposed next proof program

A productive route is to formalize **which extra constraints are genuinely required**, then see whether any can be derived from shared-qudit distillation instead of postulated:

1. **Operational causality and graph geometry:** construct distinguishable spacetime neighborhoods from the permitted read/join operations. Define a microscopic time-advance map and prove a meaningful locality bound for it.
2. **Quantum cellular automata:** classify finite-range translation-invariant unitary updates on minimal qudit networks under stated isotropy and symmetry assumptions. Recover Weyl/Dirac dispersion near appropriate modes; track *all* Brillouin/Floquet branches, not only a chosen expansion point.
3. **Fermion statistics:** prove a spin/statistics or parity-superselection principle; do not identify qudit tensor factors with Weyl particles merely because the gauge representation matches.
4. **Chiral matter and anomalies:** attempt a complete, anomaly-free multiplet derivation, including why the other Standard Model representations occur and why vectorlike additions are absent or massive.
5. **Scale-dependent two-point correlator:** construct a physical gauge-current response in momentum space and show a transverse Ward–Takahashi identity. Only then compute the coefficient of the logarithmic term, checking universality under allowed discretizations and renormalization prescriptions.
6. **Continuum/infrared regime:** state whether a limit of finite systems is admissible under the original axiom, and distinguish an exact continuum theory from a controlled finite-scale approximation. The early universe requires electroweak `g_2(μ), g_Y(μ)` rather than an unbroken electromagnetic coupling above the crossover.

**Falsification criterion:** If two extensions satisfying a fully specified set of *source-derived* constraints lead to different continuum spectra or logarithmic coefficients, those constraints are not sufficient. Conversely, merely choosing a model with the known continuum value is not a parameter-free derivation.

## 8. Reproduction

Requires Python 3, NumPy and SymPy. Run:

```bash
python chain_beta_extension_v2.py
python chain_beta_extension_v3.py
```

The v3 script performs floating-point two-qubit Clifford checks; Choi orbit; Wilson-corner enumeration; overlap corner-count and smooth-momentum tests; 4D finite `Z_5` and `Z_7` **classical background** covariance and Ginsparg–Wilson tests; one-dimensional quantum walk continuum scaling; exact chiral-anomaly checks; and vectorlike beta-shift arithmetic.

**Audit correction to v2:** The checked-in v2 script directly runs 4D `2^4`-site background tests. Its narrative also mentioned `3^4`-site tests, but these are not reproducible from that script as provided; do not treat the extra lattice size as verified by the bundled code. The v3 verification deliberately reports only its actual `2^4`-site cases.

## 9. References

Source: A. D. Wissner-Gross, *A Parameter-Free Candidate for the Fine-Structure Constant from Entangled Qudits*, 23 October 2026, supplied PDF, especially main §§“The charge”, “The object”, “Why this is the on-shell coupling”; Supplemental §§S3, S6–S8.

- H. B. Nielsen and M. Ninomiya, “A no-go theorem for regularizing chiral fermions,” *Physics Letters B* **105**, 219–223 (1981). https://doi.org/10.1016/0370-2693(81)91026-1
- E. Witten, “An SU(2) anomaly,” *Physics Letters B* **117**, 324–328 (1982). https://doi.org/10.1016/0370-2693(82)90728-6
- Overlap and Ginsparg–Wilson operator references, including an explicit convention for `m0`: https://academic.oup.com/ptep/article/2019/11/113B06/5645114
- G. M. D'Ariano and P. Perinotti, “Derivation of the Dirac equation from principles of information processing,” *Physical Review A* **90**, 062106 (2014). https://doi.org/10.1103/PhysRevA.90.062106
- P. Raynal, “Simple derivation of the Weyl and Dirac quantum cellular automata,” *Physical Review A* **95**, 062344 (2017). https://doi.org/10.1103/PhysRevA.95.062344

These cited quantum-cellular-automata constructions **introduce additional principles** (e.g. locality, unitarity, homogeneity, discrete isotropy); they are useful templates, not proofs that the source qudit axiom already implies those principles.