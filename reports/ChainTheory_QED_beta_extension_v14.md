# Chain Theory → QED beta program, v14
## A canonical CP map from record directions to physical effects: construction, exact probabilities, and nonuniqueness

**Date:** 25 October 2026  
**Primary source:** A. D. Wissner-Gross, *A Parameter-Free Candidate for the Fine-Structure Constant from Entangled Qudits*, 23 October 2026, `ChainTheoryAlpha.pdf` (especially Letter §§The count/Theorem 2 and Supplemental §§S6–S8).  
**Dependencies:** v9 (read-arm operator basis), v11 (172 record embedding), v12 (auxiliary POVM), v13 (physical Bell-join instrument).  
**Status:** New constructive theorem and exact source-state countercheck. No independent derivation of the QED beta function.

### Abstract

The source identifies the electromagnetic zero-transfer coupling *conjecturally* with a normalized trace of a positive weighted operator on a 31,104-dimensional Hilbert–Schmidt operator space. Previous versions built measurement instruments directly on that abstract space but did not identify a physical observable on the original three qutrit–qubit links.

Here we use the actual embedding of the full record space as a *linear subspace* of operators on the 216-dimensional three-link Hilbert space to define a canonical (relative to that embedding) **unital completely positive map** from matrices acting on the record space to physical measurement effects. For this embedding, the frame identity is exactly `sum_a A_a† A_a = 144 I_216`. The resulting five-outcome POVM has the paper's weight sectors, and its weighted score agrees with the proposed alpha in a maximally mixed physical state. Crucially, **it does not agree in either source pure state**: the initial tensor product of three maximally entangled links yields `0.033023429261197407...`; the normalized postselected three-pair state yields `0.034289813057854211...`. These are exact radical expressions, verified independently by SymPy with rational coefficients.

The CP map is gauge covariant but **not uniquely selected**: a convex family of distinct unital, covariant CP maps is available. It is also nonlocal with respect to the qutrit-join Alice/Bob split: its nonzero effects act on the unread link, including the remote qubit. Accordingly this is a mathematically valid *global* physical POVM, not a derivation of the paper's local distillation operation, and it does not reproduce the electromagnetic coupling on the source states.

## 1. Distinguishing three spaces

The source's physical pre-join state consists of three qutrit–qubit maximally entangled links:

\[
\mathcal H_{\rm phys}=\bigl(\mathbb C^3\otimes\mathbb C^2\bigr)^{\otimes 3},
\qquad \dim\mathcal H_{\rm phys}=6^3=216.
\]

The **readout operator-Hilbert space** is the orthogonal subspace

\[
\mathcal R=\operatorname{End}(\mathbb C^6)\otimes
\operatorname{End}(\mathbb C^6)\otimes
\bigl(\mathfrak{su}(3)_{\rm Herm}\otimes\mathfrak{su}(2)_{\rm Herm}\bigr),
\qquad \dim\mathcal R=36\cdot36\cdot24=31104.
\]

Here `su(d)_Herm` denotes traceless Hermitian `d x d` matrices, regarded as a real subspace; its complex linear span is used to define the operator-Hilbert space. The **record-weight operator** `M` acts on `R`, not directly on `H_phys`. It has nonzero eigenvalues `1, sqrt(2), sqrt(3), sqrt(6)` with respective multiplicities `79, 48, 42, 3`. The remaining 30,932 eigenvalues are zero.

An important correction in terminology: despite the paper's phrase “readout algebra,” the 24-dimensional connected unread-link factor is **not closed under operator multiplication**. For the CP construction we use only its Hilbert–Schmidt inner product and its embedding in `End(H_phys)`, not a claimed C*-algebra structure.

The v11 representative record eigenoperators have the factorized form

\[
A_r = R_r\otimes C_*, \quad
\langle R_r,R_s\rangle_{\rm HS}=\delta_{rs},\quad
\operatorname{tr}(C_*^\dagger C_*)=1,
\quad r=1,\dots,172,
\]

with `R_r` in the two read arms' operator space and `C_*` the normalized connected part of the unread pair density operator. This representation is source-motivated but not established to be the unique physical realization of the paper's 172 records.

## 2. Unital completely positive Heisenberg bridge

Fix any Hilbert–Schmidt-orthonormal basis `{A_a}` of `R` as embedded operators on `H_phys`. The relevant quadratic identities for normalized complete/traceless bases are

\[
\begin{aligned}
\sum_{a=1}^{36} B_a^\dagger B_a &=6I_6 &&(\operatorname{End}(\mathbb C^6)),\\
\sum_{a=1}^{8}L_a^\dagger L_a &=\frac83 I_3 &&(\text{traceless qutrit operators}),\\
\sum_{a=1}^{3}P_a^\dagger P_a &=\frac32 I_2 &&(\text{traceless qubit operators}).
\end{aligned}
\]

Factorization therefore gives the exact tight operator-frame identity

\[
\boxed{\sum_{a=1}^{31104} A_a^\dagger A_a
=6\cdot6\cdot\frac83\cdot\frac32 I_{216}
=144I_{216}.}
\]

Define the isometry

\[
V:\mathcal H_{\rm phys}\to\mathcal R\otimes\mathcal H_{\rm phys},
\quad V|\psi\rangle=\frac1{12}\sum_a|a\rangle\otimes A_a|\psi\rangle.
\]

Then `V† V=I_phys`. For `X` an operator on the record Hilbert space, define

\[
\boxed{\Phi(X)=V^\dagger(X\otimes I_{216})V.}
\]

This is **unital and completely positive** by Stinespring compression. For a record projector `|r)(r|` with a real-coefficient Hermitian representative `A_r` in the chosen Hermitian basis, one gets

\[
\Phi(|r)(r|)=\frac{A_r^\dagger A_r}{144}\succeq0.
\]

A different complex basis requires attention to coefficient conjugation in this formula; the definition via `V` is the unambiguous, basis-covariant construction. The paper's record operators have real coefficients in the adopted orthonormal Hermitian basis.

Let `P_w` denote the record-space spectral projectors of `M`, for `w` in `{0,1,sqrt2,sqrt3,sqrt6}`. Define the physical effects

\[
E_w=\Phi(P_w),\qquad E_w\succeq0,\qquad\sum_w E_w=I_{216}.
\]

The Lüders operations `rho -> sqrt(E_w) rho sqrt(E_w)` constitute a trace-preserving quantum instrument on the physical 216-dimensional system. This is a **global** physical instrument; it has not been shown local to Alice's qutrit join.

A useful additional identity follows by tracing the Stinespring compression:

\[
\operatorname{tr}_{216} \Phi(X)=\frac{\operatorname{tr}_{31104}X}{144}
\quad \text{for every }X,
\]

because the embedding basis is Hilbert–Schmidt orthonormal. Thus in the maximally mixed *physical* state,

\[
\operatorname{tr}\!\left(\frac{I_{216}}{216}\Phi(M)\right)
=\frac{\operatorname{tr}_{\mathcal R} M}{216\cdot144}
=\frac{79+48\sqrt2+42\sqrt3+3\sqrt6}{31104}.
\]

This is a **trace identity**. It does not by itself supply a physical reason for the system to be in the maximally mixed state.

## 3. Exact evaluation on the source's states

Define the physical input state

\[
|\Psi_{\rm in}\rangle=
|\omega_2\rangle_{Q_1q_1}
|\omega_2\rangle_{Q_2q_2}
|\omega_2\rangle_{Q_3f}.
\]

After the paper's selected `P3` join the normalized state is

\[
|\Omega\rangle
=|\omega_3\rangle_{Q_1Q_2}
|\omega_2\rangle_{q_1q_2}
|\omega_2\rangle_{Q_3f}.
\]

The selected join has probability `1/6` on the input, independently verified as a 216-component vector identity. For both these source states, the connected unread factor obeys

\[
\langle\omega_2|C_*^\dagger C_*|\omega_2\rangle=\frac34.
\]

The physical effects factor as `E_w = S_w tensor C_*†C_*/144`, where

\[
S_w=\sum_{r: w_r=w}R_r^\dagger R_r.
\]

Independent exact arithmetic gives:

| Weight | Initial unjoined state `p_w` | Postselected joined state `p_w` | Maximally mixed `p_w` |
|---|---:|---:|---:|
| 1 | 103/9216 | 55/6912 | 79/31104 |
| sqrt(2) | 1/144 | 1/108 | 48/31104 |
| sqrt(3) | 175/27648 | 25/3456 | 42/31104 |
| sqrt(6) | 1/2304 | 1/3456 | 3/31104 |
| 0 | 1685/1728 | 749/768 | 30932/31104 |

The corresponding exact mean scores are

\[
\boxed{\begin{aligned}
\alpha_{\rm initial}^{\Phi}
&=\frac{103}{9216}+\frac{\sqrt2}{144}
+\frac{175\sqrt3}{27648}+\frac{\sqrt6}{2304}\\
&=0.033023429261197407\dots,\\
\alpha_{\rm joined}^{\Phi}
&=\frac{55}{6912}+\frac{\sqrt2}{108}
+\frac{25\sqrt3}{3456}+\frac{\sqrt6}{3456}\\
&=0.034289813057854211\dots.
\end{aligned}}
\]

By comparison,

\[
\alpha_{\rm cand}=0.00729735256366226\dots.
\]

Thus **this physical realization does not reproduce the source's scalar on the source's own pure states**. No beta-function statement follows from this calculation.

The exact results are generated by `chain_beta_extension_v14_exact.py`, which uses SymPy's exact matrices and rational radicals for all 172 two-read record directions. The independent NumPy implementation reproduces the results and checks eigenvalue positivity.

## 4. Covariance and nonuniqueness

The space `R` is invariant under independent local unitary adjoint actions on the physical qudit factors. If the *state, record-space projector, and physical effect* are transported consistently, the bridge is covariant. Twelve seeded random `U(3) x U(2)` transformations on each of the three arms confirmed probability covariance at residual at most `1.11e-16`.

But covariance and complete positivity do **not** select this bridge uniquely. Consider the trace map

\[
\mathcal T(X)=\frac{\operatorname{tr}_{\mathcal R} X}{31104}I_{216}
\]

and the one-parameter family

\[
\boxed{\Phi_t=t\Phi+(1-t)\mathcal T,\qquad0\le t\le1.}
\]

Every `Phi_t` is unital, completely positive and gauge-covariant. The associated score on the initial state is

\[
\langle M\rangle_{\Phi_t,\Psi_{\rm in}}
=t(0.033023429261\dots)+(1-t)(0.007297352564\dots).
\]

At `t=0`, the original number reappears, but *only because the trace map erases dependence on the physical state*. This is not a derivation of the physical coupling. At `t=1`, the direct operator-frame map has a different score. Any intermediate value can occur without violating the tested mathematical constraints.

The construction is thus **nonunique even within the family of unital, gauge-covariant CP maps**.

## 5. Locality obstruction

The physical effects above cannot be identified with the original *local* Bell distillation at the joining qutrits. The unread factor is `C_*†C_*`, which does not act trivially on the unread remote qubit `f`.

In the six-dimensional unread-link space, with Pauli `X_f` on the qubit,

\[
\boxed{\bigl\|[C_*^\dagger C_*,I_3\otimes X_f]\bigr\|_{\rm HS}^2=\frac89.}
\]

This equality was verified in exact SymPy arithmetic. Therefore any nonzero sector effect `E_w` with `S_w !=0` fails to commute with an observable on Bob's remote unread qubit. A POVM implemented solely at the joining qutrits would have effects of the form `E_A tensor I_remote` and *would* commute with every remote operator.

This bridge does not automatically preserve the paper's selected `P3` as one of its outcomes, nor does it give Alice a permitted local implementation. It supplies a valid **global** measurement on the three links and hence is best read as a mathematical baseline for any proposed local completion.

## 6. Implications for the original fine-structure candidate and QED beta function

**What is proved:**

1. A unital completely positive operator-frame bridge exists from the source-motivated 31,104-dimensional record space to physical effects on 216 dimensions.
2. The bridge is gauge covariant when the state and effects are jointly transported.
3. The maximally mixed physical state returns the source's scalar exactly, as a trace identity.
4. The source's actual initial and postselected pure states give *different* exact weighted means under this bridge.
5. The mathematical criteria of positivity, unitality, and covariance admit a continuous family of different bridges.
6. This bridge is global and is not a local replacement for the paper's qutrit Bell join.

**What is not proved:**

- That the paper specifies `Phi`, any `Phi_t`, or a physically realized instrument.
- That the source's 172 records are a unique physical POVM.
- That the paper's scalar is an observable electromagnetic coupling, beyond its explicitly conjectural photon identification.
- That any scale-dependent coupling, propagator, fermionic spectrum, or QED beta function emerges from the axiom.

**Next high-value task (v15):** Solve or constrain the **local, matched-frame, Bell-branch-preserving** instrument problem. Explicitly require effects/Kraus operations on Alice's `Q1Q2` alone, preserve `K0=P3`, ensure frame transport without an externally chosen alignment, and seek any physical score related nontrivially to the 172 record sectors. The v14 operator-frame bridge can serve as a deliberately nonlocal baseline. If no such local bridge exists under the minimal axioms, isolate the weakest additional assumption.

## Reproducibility

- `chain_beta_extension_v14.py` — complete positivity/unitality proof checks, source states, sector probabilities, eigenvalue positivity, 12 seeded covariance tests, nonuniqueness family, exact unread locality witness.
- `chain_beta_extension_v14_exact.py` — independent **exact** SymPy computation of the 172 read-arm record directions on both source states, yielding the rational sector probabilities and radical score expressions.
- `chain_beta_extension_v14_results.json` and `chain_beta_extension_v14_exact_results.json` — machine-readable outputs.
- Prior v2–v13 files retained in the research archive. Numerical and exact checks test this explicit realization; they cannot prove the source axioms uniquely imply it.