# Chain Theory → QED beta program, v12: a covariant record instrument and its limits

**Date:** 24 October 2026  
**Primary source:** A. D. Wissner-Gross, *A Parameter-Free Candidate for the Fine-Structure Constant from Entangled Qudits*, October 23, 2026 (`ChainTheoryAlpha.pdf`).  
**Relation to earlier work:** Continues the explicit `31,104 × 172` readout-operator embedding of v11.  
**Status:** A positive construction of an **abstract** trace-preserving instrument on the readout operator Hilbert space; exact bounds and counterexamples to two possible **physical** interpretations. Neither a unique physical measurement nor the QED beta-function coefficient follows from the axiom.

## Executive technical conclusion

Let `R` denote the *operator Hilbert space* specified in the paper, of complex dimension

\[
N=(3^2-1)(2^2-1)(3\cdot2)^2(3\cdot2)^2=24\cdot36^2=31104.
\]

With the same normalized connected unread direction and 172 normalized record directions as in v11, set

\[
M=\sum_{j=1}^{172}w_j|r_j)(r_j|,
\quad w_j\in\{1,\sqrt2,\sqrt3,\sqrt6\}.
\]

The exact multiplicities are `(79,48,42,3)`; the orthogonal complement has dimension `30932`. Therefore

\[
S:=\operatorname{Tr}_{\mathcal R}M
=79+48\sqrt2+42\sqrt3+3\sqrt6
=226.97685414015095\ldots,
\qquad
\alpha_{\rm cand}=S/31104.
\]

**Five results:**

1. `M` has a well-defined spectral (Lüders) instrument on the auxiliary Hilbert space `R`, with five outcomes and Kraus projectors summing to `I_R`. For the *assumed* maximally mixed auxiliary input `tau=I_R/N`, the mean numerical score equals exactly `alpha_cand`.
2. Transporting the record projectors together with a matched local frame transformation yields exactly covariant outcome probabilities and postmeasurement branches. A fixed-frame commutator outcome can instead change from probability one to zero under a transformation that stabilizes the source-inspired three-pair Bell state. This is a coordinate/instrument test; it is not evidence of a physical signal.
3. `M` is **not** a binary POVM effect: some eigenvalues exceed one. In fact *no* effect supported only on the 172 record directions can have probability `alpha_cand` on `tau`, since its trace would be bounded by 172 rather than `S≈226.98`. A binary effect giving the same probability needs rank at least `ceil(S)=227`: at least 55 extra directions if all original 172 directions remain in its support.
4. A finite matched-frame twirl of the three commutator record projectors changes the rank-three projector into a rank-five effect with spectrum `(1, 1/2, 1/2, 1/2, 1/2, 0)` in its six-dimensional closure, while preserving its trace. Therefore gauge-averaged *outcome effects* need not retain the original rank/multiplicity interpretation even when the weighted scalar does.
5. A projective measurement on the **operator-Hilbert space** does not automatically induce a completely positive instrument on the underlying 216-dimensional physical three-pair Hilbert space. The elementary Liouville projector `|sigma_z/sqrt(2))(sigma_z/sqrt(2)|` already fails positivity when interpreted as a map on qubit density operators.

Each statement is supported by algebraic identities below and checked in `chain_beta_extension_v12.py`. The physical meaning of the abstract auxiliary instrument remains to be supplied; the original paper does not assert that `alpha_cand` is a single-outcome probability.

## 1. Source assumptions and a critical typing distinction

The source's Supplemental §S6 defines the readout dimension as the `24`-dimensional traceless–traceless operator sector of an unread qutrit–qubit link, multiplied by the full `36`-dimensional operator spaces of each of two read links. Supplemental §S7 says that the 172 record directions consist of a read-arm record operator tensored with the connected unread-link direction. It assigns one of four amplitude weights to each record, and §S8 *defines* the normalized trace as an analogue of the on-shell QED dressed-to-bare ratio. Supplemental §S11 states that the supplied Lean formalization verifies the record-weight sum and denominator without building the entire readout space itself.

The original normalized joined state is in

\[
\mathcal H_{\rm phys}=(\mathbb C^3\otimes\mathbb C^2)^{\otimes3},
\qquad \dim\mathcal H_{\rm phys}=6^3=216.
\]

Its operator space has dimension `216²=46656`. By contrast, `R` is an **operator-space subspace** of dimension 31104. It is not an independently postulated 31104-level physical system. Equipping `R` with its Hilbert–Schmidt inner product is mathematically natural, but **there is no canonical physical state** `sigma` on `R` specified by the source.

In particular, `tau=I_R/N` below is an auxiliary maximally mixed *density operator on the vector space* `R`, not the paper's original physical state `|Omega><Omega|`. This distinction is necessary to avoid treating an operator-space trace as a Born probability without justification.

## 2. Exact abstract instrument reproducing the weighted trace

Let `P_1, P_sqrt2, P_sqrt3, P_sqrt6` be the orthogonal spectral projectors of `M`, and define

\[
P_0=I_\mathcal R-P_1-P_{\sqrt2}-P_{\sqrt3}-P_{\sqrt6}.
\]

Their ranks are, respectively,

\[
(\operatorname{rank}P_0,\operatorname{rank}P_1,
\operatorname{rank}P_{\sqrt2},\operatorname{rank}P_{\sqrt3},\operatorname{rank}P_{\sqrt6})
=(30932,79,48,42,3).
\]

Define the selective maps on **auxiliary density operators** `sigma ∈ End(R)` by

\[
\mathcal I_w(\sigma)=P_w\sigma P_w.
\]

Each is completely positive, and

\[
\sum_w\operatorname{Tr}\mathcal I_w(\sigma)=\operatorname{Tr}\sigma,
\]

because the projectors sum to the identity. The nonselective map `sigma ↦ Σ_w P_w sigma P_w` is completely positive and trace preserving. This is the spectral/Lüders instrument of the Hermitian observable `M`.

For the auxiliary input `tau=I_R/N`:

\[
p_w=\operatorname{Tr}(P_w\tau)=\operatorname{rank}(P_w)/N.
\]

Numerically the outcomes are

| Score `w` | Projector rank | Probability `p_w` |
|---|---:|---:|
| 0 | 30932 | 0.9944701646090535 |
| 1 | 79 | 0.002539866255144033 |
| √2 | 48 | 0.001543209876543210 |
| √3 | 42 | 0.001350308641975309 |
| √6 | 3 | 0.000096450617283951 |

The expectation is exactly

\[
\mathbb E_\tau[w]=\sum_w wp_w
=\frac{79+48\sqrt2+42\sqrt3+3\sqrt6}{31104}=\alpha_{\rm cand}.
\]

The source's scalar is thus compatible with a legitimate **score expectation** on the auxiliary `R`. Nothing here shows it is a Born probability in `H_phys`, a photon scattering amplitude, or the physical charge renormalization.

### Why the spectral instrument is not selected uniquely

Even on `R`, choosing the spectral Lüders *update rule* adds structure: a POVM alone specifies probabilities, not the disturbance and postmeasurement state. Within each degenerate eigenspace there are many finer-grained instruments giving the same coarse probabilities and numerical expectation. In addition, the source does not specify why `tau`, instead of another density operator on `R`, is the physical input.

## 3. An exact lower bound on any binary-probability interpretation

For a binary POVM effect `E`, positivity and normalization require

\[
0\preceq E\preceq I_\mathcal R.
\]

If `E` is supported only on the selected 172 record directions, `rank(E)≤172`, and every eigenvalue is in `[0,1]`. Therefore

\[
\operatorname{Tr}E\leq172,
\qquad \operatorname{Tr}(E\tau)\leq\frac{172}{31104}.
\]

But `alpha_cand=S/31104` and `S>226`. Hence no such effect exists. More generally, any effect giving the probability `alpha_cand` at `tau` has

\[
\operatorname{rank}(E)\geq\lceil\operatorname{Tr}E\rceil
=\lceil S\rceil=227.
\]

If its support includes all the 172 original record directions, it must include at least 55 additional orthogonal directions. This bound is conditional on the **binary-probability plus maximally mixed auxiliary-input** interpretation. It does **not** constrain the source's trace-as-amplitude interpretation.

### Explicit infinite nonuniqueness

Let `P` be the rank-172 record projector. For any `t in [0,1]`, define

\[
E_t=tP+a_t(I-P),\qquad
 a_t=\frac{S-172t}{N-172}.
\]

For this range `0≤a_t≤1`, so `E_t` is a valid effect, and

\[
\operatorname{Tr}(E_t\tau)=\frac{172t+(N-172)a_t}{N}=\alpha_{\rm cand}.
\]

At `t=alpha_cand`, `E_t=alpha_cand I_R`, a completely trivial state-independent coin flip. At `t=1`, one gets a very different effect: it accepts each record direction with certainty but has nonzero weight `a_1≈0.0017773456013` on the complement. Both have the same probability on `tau`. This explicitly proves that neither the low-energy number nor generic gauge covariance selects a unique binary instrument.

## 4. A gauge-covariant instrument: explicit matched-frame test

Suppose `W` is the orthogonal/unitary representation of a matched qudit-frame change on `R`. Define the transported family

\[
P_w^{(g)}=WP_wW^\dagger,
\qquad \sigma^{(g)}=W\sigma W^\dagger.
\]

The probability covariance and instrument naturality follow identically:

\[
\operatorname{Tr}(P_w^{(g)}\sigma^{(g)})=\operatorname{Tr}(P_w\sigma),
\]

\[
\mathcal I_w^{(g)}(\sigma^{(g)})
=W\mathcal I_w(\sigma)W^\dagger.
\]

This is a **valid constructive covariance law**, but once `W` is supplied it is a consequence of unitary conjugation, not a derivation of `W` from the axiom.

To test the v11 stabilizer issue, take

\[
u_1=\operatorname{diag}(e^{i\pi/4},e^{-i\pi/4}),
\quad u_2=\bar u_1,
\quad G_i=\operatorname{diag}(\bar u_i,1),
\]

acting on the first and second qutrit–qubit read arms, with the third unread arm unchanged. This transformation fixes the source-inspired joined three-pair state. The induced read-arm action `W` rotates a selected commutator-difference direction into the *excluded* symmetric single-qubit direction.

For an auxiliary input supported on that commutator direction, tests give:

- With the original projector: `p_sqrt6=1`.
- Transforming the input while holding the projector fixed: `p_sqrt6=0`, `p_0=1` (up to numerical error).
- Transporting the projectors and input together: `p_sqrt6=1` again, with probability covariance residual `7.8×10^-16`.
- For a generic coherent test vector, all postmeasurement branch vectors obey `P_w^(g)Wpsi=WP_w psi`, with maximum Euclidean residual `6.7×10^-16`.

This validates the mathematical associated-family completion for the candidate v11 embedding. It **does not** prove which transformations of the abstract auxiliary input correspond to the *same physical experiment* when the original `|Omega>` remains unchanged. That identification is still missing.

## 5. A finite phase-stabilizer twirl changes the outcome effect's spectrum

Represent the six single-qubit directions as `q1 ⊕ q2`, with the commutator-difference projector

\[
P_- =\frac12\begin{pmatrix}I_3&-I_3\\-I_3&I_3\end{pmatrix}.
\]

Under the matched phase transformations in §4, the two Pauli-adjoint rotations are `R(+2theta)` and `R(-2theta)` about `z`. An exact average over `theta=j*pi/8`, `j=0,1,2,3`, gives

\[
\boxed{\overline P_-=
\frac12\begin{pmatrix}I_3&-Z\\-Z&I_3\end{pmatrix},
\qquad Z=\operatorname{diag}(0,0,1).}
\]

Its spectrum is

\[
\operatorname{spec}(\overline P_-)=
\{1,\underbrace{1/2,1/2,1/2,1/2}_{4},0\}.
\]

Consequently,

\[
\operatorname{rank}P_-=3,
\quad\operatorname{rank}\overline P_-=5,
\quad\operatorname{Tr}P_-=\operatorname{Tr}\overline P_-=3.
\]

The averaged operator is a legitimate **POVM effect** (spectrum in `[0,1]`), but it is no longer an orthogonal projector. Replacing the commutator contribution `sqrt(6)P_-` by `sqrt(6)\overline P_-` keeps the *weighted trace* `3sqrt(6)` exactly while changing its spectral decomposition. This provides an explicit nontrivial covariant-averaging option and confirms that scalar trace preservation need not preserve the physical outcome interpretation of 172 discrete records.

The finite average is over a designated stabilizer subgroup in this candidate embedding, not a theorem identifying the full gauge quotient of the source.

## 6. Why the auxiliary instrument is not automatically physical

The mistake to avoid is treating rank-one projections in **Liouville space** as physical quantum operations on the underlying system.

The simplest counterexample uses a qubit operator

\[
A=\sigma_z/\sqrt2
=\frac1{\sqrt2}\operatorname{diag}(1,-1).
\]

The rank-one projector `|A)(A|` acts on physical operators `X` through

\[
\mathcal S_A(X)=A\operatorname{tr}(A^\dagger X).
\]

For the valid qubit density matrix `X=|0><0|`,

\[
\boxed{\mathcal S_A(|0\rangle\langle0|)
=\frac12\operatorname{diag}(1,-1).}
\]

That output has a negative eigenvalue, so `S_A` is **not even positive**, much less completely positive. Therefore an orthogonal projector on a space of quantum operators is not automatically a permissible measurement operation on the underlying physical states. A representation of the record score using genuinely physical effects and completely positive selective maps remains to be constructed.

A physical trace-preserving local instrument would automatically obey ordinary no-signaling on averaging over outcomes, by the partial-trace identity

\[
\operatorname{tr}_A[(\mathcal E_A\otimes\mathrm{id}_B)(\rho_{AB})]
=\operatorname{tr}_A\rho_{AB}
\]

for any local completely positive trace-preserving `E_A`. This is a general theorem, not a consequence of the original Chain Theory postselection rule, which explicitly excludes its complement.

## 7. Scope of the result and source consistency

| Claim | Outcome |
|---|---|
| `M` has spectrum 0 plus weights 1, √2, √3, √6 with specified ranks | Reconstructed and verified in explicit v11 readout realization |
| Spectral projections give a complete auxiliary Lüders instrument | **Yes, constructive** |
| Weighted mean at `tau=I_R/31104` equals `alpha_cand` | **Yes, exact** |
| Gauge-covariant transported outcome probabilities and branches | **Yes, exact conjugation; numerically tested** |
| A single Bernoulli event on exactly 172 supported directions can equal `alpha_cand` at `tau` | **No, rank/trace obstruction** |
| Gauge averaging necessarily preserves record ranks/weights | **No, exact twirl counterexample** |
| A Lüders instrument on `R` is automatically a CP instrument on the original three-pair physical space | **No, elementary positivity counterexample** |
| Source axiom picks `tau`, the Kraus update, a unique support transport, or a unique measurement instrument | **Not shown** |
| Derivation of photon propagation, matter spectrum, or QED running coefficient | **Not achieved** |

The correct interpretation of v12 is a **mathematical existence result with a physical implementation gap**. It supplies a consistent abstract measurement of the weighted readout, but does not establish that `alpha_cand` is the electromagnetic coupling measured by a realistic detector. It does not contradict the source's conditional number.

## 8. Immediate next targets

1. Construct a **physical** POVM `{F_w}` and trace-preserving Kraus instrument on `H_phys`, with a source-derived map from the joined three-pair state to actual probabilities. Require compatibility with the paper's chosen `P3` branch and the pair-only ontology; do not assume an auxiliary `R` state is physically preparable.
2. Prove a **gauge-descent/naturality theorem** for that physical instrument under the state stabilizer and any independently allowed matching transformations, explicitly specifying which outcomes transform into one another.
3. If the original five score sectors cannot be realized, identify the minimal additional operational axiom and whether its unconditioned readout still equals the same `alpha_cand`. Keep frame covariance, branch selection, and physical no-signaling separate.
4. Only then pursue emergent propagation, mass scales, and the QED beta function. No progress on a physical momentum-dependent vacuum polarization follows solely from the present finite-dimensional score operator.

## 9. Reproduce

Requirements: Python 3, NumPy, SymPy, `chain_beta_extension_v9.py`, `chain_beta_extension_v11.py` (and that script's import dependencies). All included in the bundle.

```bash
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 python chain_beta_extension_v12.py
```

Results are emitted as `chain_beta_extension_v12_results.json`. The program checks the full four-sector decomposition using a compressed 175-dimensional operator-coordinate representation that extends analytically to the `31104`-dimensional readout. It uses exact SymPy arithmetic for the commutator subgroup twirl and for the physical positivity counterexample; matched-frame covariance is tested numerically to machine precision. It does not allocate a `31104 × 31104` matrix.