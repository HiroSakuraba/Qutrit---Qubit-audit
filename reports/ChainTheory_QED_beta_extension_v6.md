# Chain Theory and the QED beta function — v6
## Trace-preserving operational completions and the fate of the 172 records

**Research date:** 24 October 2026  
**Status:** Exact finite-dimensional consistency construction and nonuniqueness result. **No independent QED beta-function derivation.**  
**Primary source:** A. D. Wissner-Gross, *A Parameter-Free Candidate for the Fine-Structure Constant from Entangled Qudits*, 23 October 2026, supplied as `ChainTheoryAlpha.pdf`, especially the main text and Supplemental Material §§S1, S7, S8.  
**Reproducibility:** `chain_beta_extension_v6.py`; captured output `chain_beta_extension_v6_results.json`.

## 0. Executive result

There are four precise conclusions.

1. The source's selected qutrit Bell projection can be included as **one outcome** of a standard completely positive, trace-preserving (CPTP) instrument. With all outcomes retained and not communicated, the distant two-qubit marginal is unchanged. With only the selected branch deterministically enforced, a local intervention violates no-signalling *if* standard independent remote measurements and free intervention choices are added.
2. A conventional **nine-outcome qutrit Bell-basis measurement** is CPTP and includes the same source branch. But for the particular two `(3,2)` input links, six outcomes leave the remote qubits **unentangled**; this conflicts with treating every realized branch as a configuration of maximally entangled pairs.
3. A different, explicit **measure-and-reprepare** CPTP instrument has the source projector **exactly** as its first Kraus operator, and on that particular pair-state input *every nonzero-probability branch* is a product of a maximally entangled qutrit pair and a maximally entangled qubit pair. This construction preserves the pair-only description branch by branch and no-signalling when outcomes are uncommunicated. Its extra Kraus operations are nevertheless **forbidden by the source's definition of distillation**, which explicitly rejects replacement channels. It is not a derivation from the original single axiom.

4. A covariance stress test identifies an important limit: independently rotating one qutrit–qubit link's two-dimensional support can make the unchanged `P_3` branch produce an unentangled remote product state. Thus the fixed projector is not universally pair-preserving across all independently oriented supports. Whether the source permits these orientations as distinct physical configurations, or co-transforms the join under its matching rules, must be specified.

The original **172 record labels, sector counts (79,48,42,3), square-root weights, and numerical `alpha_candidate` survive unchanged in the original *selected* sector**: that sector uses the identical projector and normalized pair state. Their preservation in a **physical, unconditional observable** is not established, because the paper provides neither failure-branch insertion operators nor a rule for measuring its readout across an outcome register. The value is neither disproved nor proved to persist as a physical coupling. Different added branch rules give different unconditional numbers.

The missing dynamical scale and fermion operator identified in v2–v4 also remain missing. A one-loop beta coefficient cannot be inferred from this finite instrument calculus alone.

## 1. Exact input from the paper

A `(3,2)` link occupies the two-dimensional common support of a qutrit and a qubit. Its state, in fixed reference bases, is

\[
|\omega_2\rangle_{AB}=\frac{|0,0\rangle+|1,1\rangle}{\sqrt2},
\qquad A\cong\mathbb C^3,\quad B\cong\mathbb C^2.
\]

For two links with separately labelled qutrit copies `A1,A2` and remote qubits `B1,B2`, use

\[
|\Psi\rangle=|\omega_2\rangle_{A_1B_1}\otimes|\omega_2\rangle_{A_2B_2}
=\frac12\sum_{i,j=0}^1|ij\rangle_A|ij\rangle_B.
\]

Its remote marginal is `rho_B=I_4/4`. The paper's join is the rank-one qutrit Bell projector

\[
P=|\omega_3\rangle\langle\omega_3|,
\qquad |\omega_3\rangle=(|00\rangle+|11\rangle+|22\rangle)/\sqrt3.
\]

The qutrit-side contraction is

\[
(\langle\omega_3|_A\otimes I_B)|\Psi\rangle
=\frac1{\sqrt6}|\omega_2\rangle_B.
\]

Hence the Born success probability is `p=1/6` and the normalized selected state is

\[
|\omega_3\rangle_A|\omega_2\rangle_B.
\]

This reproduces the selected norm and the paper's `sqrt(6)` normalization in S7. The source excludes the complementary outcome as a configuration, so any operational completion including it is explicitly an added hypothesis.

### 1.1 The no-signalling lemma

For *any* bipartite density operator `rho_AB` and any local trace-preserving channel `E_A` acting on `A`,

\[
\operatorname{tr}_A[(\mathcal E_A\otimes\mathrm{id}_B)(\rho_{AB})]
=\operatorname{tr}_A\rho_{AB}.
\]

Proof: for arbitrary Bob observable `O_B`, by Heisenberg duality and `E_A^*(I_A)=I_A`,

\[
\begin{aligned}
\operatorname{tr}\bigl[(I_A\otimes O_B)(\mathcal E_A\otimes\mathrm{id})(\rho)\bigr]
&=\operatorname{tr}\bigl[(\mathcal E_A^*(I_A)\otimes O_B)\rho\bigr]\\
&=\operatorname{tr}\bigl[(I_A\otimes O_B)\rho\bigr].
\end{aligned}
\]

For selective, renormalized `P rho P / tr(P rho)` this theorem **does not apply** because that normalized conditioning is not an unconditional linear CPTP channel. Applied to the source input, the selected Bob state is `P_B=|omega_2><omega_2|`, with trace distance `D(P_B,I_4/4)=3/4`. This only implies operational signalling after adding physically independent Alice/Bob subsystems, free choice to apply the selected projection, and remote measurements; none is specified by the original source.

## 2. Completion A: two-outcome Lüders instrument

Use the two Kraus operators

\[
K_{\rm yes}=P,\qquad K_{\rm no}=I_9-P,
\qquad K_{\rm yes}^\dagger K_{\rm yes}+K_{\rm no}^\dagger K_{\rm no}=I_9.
\]

The unnormalized remote branch states are

\[
\sigma_{\rm yes}=\frac16 P_B,
\qquad\sigma_{\rm no}=\frac14I_4-\frac16P_B.
\]

The failure probability is `5/6`, and its normalized Bob state is

\[
\boxed{\rho_{B|\rm no}=\frac3{10}I_4-\frac15P_B.}
\]

The Bell-test probabilities are `1` on success and `1/10` on failure. On the full channel,

\[
\frac16(1)+\frac56\left(\frac1{10}\right)=\frac14
\]

and `sigma_yes+sigma_no=I_4/4`, as required.

**Source compatibility:** The `yes` branch is literally the original join. The `no` branch is ordinarily not a pure tensor product of maximally entangled pairs and is therefore not a state allowed by the paper's ontology. This completion restores standard operational no-signalling **by changing the allowed outcome space**.

## 3. Completion B: all nine qutrit Bell outcomes

Define the qutrit Bell basis, with `omega=exp(2pi i/3)`, by

\[
|\Phi_{t,s}\rangle
=\frac1{\sqrt3}\sum_{j=0}^2\omega^{tj}|j,j+s\bmod3\rangle,
\qquad t,s\in\{0,1,2\}.
\]

These nine vectors are orthonormal and complete in `C^3 tensor C^3`. The specific source Bell projector is the `(t,s)=(0,0)` outcome.

The source input has support only at `00,01,10,11` on Alice's two qutrits. Taking overlaps with the nine Bell vectors gives the exact probabilities and normalized entanglement characteristics:

| `s` | Number of phases `t` | Probability for each outcome | Remote qubits |
|---:|---:|---:|---|
| 0 | 3 | `1/6` | maximally entangled Bell pair; concurrence `1` |
| 1 | 3 | `1/12` | pure product state; concurrence `0` |
| 2 | 3 | `1/12` | pure product state; concurrence `0` |

The probabilities sum to `3/6+6/12=1`. The remote marginal averaged over all outcomes is `I_4/4`. The probability of producing an entangled remote pair **with an outcome label** is `1/2`. The other `1/2` belongs to outcomes incompatible with the pure-pair-only final ontology.

An outcome-dependent correction can align the Bell states in the three `s=0` branches. It cannot convert the `s!=0` remote product states into entangled pairs by a product of unitaries on the two distant qubits. Classical outcome communication and such corrections are additional operations.

## 4. Completion C: pair-preserving CPTP instrument with exact `P` success branch

This is the strongest constructive result of v6. It also demonstrates that the failure of B to preserve pure pair states is **not a theorem that all CPTP completions must fail**.

Define four orthonormal vectors embedded in the `4`-dimensional common-support subspace of Alice's two qutrits,

\[
\begin{aligned}
|b_0\rangle&=(|00\rangle+|11\rangle)/\sqrt2,\\
|b_1\rangle&=(|00\rangle-|11\rangle)/\sqrt2,\\
|b_2\rangle&=(|01\rangle+|10\rangle)/\sqrt2,\\
|b_3\rangle&=(|01\rangle-|10\rangle)/\sqrt2.
\end{aligned}
\]

Let

\[
|\eta\rangle=\sqrt3|b_0\rangle-\sqrt2|\omega_3\rangle.
\]

Direct calculation shows `||eta||=1`, `<eta|omega_3>=0`, and `<eta|b_j>=0` for `j=1,2,3`. Together,

\[
\mathcal V=\bigl\{\,|\omega_3\rangle,|\eta\rangle,|b_1\rangle,|b_2\rangle,|b_3\rangle,
|02\rangle,|12\rangle,|20\rangle,|21\rangle\,\bigr\}
\]

is an orthonormal basis for Alice's full **nine-dimensional** two-qutrit Hilbert space.

For each `v_j in V` construct the Kraus operator

\[
\boxed{K_j=|\omega_3\rangle\langle v_j|.}
\]

It obeys, **on the entire nine-dimensional input space**,

\[
\sum_{j=0}^{8}K_j^\dagger K_j
=\sum_j|v_j\rangle\langle v_j|=I_9.
\]

Moreover the selected first Kraus operator is exactly the source projector,

\[
\boxed{K_0=|\omega_3\rangle\langle\omega_3|=P.}
\]

For the original two-link input `|Psi>`, the nine outcomes have probabilities

\[
\boxed{(1/6,\;1/12,\;1/4,\;1/4,\;1/4,\;0,\;0,\;0,\;0).}
\]

Every one of the five nonzero outcomes produces `|omega_3>_A tensor |b_j>_B`, with `j=0` appearing in both the first and second outcome and the three other Bell states in the remaining outcomes. **Every realized branch on this specific source input is therefore a pair of maximally entangled pairs.**

The remote state averaged over uncommunicated outcomes is

\[
\begin{aligned}
\rho_B'&=\left(\frac16+\frac1{12}\right)|b_0\rangle\langle b_0|
+\frac14\sum_{j=1}^3|b_j\rangle\langle b_j|\\
&=\frac14\sum_{j=0}^3|b_j\rangle\langle b_j|\\
&=\boxed{I_4/4}.
\end{aligned}
\]

**What changed?** For `j>0`, the Kraus operators do not act as the paper's distillation projectors `P_m`. They measure an outcome and **reprepare the same qutrit Bell pair**, erasing the original qutrit-pair state. This is a replacement channel, precisely the kind excluded in the source's Corollary S2. The branches also require an explicit classical outcome distinction to support ensemble bookkeeping. Thus the extension is internally consistent in ordinary operational quantum mechanics, but it is **not a realization of the source's existing only-distillation axiom**.

**Scope warning:** Although the channel is CPTP on the full two-qutrit Hilbert space, the *all realized branches produce Bell pairs* result is proved for the particular two source links with fixed aligned common supports. It is not proved for arbitrary allowed supports, arbitrary matchings, arbitrary shared-qudit graphs, or other inputs. These are future covariance/uniqueness questions.


## 4A. A covariance stress test with independently oriented supports

The preceding pair-only result required both qutrit ends to occupy the **same** common-support plane `span{|0>,|1>}`. To test whether the operator works on all ordinary maximally entangled `(3,2)` links, rotate only the first qutrit end by the unitary permutation `|1> <-> |2>`, keeping the second link and the physical qutrit Bell projector `P` fixed. The pair inputs become

\[
|\widetilde\omega_2\rangle_{A_1B_1}=
\frac{|0,0\rangle+|2,1\rangle}{\sqrt2},\qquad
|\omega_2\rangle_{A_2B_2}=
\frac{|0,0\rangle+|1,1\rangle}{\sqrt2}.
\]

Each individual `(3,2)` link is still maximally entangled on a two-dimensional common support. But the same fixed `|omega_3>` qutrit join now contracts only the `A1=A2=0` component:

\[
(\langle\omega_3|_A\otimes I_B)|\widetilde\Psi\rangle
=\frac1{2\sqrt3}|00\rangle_B,
\qquad p_{\rm yes}=\boxed{\frac1{12}}.
\]

After normalization, the remote state is **`|00>`**, a product state with concurrence zero. This is an exact counterexample to the assertion that the *fixed* rank-one `P_3` filter necessarily maps every pair of independently oriented maximally entangled `(3,2)` links into another remote Bell pair. A different choice of other Kraus outcomes cannot repair **this already selected branch**.

**Interpretation must remain conditional:** The paper's matching operations are basis identifications, not freely selectable physical rotations. If the join projector must be co-transformed when a support is re-labelled, or independently oriented supports are prohibited as physical inputs, this example does not contradict the intended source model. But a future operational completion must specify that rule and test covariance rather than silently assuming the aligned input basis. In particular, being invariant under `U tensor conjugate(U)` acting *jointly* on a Bell state is not the same as preserving Bell-pair outputs for arbitrary *independent* changes of each input support.

## 5. What happens to the 172 records and the proposed `alpha(0)`?

The source constructs an algebraic **conditional readout** with

\[
N=24\cdot36\cdot36=31104,\quad
(n_1,n_{\sqrt2},n_{\sqrt3},n_{\sqrt6})=(79,48,42,3),
\]

and hence

\[
\boxed{
\alpha_{\rm cand}=
\frac{79+48\sqrt2+42\sqrt3+3\sqrt6}{31104},\qquad
\alpha_{\rm cand}^{-1}=137.035999189566155\ldots
}.
\]

**Exact selected-branch preservation:** Because the first Kraus operator in both completions is `P`, with the same selected probability `1/6` and the same normalized post-projection pair state, the original selected-sector record construction can still be applied exactly as before. Its original 172 labels, sector counts, amplitudes, and readout are unchanged **when one conditions on that first result**. This is an identity of a subexperiment, not evidence that the modified full theory still predicts a physical `alpha(0)`.

**Unconditional value requires new definitions:** A nonselective instrument is a quantum operation with an outcome register. The paper defines no readout operator on this enlarged operational space, no failure-branch sector, and no rule relating outcomes to a physical zero-momentum photon residue. For example, if one *adds* a common outcome-dependent record space of dimension `N` with insertion `M_j` and defines its readout as

\[
\overline\alpha=\sum_j p_j\frac{\operatorname{tr}_R M_j}{N},
\]

then in the simple two-outcome case the following choices are all algebraically available:

| Added choice of `M_no` | Unconditional readout | Numerical value |
|---|---|---:|
| `0` | `alpha_cand/6` | `0.001216225427277` |
| `M_yes` | `alpha_cand` | `0.007297352563662` |
| `I_R` | `alpha_cand/6+5/6` | `0.834549558760610` |

**These rows are not competing predictions**: they demonstrate underdetermination because the source does not define `M_no` or justify this added averaging observable. On the selected branch all choices preserve the identical original value.

The above direct-sum/flagged example is not itself a derived photon propagator. Physical interpretation of any unconditioned readout, and its relation to `Z_3`, still needs a separate proof. One cannot conclude that no-signalling forces a particular shift in `alpha`, or that the new measure-and-reprepare instrument automatically preserves all record weights on every branch. Neither conclusion follows from the test.

## 6. Formal obstruction and what is *not* ruled out

A CPTP channel cannot consist only of the single Kraus operator `P`, because `P^dagger P=P!=I_9`. Normalizing the selected outcome yields a nonlinear, trace-preserving **state update on the restricted domain where success has nonzero probability**, not a linear CPTP operation. Operational locality and free intervention then lead to the remote trace-distance witness from v5.

This does **not** prove the source's metaphysical pair ontology is inconsistent. It instead establishes that a standard free-choice, local, no-signalling quantum instrument cannot be identified with an always-successful selected rank-one projection *without some additional restriction or other outcomes*. Operational completions exist, as exhibited here, but they modify the original act rules or the interpretation of interventions.

Crucially, it is possible to preserve admissible pair-states branch-by-branch on this single fixed input **and** no-signalling, so a stronger claim that every trace-preserving completion must produce non-pair outcomes is false. The allowed operations and covariance under all matching transformations are the actual unresolved issues.

## 7. Relevance to the QED beta-function problem

All of these constructions are finite-dimensional and scale-free. They introduce no physical `mu`, no four-momentum, no charged fermionic determinant, no propagating photon, and no scheme-independent logarithmic response. Thus none determines

\[
\beta(\alpha)=\frac{2}{3\pi}\alpha^2+O(\alpha^3).
\]

The realistic next tasks, in decreasing order of foundational importance, are:

1. **Covariant extension check:** Can a *single fixed linear instrument*, defined without choosing the reference `|0>,|1>` support or an outcome-dependent reset rule, be covariant under the paper's whole matching group and keep all nonzero branches within the pair-only ontology for **every** admissible input? Counterexamples should use two or more independent matching rotations rather than one fixed basis.
2. **Record-instrument compatibility:** Define explicit 31,104-dimensional record-insertion operators for each outcome, not just copied eigenvalue lists. Prove whether an outcome-averaged operator can satisfy the same exchange/readout rules and preserve the 172 weighted records.
3. **Causal calculus:** Define what interventions and remote observers are physically admissible while preserving a no-signalling principle; avoid equating formal contraction order with time.
4. **Scale and propagation:** Only after 1–3 are settled, specify a physical scale transformation and derive (rather than assume) a fermionic spectrum and vacuum-polarization response. The beta coefficient remains the final diagnostic.

## 8. Reproducibility and test matrix

Run from this directory:

```bash
python chain_beta_extension_v6.py > chain_beta_extension_v6_results.json
```

Dependencies: NumPy and SymPy (standard Python 3.10+). Tests include:

- Exact `SymPy` construction of the original two `(3,2)` link input as a `9x4` amplitude matrix.
- Exact full-`9x9` orthonormality of the instrument basis and Kraus completeness.
- Exact equality of the first Kraus matrix to `|omega_3><omega_3|`.
- Exact branch probabilities and four nontrivial Bell-state conditional states.
- Exact remote marginal preservation `I_4/4` and Lüders failure state.
- Independent numerical 9-outcome qutrit Bell-basis completeness, probability, concurrence, and no-signalling checks.
- An exact independently rotated support counterexample: fixed `P_3` probability `1/12`, remote product state, concurrence `0`.
- Symbolic verification of the paper's four record-sector sums, original `N=31104` and numerical candidate.

All tests passed. No continuum limit or perturbative QED coefficient was computed, and no new physical interpretation of `alpha` is asserted.

## References

- A. D. Wissner-Gross, *A Parameter-Free Candidate for the Fine-Structure Constant from Entangled Qudits*, 23 October 2026, source PDF supplied by the user, esp. §§S1, S7, S8.
- N. Bao, A. Bouland and S. P. Jordan, *Grover Search and the No-Signaling Principle*, Phys. Rev. Lett. **117**, 120501 (2016), https://doi.org/10.1103/PhysRevLett.117.120501. Background on signalling consequences of postselection; not used to derive the formulas.
- RFC 9340, *Architectural Principles for a Quantum Internet*, §4.4.2, https://www.rfc-editor.org/rfc/rfc9340.html. Operational background on Bell-measurement outcomes and classical messaging.

**One-sentence conclusion:** The selected 172-record value can coexist algebraically with a trace-preserving instrument and a pair-only realization exists on the aligned tested input; that realization is not universally pair-preserving under independent support rotations, and standard no-signalling requires new operations or restrictions, leaving unconditional `alpha` and its beta function undetermined.