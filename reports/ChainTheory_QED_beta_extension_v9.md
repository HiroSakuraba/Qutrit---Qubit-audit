# Chain Theory → QED beta function: v9 gauge-descent audit of all 172 records

**Research date:** 2026-10-24  
**Primary source:** A. D. Wissner-Gross, *A Parameter-Free Candidate for the Fine-Structure Constant from Entangled Qudits* (2026-10-23), uploaded as `ChainTheoryAlpha.pdf`.  
**Reproducibility:** `chain_beta_extension_v9.py`, `chain_beta_extension_v9_results.json`. Earlier v2–v8 test suites rerun without failure.  
**Status:** Exact finite-dimensional obstruction to a **fixed** 172-dimensional record-space representation under independent matched local basis transformations; this is **not** a disproof of a covariant-family interpretation of the paper and **not** an independent derivation of the QED beta function.

## Executive finding

The paper's 172-record basis consists of **169 two-slot products plus 3 antisymmetric qubit commutator differences**. The three differences occupy a 3D subspace inside the 6D space of one-slot qubit letters from the two read links. An allowed pair-preserving local matching on the first read link can rotate two of those differences into the three *excluded* symmetric-sum directions. Thus the paper's selected 172D subspace is **not an invariant subspace of independent link matchings**, even though its normalized trace can be transported covariantly.

This gives a precise dichotomy: (i) a *diagonal/correlated* identification of the two read frames preserves the 172-record operator; or (ii) a *covariant family* of record subspaces, related by explicit transport, is needed. Treating the original 172D subspace as fixed across independent matching frames is incompatible with the finite operator algebra.

## 1. Source input and canonical operator embedding

Supplemental §S7, Corollary S6, Theorem S6, Eqs. (S3)–(S5), give slots `(8,3,8,3)`—two qutrit read slots `Q1,Q2` and two qubit slots `q1,q2`. The 169 distinct two-slot words consist of:

| Slot pair | Count |
|---|---:|
| Q1-Q2 | 64 |
| q1-q2 | 9 |
| Q1-q1, Q1-q2, Q2-q1, Q2-q2 | 4 × 24 = 96 |
| **Two-slot total** | **169** |

The three additional records are explicitly the qubit commutator **differences** `[s_μ,s_ν]_{q1} - [s_μ,s_ν]_{q2}`. Since `[σ_i,σ_j]=2iε_{ijk}σ_k`, after normalized basis changes they span exactly the differences of corresponding Pauli letters.

We build a transparent 175D ambient subspace of the read-arms' **operator** space:

\[
\mathcal V_{175}=\mathcal P_{169}\oplus (\mathcal V_{q1}\oplus\mathcal V_{q2}),
\qquad \dim\mathcal V_{q1}=\dim\mathcal V_{q2}=3.
\]

The explicit operators act on `C3(Q1) ⊗ C2(q1) ⊗ C3(Q2) ⊗ C2(q2)`, a 36D Hilbert space whose full operator space has dimension 1296. Normalized Gell-Mann and Pauli letters are used with normalized identities on unused slots; all 175 basis operators are **Hilbert–Schmidt orthonormal**, verified by their explicit Gram matrix.

To connect with the paper's 31,104D readout algebra, these are the read-arm components of the record operators, with a fixed normalized connected/unread-link factor in the 24D unread sector. The 175D analysis is sufficient to test closure under transformations acting only on the read arms; it is **not** an explicit construction of dynamical gauge fields or all interactions in the source.

The paper's selected commutator space is

\[
\mathcal V_- = \{(v,-v):v\in\mathbb R^3\},\qquad
P_- = \frac12\begin{pmatrix}I_3&-I_3\\-I_3&I_3\end{pmatrix},
\]

and the 172D record-space projector is

\[
P_{172}= I_{169}\oplus P_-.
\]

The complementary 3D symmetric space is `V_+={(v,v)}` and is **not** one of the paper's commutator records.

## 2. Reconstructing the complete weighted readout

The paper defines two exchanges `X,Y`. On two letters at the same end, `X` reads minus their swap parity (`X=-S`), with `Y=+1`; on mixed-end letters, `X=+1` and `Y` reads minus the qubit-slot swap parity. The commutator differences have `X=Y=-1`. The paper's amplitudes are therefore exactly reproduced by the following operator on the 175D ambient space (zero on `V_+`):

- On the 64 qutrit-qutrit words, `M_QQ=√3 P_sym + P_antisym`; symmetric/antisymmetric dimensions 36/28.
- On the 9 qubit-qubit words, `M_qq=√3 P_sym + P_antisym`; dimensions 6/3.
- On the 96 mixed words, `M_mixed=√2 P_sym + P_antisym` for the exchange of qubit slots, 48 dimensions apiece.
- On the six one-slot qubit letters, `M_comm=√6 P_-` (three nonzero directions, three zero directions).

The complete spectrum on the 175D ambient subspace is:

| Eigenvalue | Multiplicity |
|---|---:|
| 0 (extra excluded directions) | 3 |
| 1 | 79 |
| √2 | 48 |
| √3 | 42 |
| √6 | 3 |

This reconstructs the source's 172 nonzero record weights and provides a concrete `175×175` representation of the record insertion. Its trace is **exactly**

\[
\operatorname{tr}M=79+48\sqrt2+42\sqrt3+3\sqrt6.
\]

Using the source's normalization (`N=24·36·36=31104`),

\[
\alpha_{\rm cand} = \frac{\operatorname{tr}M}{31104}
=0.00729735256366226\ldots.
\]

No other numerical normalization is introduced in the audit.

## 3. An exact matching that fixes the physical link ket but moves two records out of the subspace

Choose an element of `SU(2)` on `q1` together with a matching `SU(3)` element on `Q1`:

\[
V_{q1}=\begin{pmatrix}i&0\\0&-i\end{pmatrix},\qquad
V_{Q1}=\operatorname{diag}(-i,i,1).
\]

Both determinants are +1. On the two-dimensional occupied support, the qutrit operation is the complex conjugate of the qubit operation. Therefore the two-qudit pair **ket itself**, not merely its ray, is fixed:

\[
(V_{Q1}\otimes V_{q1})|\omega_2\rangle_{Q1q1}=|\omega_2\rangle_{Q1q1}.
\]

This is a matched local basis transformation of the form used in the source's §S1: it changes neither the entangled pair state nor its Schmidt support. The second read link is untouched.

However, adjoint conjugation of the first qubit's Pauli letters gives

\[
(\sigma_x,\sigma_y,\sigma_z)_{q1}\mapsto
(-\sigma_x,-\sigma_y,+\sigma_z)_{q1}.
\]

Consequently, the commutator-difference records obey

\[
\begin{aligned}
\sigma_x(q1)-\sigma_x(q2)&\mapsto-\sigma_x(q1)-\sigma_x(q2),\\
\sigma_y(q1)-\sigma_y(q2)&\mapsto-\sigma_y(q1)-\sigma_y(q2),\\
\sigma_z(q1)-\sigma_z(q2)&\mapsto \sigma_z(q1)-\sigma_z(q2).
\end{aligned}
\]

The first two outcomes lie **orthogonally outside** the original commutator-difference subspace. Hence the original 172D record space is not closed under this particular **independent** matched link transformation. The accompanying `SU(3)` operation rotates qutrit-letter product sectors among themselves and does not repair the lost single-qubit directions.

Writing `W` for the 175D induced adjoint action, this produces the exact identities

\[
\boxed{\left\|P_{172}-WP_{172}W^\dagger\right\|_{\rm HS}^2=4},
\]

and

\[
\boxed{\left(\|P_{172}Wr_x\|^2,\|P_{172}Wr_y\|^2,\|P_{172}Wr_z\|^2\right)=(0,0,1)}.
\]

The `6×6` qubit-sector identities have an independent exact SymPy verification. For the full 175D action, direct conjugation of **all 175 explicit operators** on the underlying 36D Hilbert space agrees with the constructed matrix `W` to Frobenius residual `3.44×10^-15`.

## 4. General formula for all relative qubit-frame rotations

Under local adjoint rotations `R1,R2∈SO(3)` on the two qubit slots, the one-letter sector acts by

\[
T_{R_1,R_2}=\begin{pmatrix}R_1&0\\0&R_2\end{pmatrix}.
\]

The difference space `V_-` is invariant **if and only if** `R1=R2`. Indeed, the transform of `(v,-v)` is `(R1v,-R2v)`, which is again a difference for every `v` exactly when the rotations coincide. This shows the obstruction for the full continuous family, not only the single example above.

Define `R=R1^T R2` and let `θ` be its standard SO(3) rotation angle, so `tr R=1+2 cos θ`. Direct block multiplication gives

\[
\begin{aligned}
\operatorname{tr}(P_- T P_- T^T)&=\frac{3+\operatorname{tr}R}{2}=2+\cos\theta,\\
\left\|P_- - TP_-T^T\right\|_{\rm HS}^2&=3-\operatorname{tr}R=2(1-\cos\theta).
\end{aligned}
\]

In the 175D ambient record algebra, the 169 product directions are closed under these slotwise adjoint rotations, so the same projector-distance formula applies to the 172D selected space. Its minimal invariant linear closure under independently admissible qubit rotations adds the **three symmetric-sum directions**: `169 + 6 = 175`. This claim is about the stated finite operator model; it does not assert uniqueness of a full physical completion.

Random matched `SU(2)` tests (20 seeded Haar draws) verified the general distance formula with maximum absolute discrepancy `3.11×10^-15`. The matching transformations kept the pair ket fixed with worst residual `5.57×10^-16`.

## 5. Three inequivalent readout prescriptions: only one belongs to the original paper

**Covariantly transported operator.** With `M(g)=W M W†` and the transformed 172D subspace `P(g)=WP172W†`, the trace is unchanged:

\[
\boxed{\operatorname{tr}(WMW^\dagger)=\operatorname{tr}M}.
\]

The source's numerical candidate is therefore unchanged **if** its record spaces and observable are transported consistently. An explicit prescription for this transport across the join and unread link is still required for a physical interpretation.

**Fixed-coordinate re-projection (new test, not a prediction).** If the basis is transformed and the resulting operator is projected **back** onto the untransformed fixed 172D record subspace, the three commutator directions are not all retained. Precisely,

\[
\operatorname{tr}\big(P_{172}WMW^\dagger P_{172}\big)
=\operatorname{tr}M-\sqrt6(1-\cos\theta).
\]

The corresponding *hypothetical* normalized discrepancy is

\[
\boxed{\Delta\alpha_{\rm fixed}=\frac{\sqrt6(1-\cos\theta)}{31104}}.
\]

For the `π` rotation, `Δα=2√6/31104 = 0.0001575031984814...`, and the original fixed-basis readout would shift from `0.00729735256366226` to `0.00713984936518083`, a 2.16% change. **This is not an observed change of the fine-structure constant.** It is an explicit demonstration that frame-handling prescriptions change answers unless a gauge-covariant one is specified. The original paper does not prescribe this re-projection operation.

**Independent local twirl (another new test).** Averaging the three difference-record projectors over four matched `SU(2)` Pauli rotations acting on one read link gives

\[
\boxed{\frac14\sum_{g\in\{I,i\sigma_x,i\sigma_y,i\sigma_z\}}
T_g(\sqrt6P_-)T_g^\dagger=\frac{\sqrt6}{2}I_6.}
\]

Its trace stays `3√6`, but all six one-slot directions receive equal weight `√6/2` instead of three directions having weight `√6` and three having zero weight. This **does not preserve the original record content**; it is only an example of an invariant completion of this sector.

## 6. What correlated matching accomplishes

If the **same** qubit-frame transformation acts on both read qubits, the difference records remain differences. If the **same** qutrit-frame transformation also acts on both read qutrit copies, the same-end exchange swaps and mixed qubit-slot swaps commute with the induced action. In the explicit 175D reconstruction,

\[
W_{\rm diag}P_{172}W_{\rm diag}^\dagger=P_{172},\qquad
W_{\rm diag}MW_{\rm diag}^\dagger=M.
\]

Both identities passed the numerical tests for the matched `π` operation (residual < `10^-12`). They hold algebraically for simultaneously applied identical adjoint rotations because swaps commute with diagonal tensor actions and the antisymmetric qubit-letter subspace is invariant.

**Qualification:** These equations concern the **two read arms**. A complete three-copy gauge action should also specify what happens on the unread link's 24D connected factor and on the qutrit–qutrit join `P3`. We have **not** proved that the original paper enforces diagonal frame identification on all three copies. Indeed, its §S1 speaks of independent matchings per pair, while its §S7 uses one shared qutrit for link exchange. Reconciling these roles is part of the unresolved theory.

## 7. Consequences for the photon conjecture and beta-function research

This v9 result is stronger than merely finding a noninvariant matrix element in the qutrit QQ sector (v8): for independently matched link frames, the **selected linear record subspace itself** need not be invariant. The smallest obvious repairs are distinct theories or interpretations:

1. **Relational/diagonal matching:** postulate a connection identifying the two read frames (and consistently the third unread frame), under which the source's selected records are globally well-defined. This is additional structure unless it is derived from the shared-qudit axiom.
2. **Covariant record bundle:** allow a frame-dependent 172D record subspace with a unitary transition rule and specify gauge-invariant expectation values and join compatibility. The numerical trace can remain fixed, but the physical readout needs an operational definition.
3. **Enlarge and average:** include all six one-slot qubit directions, then build an independently invariant insertion. This changes the paper's record selection and potentially its interpretation, even if trace matching can be engineered.

**What this does not prove:** It does not show that the photon conjecture is false, that the original numerical expression must be changed, or that the QED beta function is derived. The original v2–v8 underdetermination remains: a momentum scale, propagating charged fermions, and suitable continuum dynamics have not followed from the single qudit axiom.

**Next discriminating project (v10):** formulate a **three-copy relational gauge connection** as a discrete structure, calculate its transport around loops of shared qudits, and test whether its holonomies are determined by the pair states or remain free physical/gauge data. If the connection is not fixed, produce an exact counterexample with identical pair-state data but different 172-record physical readouts. Distinguish a transported scalar trace (automatically invariant) from a gauge-invariant measurement procedure.

## 8. Computational audit and test scope

- The 175 underlying read-arm matrices (dimension 36 each) were explicitly built and verified Hilbert–Schmidt orthonormal to <`1e-13`.
- Eigenvalue multiplicities exactly matched the integer sector counts `79,48,42,3`, plus `3` zero directions. The algebraic form of all weights follows from paper §S7; the floating-point eigenspectrum independently checks the implementation.
- The explicit matched local `SU(3)×SU(2)` operation fixes the original two-qudit pair ket exactly; both matrices have determinant one.
- Its conjugation of all 175 physical read-arm operator directions matched the induced 175D representation with Frobenius error `3.44e-15`.
- The 6D commutator-subspace leakage, fixed-projection commutator trace, and complementary ranks were checked **exactly in SymPy**.
- Twenty seeded random matched `SU(2)` tests and four continuous-angle values verified the general relative-rotation formulas.
- v2, v3, v4, v5, v6, v7, and v8 scripts were rerun successfully in the same container session.
- All computations concern finite matrices, representation transport, and the paper's stated record-counting rules—not photon scattering amplitudes, early-universe time variation, or a first-principles QED renormalization-group calculation.