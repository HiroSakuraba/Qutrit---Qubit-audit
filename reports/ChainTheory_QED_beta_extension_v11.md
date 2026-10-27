# Chain Theory → QED beta program, v11: full readout embedding and gauge-descent audit

**Date:** 2026-10-24  
**Primary source:** A. D. Wissner-Gross, *A Parameter-Free Candidate for the Fine-Structure Constant from Entangled Qudits* (2026-10-23), attached as `ChainTheoryAlpha.pdf`.  
**Status:** Explicit finite-dimensional realization and an exact stabilizer/record-invariance mismatch. **Not** a derivation of spacetime, electrodynamics, or the first QED beta-function coefficient.

## 0. Summary and what changed since v10

The next unresolved step in v10 was to reconstruct all 172 record operators with their unread-link dependence in the full `N=31104`-dimensional readout and test them under the original three-copy matching transformations. This report supplies **one explicit, natural realization**:

\[
\mathcal R=\underbrace{\operatorname{Herm}_0(\mathbf C^3)\otimes
\operatorname{Herm}_0(\mathbf C^2)}_{24\;\mathrm{dimensions\ at\ unread\ }Q_3f}
\otimes
\underbrace{\operatorname{End}(\mathbf C^3\otimes\mathbf C^2)}_{36\;\mathrm{at\ read\ }1}
\otimes
\underbrace{\operatorname{End}(\mathbf C^3\otimes\mathbf C^2)}_{36\;\mathrm{at\ read\ }2}.
\]

The Hermitian/traceless real factors are understood as complexified when defining the operator-space Hilbert–Schmidt representation. They are *operator spaces*, **not** the physical Hilbert space of a 31,104-level particle. With a normalized connected component `c` of the unread state and 172 orthonormal read-arm record directions `r_j`, write

\[
|R_j)=|c)\otimes|r_j),\qquad
M_{\rm v11}=\sum_{j=1}^{172}w_j|R_j)(R_j|.
\]

The script constructs the columns `|R_j)` explicitly as a `31104 × 172` complex array, checks the column norms and selected Gram entries, verifies the whole 172-dimensional read-arm eigensystem, and shows that the trace of this low-rank insertion reproduces the **paper's** number with no fitted parameters:

\[
\operatorname{tr}_{\mathcal R} M_{\rm v11}
=79+48\sqrt2+42\sqrt3+3\sqrt6,
\qquad
\frac{\operatorname{tr}M_{\rm v11}}{31104}=0.007297352563662260\ldots.
\]

**Most important new result:** A nontrivial matched `SU(2)` transformation can leave the entire **normalized joined three-pair state** invariant yet move the **fixed** 172-record projector of this embedding. In a concrete exact example its Hilbert–Schmidt displacement squared is 4. The trace stays unchanged under covariant transport.

**Key qualification and a correction to the v9 interpretation:** Supplemental §S7 expressly says “under every matching” attaches to **reads**, not necessarily to each individual record remaining fixed. Consequently, failure of the *fixed-coordinate record projector* to descend is not by itself a contradiction of the source. The result reveals what a physical equivalence/covariance prescription would need to specify.

## 1. Precise source inputs

The source defines:

1. A qutrit–qubit link with two-dimensional common support and state `|ω₂⟩`.
2. A three-arm object consisting of two read links plus an unread state link, with the two read copies joined by a qutrit Bell projection `P₃` (Supplemental §S7, Eq. S2).
3. Two read slots at each of the two read links, with traceless-letter dimensions `(8,3,8,3)` (Supplemental §S7).
4. The `169` two-slot products plus `3` antisymmetric qubit commutator differences (Theorem S6).
5. The exchange-sector sizes `79,48,42,3` and weights `1,√2,√3,√6` (Lemma S10, Theorem S7).
6. The unread 24-dimensional traceless-traceless operator factor and two full 36-dimensional read-arm operator factors (Lemma S8).
7. “The record operator at the two read links times the connected part of the unread link” as a direction in `\mathcal R` (Supplemental §S7, after its record counting).

The text does **not** prescribe a momentum scale, causal dynamics, full gauge-field connection, or how to compare every independently chosen virtual-copy frame as a physical instrument. These remain open.

**Formalization boundary:** In Supplemental §S11 the author expressly notes that the Lean file computes the record-weight sum and the readout dimension **without constructing the readout space itself**. The present `31104×172` construction therefore tests an additional realization that the paper describes mathematically but does not formally build in Lean. It is an independent consistency exercise, not a claim to reproduce a previously checked Lean construction.

## 2. Exact connected state on the unread qutrit–qubit arm

Use ordered computational basis `|00⟩,|01⟩,|10⟩,|11⟩,|20⟩,|21⟩` on `\mathbf C^3\otimes\mathbf C^2`. Let

\[
|\omega_2\rangle=(|00\rangle+|11\rangle)/\sqrt2,\qquad
\rho=|\omega_2\rangle\langle\omega_2|.
\]

Define the unique orthogonal projector onto the traceless-traceless sector by subtracting the marginals and restoring the scalar part:

\[
C = \rho-\rho_3\otimes\frac{I_2}{2}
-\frac{I_3}{3}\otimes\rho_2+\frac{I_6}{6},\qquad
\rho_3=\mathrm{tr}_2\rho,\quad \rho_2=\mathrm{tr}_3\rho.
\]

Both partial traces of `C` vanish identically. In the chosen basis,

\[
C=
\begin{pmatrix}
\tfrac14&0&0&\tfrac12&0&0\\
0&-\tfrac14&0&0&0&0\\
0&0&-\tfrac14&0&0&0\\
\tfrac12&0&0&\tfrac14&0&0\\
0&0&0&0&0&0\\
0&0&0&0&0&0
\end{pmatrix},
\qquad
\mathrm{tr}(C^2)=\frac34.
\]

These are **exact SymPy identities** (the last identity is independent of any floating-point tolerance). Define `c=C/√(3/4)`. Expand `c` in the `24` Hilbert–Schmidt orthonormal tensor products of the eight Gell-Mann and three Pauli directions. The recovered 24-vector `c_vec` has norm one; the reconstruction error is below machine precision.

For a local basis transformation `G⊗u`, the projection defining `C` is equivariant, since taking partial traces commutes with local conjugation:

\[
C((G\otimes u)\rho(G\otimes u)^\dagger)
=(G\otimes u)C(\rho)(G\otimes u)^\dagger.
\]

For the paper's matched local stabilizer `G=diag(\bar u,1)` (with `u∈SU(2)`), the pair state, and hence `C`, is individually fixed. The numerical covariance residual for an additional non-stabilizing `SU(3)×SU(2)` basis change is `≈10⁻¹⁶`.

## 3. All 172 read directions, without a full 31104×31104 matrix

Use orthonormal Gell-Mann and Pauli operator bases at the four read slots `Q1,q1,Q2,q2` and normalized identity operators at unoccupied slots. The read-arm physical Hilbert space has dimension `3·2·3·2=36`; its operator Hilbert space has dimension `36²=1296`.

Exactly `175` elementary mutually orthonormal directions are convenient: the `169` two-slot products plus **six** single-qubit directions (three at each read qubit). The original `172` records are a subspace of this 175-dimensional ambient space, excluding the three symmetric sums of the single-qubit directions.

We construct an explicit orthonormal 172-column matrix `E` in this ambient space:

- `Q1Q2`: symmetric `36` directions of weight `√3` and antisymmetric `28` of weight `1`.
- `q1q2`: symmetric `6` directions of weight `√3` and antisymmetric `3` of weight `1`.
- Two pairs of mixed-slot blocks: `48` symmetric directions of weight `√2` and `48` antisymmetric directions of weight `1`.
- Qubit commutator differences: `3` directions of weight `√6`.

This gives `42` at weight `√3` and `79` at weight `1`, with the source's other two sectors. A direct numerical check verifies that each of the 172 columns is an eigenvector of the v9 weight matrix with the claimed eigenvalue, and `E†E=I₁₇₂`.

Embed these actual read-arm **operators** into the 1296-dimensional operator space as columns `r_j`; then form their tensor products with the unread connected direction:

\[
V=[\,c\otimes r_1\quad\cdots\quad c\otimes r_{172}\,]
\in\mathbf C^{31104\times172},\qquad V^\dagger V=I_{172}.
\]

The script instantiates all `5,349,888` complex entries of `V` (a `31104×172` array), verifies all 172 column norms and a well-distributed eight-column Gram sample. The algebraic tensor-product factorization and the earlier full 172-vector Gram test in the read sector together prove the entire Gram identity, without multiplying a huge 31,104-square matrix.

The resulting finite-rank insertion and subspace projector are

\[
M=V\operatorname{diag}(w_1,\ldots,w_{172})V^\dagger,
\qquad P=VV^\dagger.
\]

The number of additional orthogonal directions in `\mathcal R` is `31104−172=30932`. Thus the source's denominator is reproduced at its stated operator-space level, not chosen after the calculation.

**Scope:** Taking the same normalized connected state direction for all 172 read operators is the most direct implementation of the source's “times the connected part” wording. The source does not uniquely formalize a gauge-dependent family of all these operator directions; this is an **explicit candidate realization**, not an independent proof that the axiom fixes it.

## 4. What matching invariance means, and what it cannot mean

The source's own Supplemental §S7 (paragraph “Why the exchange eigenbasis”) distinguishes:

- **Individual fixed record directions** under every matching, which it explicitly declines to require;
- **Reads** understood under every matching, which it chooses so that the traceless basis directions remain meaningful.

This distinction must be respected. Transporting `P` and `M` along with the chosen operator bases is always covariant under the induced unitary/orthogonal action `W`:

\[
P\mapsto WPW^\dagger,\quad M\mapsto WMW^\dagger,
\quad\operatorname{tr}(WMW^\dagger)=\operatorname{tr}M.
\]

Such an equality of traces **does not prove** that an invariant gauge-quotient operator exists: an operator on the gauge quotient would need the appropriate independence from the chosen representative, or a definition as a section of an associated bundle together with a rule for comparing measurement outcomes.

There is therefore no contradiction between the source's invariance claim *for the numerical readout* and the following failures of fixed-coordinate projector invariance. The latter are useful tests for an operational completion, not proofs against the original scalar construction.

## 5. Full three-pair state stabilizer versus fixed commutator-record projector

This is the most substantial new comparison.

After the source's selected qutrit projection, its normalized three-pair state is

\[
|\Omega\rangle=
|\omega_3\rangle_{Q_1Q_2}\otimes
|\omega_2\rangle_{q_1q_2}\otimes
|\omega_2\rangle_{Q_3f}.
\]

It is the `12`-term state given by the normalization of Eq. (S2); the original read links enter through the *pre-projection* state. We now test transformations that are **matched pair stabilizers of the pre-projection read links**, and whether they also preserve `|Ω⟩` and the **fixed** record subspace.

Take a qubit `SU(2)` transformation `u_i` at `q_i`, with corresponding qutrit matching `G_i=diag(\bar u_i,1)` at `Q_i`. Each such pair transformation fixes its corresponding pre-projection `|ω₂⟩` exactly. The following two prescriptions preserve different structures:

**A. Identical qubit frames:** `u₂=u₁`.

The Pauli adjoint rotations are identical, so the three commutator differences `σ_a(q1)−σ_a(q2)` transform within their own 3D space. Consequently the entire fixed 172-record projector is preserved by this pair of matchings.

But the qutrit join `|ω₃⟩` need **not** remain invariant, because `G₁⊗G₂` need not have the `G⊗\bar G` form required to stabilize that Bell state. For the exact choice

\[
u_1=u_2=\operatorname{diag}(i,-i),
\qquad G_1=G_2=\operatorname{diag}(-i,i,1),
\]

we obtain

\[
\big|\langle\omega_3|(G_1\otimes G_2)|\omega_3\rangle\big|^2
=\frac19,
\quad
\lVert P_{172}-WP_{172}W^\dagger\rVert_{\rm HS}^2=0.
\]

The join is genuinely changed by this transformation when held in the same physical coordinate frame, even though it preserves the selected commutator space.

**B. Conjugate qubit frames:** `u₂=\bar u₁`.

This preserves **both** the qutrit Bell join and the remote qubit Bell pair. The unread `Q₃f` pair can be left untouched. Therefore it leaves the entire normalized `|Ω⟩` exactly invariant.

Choose

\[
u_1=\operatorname{diag}(e^{i\pi/4},e^{-i\pi/4}),
\qquad u_2=\bar u_1.
\]

The qubit Pauli adjoint rotations differ by a relative 180-degree rotation. As shown in v9, the fixed three-dimensional *difference* subspace then rotates in the ambient six-dimensional single-qubit sector, yielding the **exact** result

\[
\boxed{
(G_1\otimes G_2\otimes u_1\otimes u_2\otimes I_6)|\Omega\rangle=|\Omega\rangle,
\qquad
\|P_{172}-WP_{172}W^\dagger\|_{\rm HS}^2=4.
}
\]

The code directly constructs the `216`-component normalized source-inspired state in the tensor order `Q1,Q2,q1,q2,Q3,f`. The state-stabilizer norm error was `1.31×10⁻¹⁷`, while the fixed-record projector shift was `3.9999999999999982`.

Because the matched transformation fixes the unread link as well, the very same subspace displacement holds for our full `31104`-dimensional embedding, using `P_full=|c)(c|⊗P_read`.

**Inference:** Even a transformation stabilizing the **entire** normalized joined state does not necessarily stabilize this particular fixed-coordinate record subspace. This is a mathematically explicit gauge-descent issue for a fixed physical-instrument interpretation.

**Not inferred:** It does *not* show that the source's scalar numerical coupling is gauge-dependent; the covariantly transported weighted trace is unchanged. Nor does it rule out a coherent operator-bundle interpretation in which the record subspace itself transforms along with the frame.

## 6. Why the two natural frame correlations differ

There are two independent matching relationships in the fixed coordinates:

- To keep the qubit commutator-difference subspace fixed, the two Pauli *adjoint* rotations must coincide: `Ad(u₁)=Ad(u₂)`, or `u₂=±u₁` in `SU(2)`.
- To preserve the qutrit Bell join under the pre-read-link matching relation `G_i=diag(\bar u_i,1)`, one needs `G₂=\bar G₁`, equivalent to `u₂=\bar u₁`.

These conditions do not coincide for general `u∈SU(2)`. Their continuous intersection includes the real `SU(2)` matrices, forming an `SO(2)` subgroup; the v11 script tests a nontrivial real rotation and verifies preservation of both structures.

This does **not** establish that the physical gauge group reduces to `SO(2)`. It establishes only that treating **both** the qutrit join and the individual commutator-difference projector as separately fixed in one coordinate chart is more restrictive than allowing covariant transport.

## 7. All-three-arm covariance of the constructed 31104D realization

With arbitrary independent local basis changes on the three qutrit–qubit arms, let

\[
W_{\mathrm{unread}}=\mathrm{Ad}(G_3)\otimes\mathrm{Ad}(u_3)
\]

act on the 24-dimensional unread connected operator space, and let `W_read` be the adjoint action of the transformations at `Q1,q1,Q2,q2` on the read-arm operator space. The explicit candidate's full record columns transform as

\[
V'=\left(W_{\mathrm{unread}}\otimes W_{\mathrm{read}}\right)V.
\]

This yields `M'=WMW†` and exactly invariant trace and rank. The program checks a representative **non-stabilizing** local change on the unread arm, an independent paired matching on the first read arm, and direct conjugation on all `175` elementary read-arm operator directions. It directly compares the `31104`-entry tensor coordinates for 12 representative full record columns by two independent computational routes:

1. Transform the underlying physical qutrit/qubit operators by direct matrix conjugation, then expand in the operator basis.
2. Apply the induced adjoint representation matrices to the `24⊗1296` record coordinates.

The maximum absolute column discrepancy across the 12 selected records is `6.86×10⁻¹⁶`. All 172 record columns are covered algebraically because they are explicit linear combinations of the 175 elementary directions, each of which is checked, and the unread factor is common.

This confirms the **covariance of the constructed family**. It does not establish the existence of a canonical physical comparison map between different pair-copy frames in the original theory.

## 8. Audit table: proven, conditional, open

| Statement | Status |
|---|---|
| `tr(C²)=3/4`, `C` has zero partial traces | Exact symbolic result |
| Explicit 172 orthonormal record directions inside the `24·36·36` algebra | Constructed; factorization proof and numerical checks |
| Weighted multiplicities and numerical `α_candidate` | Reproduced from source's definitions |
| Full three-pair state invariant under conjugate matched read-link rotations | Exact algebra; numerical check |
| Same matched rotations necessarily fix the chosen 172-record projector | **False** for the fixed-coordinate embedding; exact counterexample |
| Record projectors and weighted insertion transport covariantly | True in the constructed operator-bundle realization |
| Gauge-independent scalar trace under such transport | True; follows from cyclicity |
| Unique physical identification of the independent frames | Not derived |
| Uniqueness of the selected 172-record instrument from original axiom | Not derived |
| Complete no-signaling measurement instrument or propagator | Not derived |
| Physical spacetime/energy scale and one-loop QED coefficient | Not derived |

## 9. What to do next (v12)

1. **Specify the quotient.** Formalize precisely whether independent pair matchings should be viewed as passive changes of the local operator frame, physical symmetries of the postselected state, or restrictions on allowed interventions. These are mathematically different uses of “gauge.”
2. **Test an instrument-valued naturality condition.** Define a measurement instrument `\mathcal I_\rho` or a family of readout maps depending functorially on each state-and-frame configuration. Require `\mathcal I_{g\rho g^\dagger}=g\circ\mathcal I_\rho\circ g^{-1}` (with the appropriate domain and outcome labels). This is stronger than invariance of a scalar trace.
3. **Resolve the stabilizer.** Determine whether state-preserving `SU(2)` rotations of the normalized joined object may rotate the commutator outcome labels, and if so construct the admissible associated-bundle representation of the records. Identify which transformations correspond to physically distinct record outcomes.
4. **Prove or refute uniqueness of the frame connection.** With the physical equivalences fixed, classify possible relational transports; do not assign arbitrary Wilson loops to the three-copy comparison graph and call them emergent gauge curvature.
5. **Return to the beta coefficient only after dynamics.** Neither an exact `α(0)` nor a covariant finite-dimensional readout supplies a mass scale, electron field, four-dimensional propagation, or the universal logarithmic vacuum polarization.

## 10. Reproducibility

From the folder containing the bundle files, run:

```bash
OPENBLAS_NUM_THREADS=2 python chain_beta_extension_v11.py
```

Dependencies: NumPy, SymPy, and `chain_beta_extension_v9.py` (included). Output: `chain_beta_extension_v11_results.json`.

The script allocates a `31104×172` complex array, not a dense `31104×31104` array. It verifies the exact connected-state projection, the complete read-arm eigensystem and orthonormality, norms of every full record column, a Gram sample, covariance of the unread arm, direct physical transformation of all 175 read-arm elementary operators, full three-arm covariance for 12 complete record directions, and both join/commutator stabilizer tests.

All suite scripts **v2 through v11** were rerun without failures in the same Python environment. Their individual regression logs are included in the bundle. The previous test suites exercise their own stated finite-dimensional constructions and should not be interpreted as independent verification of claims beyond those constructions.