# Chain Theory → QED beta function: v10, the three-copy connection and holonomy audit

**Date:** 2026-10-24  
**Primary source:** A. D. Wissner-Gross, *A Parameter-Free Candidate for the Fine-Structure Constant from Entangled Qudits* (2026-10-23), `ChainTheoryAlpha.pdf`.  
**Immediate predecessors:** `ChainTheory_QED_beta_extension_v8.md`, `ChainTheory_QED_beta_extension_v9.md`.  
**Reproduction:** `python chain_beta_extension_v10.py`, producing `chain_beta_extension_v10_results.json`.  
**Status:** A consistent *additional* relational gauge-connection calculus is built and tested. Its holonomy is **not determined by the original finite pair data**. This does not constitute a no-go theorem against all possible dynamical completions and does not derive QED running.

## Abstract

The source's qutrit–qubit link is a rank-2 maximally entangled state in `C^3 ⊗ C^2`, and at the photon-candidate object the qutrit is represented by three virtual copies, two read and one unread. The source grants each pair its own basis matching. Its 172-record count is computed with specific exchange identifications. v9 showed that the three commutator-difference record directions fail to remain a **fixed** subspace under independently applied matched pair symmetries, although they transform covariantly. We now introduce explicit relative transports between three copy frames, prove that the relevant commutator subspace transports covariantly, and test closed products of those transports. Pure-gauge transports have trivial loop holonomy, whereas nonflat transports have gauge-invariant, different Wilson traces, without changing the underlying three-link pair state. Because the comparison transports do not appear in the source's axiom or pair-state tensors, neither their presence nor their loop holonomy is fixed by that input. A support-overlap construction provides a candidate partial transport but is nonunique as a *unitary* at rank-deficient overlaps. The 172-record trace is preserved when transported consistently; its invariant value alone does not prove that the original record operator is a well-defined gauge-independent observable on physical states.

## 1. Original model: precisely what is available

The source specifies:

1. Links of type `(3,2)` in the maximally entangled state `|ω₂⟩ = (|00⟩+|11⟩)/√2`, with the third qutrit direction unoccupied on each link.
2. Three *virtual copies* of a shared qutrit, one for each link; the two read links are joined while the third link is unread (Supplemental S1 and S7).
3. Local matched basis actions `U ⊗ \bar U` that preserve equal-dimensional Bell states and corresponding support-preserving actions on a `(3,2)` pair.
4. A normalized trace of 172 record directions inside a readout algebra of dimension `N=31104`. Its value is

   `tr M = 79 + 48√2 + 42√3 + 3√6`,  
   `αcand = tr M / 31104 = 0.007297352563662262…`.

The source does **not** construct a physical spacetime, a propagation Hamiltonian, a gauge-link connection, an energy-dependent photon propagator, or a loop holonomy. In particular its three readout arms do not form a triangle of physical network links. Those are scope limitations, not numerical disagreements with the paper.

## 2. Each pair admits an independent, exact stabilizer

Let `u ∈ SU(2)` act on the qubit and let

`G₃(u) = diag(\bar u, 1) ∈ SU(3)`

act on the qutrit with occupied subspace spanned by `|0⟩,|1⟩`. Directly,

`(G₃(u) ⊗ u)|ω₂⟩ = |ω₂⟩`.

For three copies the operations are independent: `g₁,g₂,g₃ ∈ SU(2)` give the identical three-link product ket. Sixty random Haar tests establish this to norm error below `1.2×10⁻¹⁵`.

This stabilizer preserves the *fixed physical pair ket* while rotating the basis of qubit letters. A record such as `σₓ(q₁) − σₓ(q₂)` is therefore not invariant under separate transformations of the two links. That is the exact v9 obstruction, not evidence that the Bell pair itself changes.

## 3. A covariant connection resolves the commutator-frame problem

A single inter-read **adjoint** transport `T ∈ SO(3)` identifies the three real Pauli-letter directions. Define

\[
P_-(T)=\frac12\begin{pmatrix}
 I_3 & -T\\ -T^{\mathsf T} & I_3
\end{pmatrix}.
\]

It projects onto the three-dimensional subspace `{(v, −Tᵀv):v∈R³}`. By `TTᵀ=I`,

\[
P_-(T)^2=P_-(T),\qquad\operatorname{tr}P_-(T)=3.
\]

If independently chosen local frames act by `R₁,R₂ ∈ SO(3)`, the transport transforms as `T → R₁TR₂ᵀ`. Therefore

\[
\operatorname{diag}(R_1,R_2)P_-(T)\operatorname{diag}(R_1,R_2)^{\mathsf T}
=P_-(R_1TR_2^{\mathsf T}).
\]

Thus a **family** of commutator record subspaces can be gauge covariant. The amplitude-weighted insertion `M_comm(T)=√6 P_-(T)` has fixed trace `3√6`. Seventy independent random rotation tests had maximum covariance error `1.25×10⁻¹⁵` and projector idempotence error `1.30×10⁻¹⁵`.

Crucial distinction: this *defines* a relative transport. The original qudit pair data do not independently give the matrices `T`; the formula solves covariance once `T` has been supplied.

## 4. Three-copy comparison transports and loop holonomy

The next section is explicitly a **proposed relational completion**. Place `SU(2)` transports on three oriented comparison edges between the occupied support frames:

\[
T_{12},T_{23},T_{31}\in SU(2),\qquad
T_{ij}\longmapsto g_iT_{ij}g_j^\dagger.
\]

The ordered comparison-loop product

\[
H=T_{12}T_{23}T_{31}
\]

transforms as `H → g₁Hg₁†`. Thus the Wilson trace `W=tr H` is invariant under independently selected gauges at all three copies. This is standard parallel-transport algebra, **not** a prediction of a physical Wilson loop from the original model.

If the transports arise *entirely* from three frame choices, `Tij=g_i g_j†`, they are pure gauge, and `H=I₂` identically. This is a useful negative result: changing individual matching frames, by itself, **cannot produce gauge-invariant nontrivial holonomy**.

But the same three Bell-pair state can be furnished with an additional, nonflat transport assignment. Compare:

| Proposed comparison transports | Holonomy | `tr_SU2 H` | `tr_SO3 Ad(H)` |
|---|---|---:|---:|
| `(I,I,I)` | `I` | `2` | `3` |
| `(I,I,exp(−iπσ_z/4))` | quarter-turn | `√2` | `1` |

The difference is gauge-invariant within the proposed extended theory. Both assignments leave all three original maximally entangled pair kets unchanged *because they introduce extra relational labels*. The numerical tests used 60 independent sets of local `SU(2)` matchings, and verified pure-gauge flatness, nonflat holonomy covariance and Wilson traces to errors at or below `1.7×10⁻¹⁵`.

**Interpretation:** The original pair-state specification does not distinguish the assignments. Hence *if* a connection is added as an independent object, its curvature/holonomy must be restricted by further dynamics or axioms. The result is **not** a proof that no mathematical construction could derive a preferred flat connection from stronger definitions of sharing.

## 5. Uniqueness stress test: a candidate support-overlap connection can be singular

A possible idea is to identify support frames of two links with their overlap matrix. For qutrit-side isometries `V_i:C²→C³`, define

\[
S_{12}=V_1^\dagger V_2.
\]

The polar decomposition of a full-rank `S12` gives a uniquely defined unitary polar factor. But rank deficiency occurs for **allowed pairs** of two-dimensional qutrit supports. For example,

\[
V_1=(e_0,e_1),\quad V_2=(e_1,e_2)
\]

have

\[
S_{12}=\begin{pmatrix}0&0\\1&0\end{pmatrix},\qquad
\operatorname{sv}(S_{12})=(1,0).
\]

Its polar **partial isometry** is unique on the supported subspace, but completing it to a *full unitary frame transport* leaves a phase on the missing dimension:

\[
U_\phi=\begin{pmatrix}0&e^{i\phi}\\1&0\end{pmatrix},\quad \phi\in[0,2\pi).
\]

All `Uφ` are unitary and agree on the first column while differing on the second. The script tests four distinct completions and verifies the rank, singular values, and equal action on the supported column. Thus the overlap/polar approach cannot define a **unique, everywhere-valid** full `U(2)` connection without an additional rule at singular overlaps. A partial isometry may still be the correct physical object; that alternative requires its own network composition law.

There is an additional scope condition: computing `Vi†Vj` assumes that the two qutrit copies already share a common ambient `C³` identification. The original text describes independent matchings of the copies; a physical definition of cross-copy inner products is precisely what is under investigation.

## 6. Full 172-record readout and the unread third arm

The earlier v9 code builds a 175-direction **ambient, source-inspired read-arm embedding** with 169 two-slot product operators, three antisymmetric commutator differences, and three excluded symmetric sums. In it,

\[
M_{175}=M_{169}\oplus \sqrt6P_-(I),\quad
\operatorname{spec}(M)=\{1^{\times79},\sqrt2^{\times48},\sqrt3^{\times42},\sqrt6^{\times3},0^{\times3}\}.
\]

Covariantly transporting both read arms leaves the rank 172, spectrum, and trace unchanged; the v10 test gives trace difference `1.42×10⁻¹³` (floating-point). It therefore reproduces the **numerical** record calculation under a consistent transported-bundle prescription.

The unread link contributes a 24D connected traceless operator sector `(su(3) ⊗ su(2))` in the source. v10 includes its independent local matched action `Ad(G₃)⊗Ad(u₃)`. An *illustrative embedding* of a particular record operator with one chosen normalized unread-sector factor changes under this independent transformation: for the seeded example, the original and transported factors have squared overlap `0.2260144667`, giving a full tensor-operator Hilbert–Schmidt displacement squared `493.802770239` in the illustrative model.

**Do not overinterpret this numerical value.** The illustrative single-factor embedding is not a reconstruction of the source's **individual unread-link record operators**; the source's 172 directions were not individually shown to descend to the quotient by independent gauges on all three copies. This remains an unresolved mathematical task. Nor was the original 31,104D readout rewritten as a completely gauge-invariant physical instrument.

## 7. What is now proved, and what remains open

**Proved within the explicit finite-dimensional comparison model:**

- Each original `(3,2)` pair possesses an exact matched local `SU(2)` stabilizer; three such stabilizers are independent on the original product of three virtual links.
- The commutator-difference rank-3 record family becomes covariant when an `SO(3)` relative transport is *specified*.
- Purely frame-derived comparison connections are flat (`H=I`).
- Additional nonflat comparison connections possess distinct gauge-invariant Wilson traces while leaving original local Bell-pair data unchanged.
- A proposed support-overlap/polar construction does not deliver a unique full unitary transport on all allowed rank-2 qutrit support pairs.
- The source-inspired 172-record spectral trace survives covariant transport. This does not imply a unique operational instrument or loop holonomy.

**Unproved:**

- That a physical, spacetime-associated loop exists in the original three-arm object.
- That shared-copy identity enforces a unique connection (flat or curved) across all admissible links.
- That the 172 *individual* record operators, including the unread-arm dependence, descend to well-defined physical observables under every independent matching.
- That a unique no-signaling operational completion is selected, or that this induces causal evolution.
- That charged propagating fermions, energy/momentum, vacuum polarization, or the QED one-loop beta coefficient emerge from these axioms.

## 8. Recommended v11 research

1. Reconstruct each of the 172 *full* record operators in the entire `N=31104` readout (or prove an equivalent representation-theoretic descent), not only the read-arm 175D embedding.
2. Specify the gauge quotient before invoking local overlap matrices. Examine whether the shared-qudit identification reduces `SU(2)^3` to diagonal `SU(2)` or instead requires a connection.
3. If a relational connection is added, derive a composition/cocycle law and compare flat and nonflat sectors using gauge-invariant Wilson traces, including singular-support transitions.
4. Only then investigate a local transfer operator with well-defined scale and spectrum. The QED beta function requires a propagator or equivalent scale-dependent response, none of which is present yet.

## External context

The covariance law for parallel transport and the gauge invariance of Wilson traces are standard gauge-theory results; see C. Rovelli, *Gauge theory and gravity in the loop formulation* (CERN hosted PDF), and the 2026 review *Foundations of Wilson Lines and Loops From Basic Concepts to Gauge Field Reconstruction* (Fortschritte der Physik, DOI 10.1002/prop.70156). Their mathematics motivates this extension but is **not** evidence that the source paper already defines gauge holonomy.

## Reproducibility notes

Run `python chain_beta_extension_v10.py`. It imports `chain_beta_extension_v9.py` and writes `chain_beta_extension_v10_results.json`. Its exact SymPy checks include the rank-3 projector and `SU(2)` Wilson trace `√2`. Its numerical checks use deterministic seeds and evaluate 60 independent three-copy gauge transformations, 70 random adjoint projector-covariance tests, the transported 172-record weighted operator, the independent unread-link action, and four distinct unitary completions of a singular support-overlap matrix. It does not assume any experimental value of `α` to select a connection.