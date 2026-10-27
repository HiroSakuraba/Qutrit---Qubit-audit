# Chain Theory and the QED beta function — v8
## Gauge descent of the three-copy object: a phase-stabilizer obstruction, a partial remedy, and an exact record-sector calculation

**Date:** 24 October 2026  
**Status:** Exact finite-dimensional results and reproducible numerical checks; no independent derivation of the QED beta-function coefficient.  
**Source:** A. D. Wissner-Gross, *A Parameter-Free Candidate for the Fine-Structure Constant from Entangled Qudits* (23 October 2026), supplied as `ChainTheoryAlpha.pdf`. Consult especially the main text pp. 1–5 and Supplemental §§ S1, S3, S6–S8.  
**Companion files:** `chain_beta_extension_v8.py` and `chain_beta_extension_v8_results.json`. All v2–v7 scripts were rerun successfully in this session.

## Abstract

This note follows v7 by testing whether the proposed three-link qutrit object's Bell join and the record-exchange decomposition descend to physical input **rays** when independent link orientations are regarded as unobservable matching freedom. The answer depends on which level of structure is demanded. The full rank-one qutrit Bell projector fails to descend under an explicit one-parameter `SU(3)` stabilizer of a qutrit–qubit input ray, even though its restriction to the two occupied directions is unchanged. A phase-stabilizer twirl creates a *mixed* rank-two effect with exactly the same selected-branch probability `1/6` and remote Bell state for the aligned input, so the full-projector ambiguity does **not** alone falsify the paper's selected-branch result. A rank-two support-adapted pure join is also available, with success probability `1/4` rather than `1/6`, but requires a relative identification of the two qubit Schmidt bases. Most sharply, in the 64-dimensional qutrit–qutrit record sector, an independent ray-preserving `SU(3)` transformation changes the weighted swap readout operator, while leaving its trace `28+36 sqrt(3)` invariant. The source does not prescribe a complete three-copy gauge quotient or operator-transport law, so this is a **conditional obstruction to operator descent, not a disproof of the 172-record scalar value**.

The beta function remains undetermined because no spacetime, charged Dirac field, propagator, or scale operation has been derived from the source axiom.

## 1. What is being audited

The source's qutrit–qubit link has support dimension two. Its object comprises three links at a shared qutrit, represented as virtual copies `Q1,Q2,Q3`: `Q1,Q2` participate in two reads joined by the maximally entangled qutrit projection, and the third link `Q3-f` is unread. Its record space has dimension `N=(3^2-1)(2^2-1)(3*2)^4=24*36*36=31104`. The source assigns 172 orthogonal record directions and computes

\[
\alpha_{\mathrm{cand}}=\frac{79+48\sqrt2+42\sqrt3+3\sqrt6}{31104}.
\]

The paper **explicitly** states that each pair has its own matching and that the bases of shared-qudit copies need not align. It also treats a joint qutrit projection as always available, with no failure branch. These statements motivate an operational consistency question, not a proven inconsistency: *what is the group of admissible simultaneous transformations of the three virtual copies, their far qubits, the join projector, and the readout operators?*

Two different symmetries must be distinguished throughout:

- The **correlated Bell symmetry**, `(U ⊗ conjugate(U))|omega3> = |omega3>`, which exactly fixes the join; and
- The **independent input-ray stabilizers**, which leave individual `(3,2)` pair states unchanged as physical rays but need not fix that particular join on the complete `3×3` space.

The first is explicitly part of the paper's pair-matching mathematics. Whether the second should be **quotiented as physical gauge** when constructing the shared object is the conditional assumption under test. It is not legitimate to replace the first by the second without explanation.

## 2. A concrete stabilizer of the input pair ray

Work in an aligned qutrit frame and let

\[
|\psi\rangle_{Aq}=\frac{|0,0\rangle+|1,1\rangle}{\sqrt2},\quad
G_\phi=\operatorname{diag}(e^{i\phi},e^{i\phi},e^{-2i\phi})\in SU(3).
\]

Then

\[
(G_\phi\otimes I_2)|\psi\rangle=e^{i\phi}|\psi\rangle.
\]

Thus `G_phi` fixes the physical **ray** of the first qutrit–qubit link (it changes only its global phase). Apply it to virtual copy `Q1`, with `Q2,Q3` unchanged. The other two input link rays remain unchanged. The normalized qutrit join vector is

\[
|J_3\rangle=|\omega_3\rangle=\frac{|00\rangle+|11\rangle+|22\rangle}{\sqrt3}.
\]

Its transform is

\[
|J_3(\phi)\rangle=\frac{e^{i\phi}(|00\rangle+|11\rangle)+e^{-2i\phi}|22\rangle}{\sqrt3}.
\]

The overlap and trace distance of the two rank-one projectors are exactly

\[
\langle J_3|J_3(\phi)\rangle=\frac{2e^{i\phi}+e^{-2i\phi}}3,\qquad
D_{\rm tr}\big(P_3,P_3(\phi)\big)=\frac{2\sqrt2}{3}\left|\sin\frac{3\phi}{2}\right|.
\]

At `phi=pi/3`, this gives `2 sqrt(2)/3`, reproducing and sharpening the v7 example.

**Conditional operator-descent lemma.** If a physical operation on the two read qutrits must be a unique operator-valued function of the three **individual input pair rays**, independent of their ray-preserving transformations, the rank-one projector `P3` does not descend to that quotient. The proof is the explicit `G_phi` calculation above. It is not a contradiction if the join carries additional frame/connection data or if the source's gauge group constrains the three copy transformations.

**New qualification: the aligned selected branch is insensitive.** Define the occupied-support projector `S=(|0><0|+|1><1|)^{⊗2}`. Direct multiplication gives

\[
\boxed{SP_3(\phi)S=SP_3S=\tfrac23|\omega_{2,\mathrm{occ}}\rangle
\langle\omega_{2,\mathrm{occ}}|}\quad\text{for every }\phi.
\]

The global full-space projectors differ, but **their effects on the original aligned pair input do not**. The selected branch therefore still has probability `1/6` and leaves Bob with `|omega2>`. The operator-descent defect concerns the unused qutrit direction, other support orientations, and the extended record algebra. It cannot be advertised as evidence that the original selected numerical answer changes.

## 3. An invariant mixed effect preserves the original selected branch

Averaging the full qutrit Bell projector over the above ray stabilizer phase removes interference between the occupied two-plane and the unused one-dimensional complement. It suffices to average the three phases `phi=0, 2*pi/9, 4*pi/9`. Denote `|b>=|omega_{2,occ}>` and `|c>=|22>`. Exactly,

\[
\boxed{E_{\rm stab}=\tfrac23|b\rangle\langle b|+\tfrac13|c\rangle\langle c|.}
\]

This is positive, has trace one, rank two, and purity `tr(E_stab^2)=5/9`. It is a valid measurement **effect** because `0≤E_stab≤I9`; it is **not** a rank-one qutrit Bell projection. On the occupied support,

\[
S E_{\rm stab} S=S P_3 S.
\]

Let `Psi` be the `9×4` coefficient matrix of the two original aligned `(3,2)` pair states. The unnormalized remote state of the corresponding postselected local measurement effect is

\[
\boxed{\Psi^\dagger E_{\rm stab}\Psi
=\tfrac16|\omega_2\rangle\langle\omega_2|.}
\]

Thus the selected branch has the **same** probability `1/6` and **same** conditional remote Bell state as the original projector. The finite three-phase average was independently verified numerically with Frobenius residual `1.14×10^{-16}`.

**Limit:** This success concerns the original aligned support. For supports `span(e0,e2)` on `Q1` and `span(e0,e1)` on `Q2`, the same effect instead yields

\[
\Psi_{\rm rotated}^\dagger E_{\rm stab}\Psi_{\rm rotated}
=\tfrac1{12}|00\rangle\langle00|,
\]

which is a product remote state. The effect has therefore **not** solved the arbitrary-support transport problem. Nor has a full trace-preserving, pair-preserving instrument or a rule for the 172 record weights of its other outcomes been derived. This is a conditional replacement measurement effect, not a correction implicitly required by the original paper.

## 4. The support-only pure join: another option with a different normalization

Let `Vi=[vi0,vi1] : C2 -> C3` be the isometries encoding two input link states in specified qubit Schmidt bases. A natural pure join on the **occupied two-planes** is

\[
|J_2(V_1,V_2)\rangle=\frac{v_{10}\otimes v_{20}+v_{11}\otimes v_{21}}{\sqrt2}.
\]

This construction is covariant under independent `SU(3)` transformations of the qutrits. If the input representatives change by overall phases `Vi -> exp(i phi_i) Vi`, then

\[
|J_2\rangle\mapsto e^{i(\phi_1+\phi_2)}|J_2\rangle,
\]

so its *projector* is phase-independent. Direct contraction of the two `(3,2)` links gives

\[
\boxed{\langle J_2|\psi(V_1)\psi(V_2)\rangle=	frac12|\omega_2\rangle},\quad
p_2=\boxed{\tfrac14}.
\]

Compare the transported full qutrit join,

\[
\langle J_3(V_1,V_2)|\psi(V_1)\psi(V_2)\rangle=\frac1{\sqrt6}|\omega_2\rangle,
\quad p_3=\boxed{\tfrac16}.
\]

Both produce a maximally entangled remote qubit pair, but they are physically different conditional operations. The support-only join has Schmidt rank two, not three, and cannot inherit the original three-direction `sqrt(3)` record weights without **recomputing the entire readout**.

**Crucial restriction:** the formula for `J2` pairs `vi0` with `v20` and `vi1` with `v21`, and hence presupposes identification of the *qubit Schmidt bases*. Independent changes of those bases generally change this projector. Projective invariance here has been proven for overall phases of the pair rays in fixed qubit coordinate systems, not for independent arbitrary internal Schmidt-frame relabelings. This is a *partial* solution, not a unique ray-only operational construction in the broadest sense.

128 seeded random complex `SU(3)` support pairs were tested, including arbitrary overall phases of the input link representatives. The maximum probability/normalized-output deviation was `1.71×10^{-15}` and the maximal pair-ray phase-equivariance residual was `5.22×10^{-16}`.

## 5. Central charge does not resolve the relative unused-direction phase

The source's generator on the direct sum is

\[
Q=i\operatorname{diag}(3,3;-2,-2,-2)=6iY.
\]

On the six-dimensional `(3,2)` tensor link it acts as `i I6`, and its normalized left-action trace form gives `B(Q,Q)=1`. Thus it produces a global phase on that link; it does **not** produce the relative phase between the occupied qutrit two-plane and its unused complement appearing in `G_phi`. In particular, scalar central transformations cannot convert `P3(phi)` into `P3(0)` while leaving the rest of a fixed full-space join untouched.

There is a possible *correlated* compensation: allow phases `phi,psi` on qutrit copies `Q1,Q2`. The rank-one Bell projector is then preserved precisely when

\[
\boxed{e^{3i(\phi+\psi)}=1.}
\]

The continuously connected choice `psi=-phi` is an instance of the correlated `U ⊗ conjugate(U)` Bell stabilizer. It would amount to a specified relationship between the phases of the read copies. The source does not provide an explicit three-copy gauge connection establishing this relationship for independently oriented input pairs, nor does it prove that its global conserved count enforces it.

## 6. A direct calculation in 64 of the 172 record directions

The first two qutrit read slots each have eight traceless Gell-Mann letters. The sector with one letter in `Q1` and one in `Q2` therefore spans

\[
\mathcal H_{QQ}\cong\mathbb C^8\otimes\mathbb C^8,
\quad\dim=64.
\]

The read-exchange operation in this sector is the factor swap `S(a⊗b)=b⊗a`. Its symmetric and antisymmetric eigenspaces have dimensions

\[
\dim \mathcal H_+=\frac{8(8+1)}2=36,
\qquad
\dim \mathcal H_-=\frac{8(8-1)}2=28.
\]

Following the source's survival-weight prescription (main text Eq. 5), the `X=-1` directions receive amplitude weight `sqrt(3)` while the `X=+1` directions receive weight one in this same-end, `Y=+1` sector. **Importantly, Supplemental §S7 sets `X=-S` on these same-end products, not `X=S`.** Thus the symmetric 36-dimensional subspace has weight `sqrt(3)` and the antisymmetric 28-dimensional subspace has weight 1. The corresponding finite-dimensional **sector model** is

\[
P_\pm=\frac{I\pm S}2,\qquad
M_{QQ}=\sqrt3P_++P_-,\quad
\operatorname{tr}M_{QQ}=28+36\sqrt3.
\]

Take the same input-ray stabilizer `G_{pi/3}` on `Q1`. Its adjoint action on the eight Gell-Mann letters, in the explicitly tested orthonormal basis ordering, is

\[
R=\operatorname{diag}(1,1,-1,-1,-1,-1,1,1),\quad \operatorname{tr}R=0.
\]

On the two qutrit-letter slots its action is `W=R⊗I8`. Define `M'_{QQ}=W M_{QQ} W†`. The trace is unchanged, as it must be under simultaneous conjugation:

\[
\boxed{\operatorname{tr}M'_{QQ}=\operatorname{tr}M_{QQ}=28+36\sqrt3.}
\]

But the operator itself **is not unchanged**. Exactly,

\[
\begin{aligned}
\|M'_{QQ}-M_{QQ}\|_{\mathrm{HS}}^2
&=(\sqrt3-1)^2\|W P_- W^\dagger-P_-\|_{\mathrm{HS}}^2\\
&=(\sqrt3-1)^2\frac{8^2-|\operatorname{tr}R|^2}{2}\\
&=\boxed{32(\sqrt3-1)^2}.
\end{aligned}
\]

There is an explicit expectation-value illustration. For the **fixed coordinate probe** `v=a⊗a`, with `a=(e_1+e_3)/sqrt(2)` in the eight-letter basis,

\[
\langle v|M_{QQ}|v\rangle=\sqrt3,\qquad
\langle v|M'_{QQ}|v\rangle=\frac{1+\sqrt3}{2}.
\]

If the probe is simultaneously transported `v -> Wv`, its expectation remains `sqrt(3)`. **This is precisely the distinction between covariance and invariance.** The difference is not an operational prediction for a gauge-invariant probe unless the source defines what stays fixed during its proposed frame transformation.

By contrast, the correlated sign action `R⊗R` commutes with the swap `S` and preserves this sector operator. This explicit case is compatible with the paper's correlated matching symmetry.

**Scope of the computation:** The `64×64` operator is a direct model of the paper's qutrit–qutrit same-end sector, using its exchange and weight rules. It does **not** reconstruct all 172 orthonormal record vectors in the full 31,104-dimensional readout space. We have **not** proven that the fully transported 172-record operator fails to be defined, or that the paper's scalar normalized trace changes. The source assigns other sectors and commutator records, for which it supplies counts and weights but no physical three-copy matching/measurement prescription of the kind audited here.

## 7. The required three-copy gauge-descent diagram

Let `g=(g1,g2,g3)` act on the qutrit virtual copies, with corresponding matched transformations (where defined) on the far qubit endpoints. Write `S_input` for the physical equivalence class of the three link rays and `M(S_input)` for the proposed readout insertion.

A strong formulation of physical operator descent would require

\[
S_{\rm input}(g\cdot\Psi)=S_{\rm input}(\Psi)
\quad\Longrightarrow\quad
M(g\cdot\Psi)=M(\Psi)
\]

when both operators are expressed in the same fixed physical reference frame. A weaker **covariant-family** requirement instead allows

\[
M(g\cdot\Psi)=W_{\mathcal R}(g)M(\Psi)W_{\mathcal R}(g)^\dagger
\]

with a specified action `W_R` on every one of the full readout space's tensor factors and insertions. In the latter case, the scalar normalized trace is automatically invariant by cyclicity, but that property alone supplies neither an operationally gauge-invariant **observable** nor the quantum-field-theory photon propagator.

The paper's distinction between independent link matchings, correlated pair symmetries, shared-copy identification, and the conserved charge must determine which of these requirements is appropriate. The v8 calculations prove that the distinction matters **already in the qutrit-sector swap**, before any attempt to derive a scale-dependent beta function.

## 8. What was verified and what remains open

| Test | Result | Meaning |
|---|---|---|
| Input-pair ray stabilizer under `G_pi/3` | Exact | A nontrivial subgroup acts trivially on the link *ray* |
| Full qutrit Bell join on that quotient | Projector trace distance `2 sqrt(2)/3` | Not uniquely defined as a fixed full-space operator without further relational structure |
| Occupied restriction of the join | Identical exactly | Original aligned selected branch unaffected |
| Stabilizer-twirled effect | Rank two, purity `5/9`, `p=1/6` | Mixed effect preserves selected aligned branch, but not an original rank-one join |
| Support-only pure join | `p=1/4` | Projective-phase remedy with changed normalization and fixed Schmidt-basis identification |
| Qutrit QQ exchange sector (64 directions) | `28+36 sqrt(3)` trace; operator changes by HS²=`32(sqrt(3)-1)^2` | Trace covariance does not imply operator invariance |
| Correlated qutrit transformation | Preserves 64-sector swap | Correlated and independent frame operations must not be conflated |
| Full 172 directions under complete three-copy gauge action | Not established | Needs a formal transport law and definition of physical operator equivalence |
| QED one-loop beta function from the axiom | Not derived | Requires spacetime, fermion dynamics, physical scale |

## 9. v9 research agenda: a genuinely discriminating test

1. **Specify the admissible groupoid:** Define independent link matchings, the shared-copy identification, phase/ray equivalence, and transformations of the two read links plus unread state link as explicit objects and arrows. Do not assume all `SU(3)^3` operations are automatically gauge; the source's ``one qudit, multiple virtual copies'' language might restrict them.
2. **Construct the induced action on all 172 record operators, not just their count.** Prove either equivariance under the actual allowed groupoid or produce one complete counterexample that changes an operationally gauge-invariant expectation. The `64×64` obstruction supplies an exact test fixture.
3. **Check observationally equivalent completions:** Try the support-only rank-two join, the ray-stabilizer mixed effect, and the transported rank-three join as competing physically valid outcomes. Derive **all** readout insertions and probabilities from each instrument, rather than substituting `2` for `3` in the old counts.
4. **Require causal and spectral structure:** Even full operator descent would not give a beta function. A separate derivation of propagating charged degrees of freedom and an RG scale is indispensable.

**Evidence standard:** Agreement with `alpha(0)` from the original 172-record count is not enough; any preferred completion must be selected by the original axiom or an explicitly declared new principle, and the running coefficient must not be smuggled in via a conventional Dirac determinant.

## 10. Reproducibility

Execute from this folder:

```bash
python chain_beta_extension_v8.py > chain_beta_extension_v8_results.json
```

Dependencies: Python 3, NumPy, SymPy. The file uses exact SymPy matrix arithmetic for the stabilizer, phase-averaged effect, selected remote state, and 64-sector swap calculations. NumPy checks the finite phase average, the actual Gell-Mann adjoint, and 128 seeded independent complex `SU(3)` support pairs. All v2–v7 scripts were rerun successfully on 24 October 2026.

### Background (not an input to the calculation)

- Correlated `U⊗conjugate(U)` is the isotropic-pair symmetry, distinct from independent local twirling. See e.g. [Families of invariant states, Sec. 6](https://www.felixleditzky.info/teaching/FT25/ch6-invariant-states.pdf).
- For context on tensor-network gauge constraints and gauge-invariant operators, see [Gauged Gaussian projected entangled pair states, Phys. Rev. D 110, 054511 (2024)](https://doi.org/10.1103/PhysRevD.110.054511).

No result from these external sources is used to insert a beta coefficient or to amend the source paper's definitions.