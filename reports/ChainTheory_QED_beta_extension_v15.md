# Chain Theory → QED beta program, v15
## A local covariant physical instrument, its exact swap-algebra reduction, and the remaining gauge/normalization obstruction

**Date:** 25 October 2026  
**Primary source:** A. D. Wissner-Gross, *A Parameter-Free Candidate for the Fine-Structure Constant from Entangled Qudits*, 23 October 2026 (`ChainTheoryAlpha.pdf`), especially Letter Eqs. (2)–(6) and Supplemental §§S1, S6–S8, S11.  
**Immediate dependencies:** v9/v11 explicit 172-dimensional representative record embedding; v13 literal `P3` branch on the 216-dimensional three-pair system; v14 unital CP Stinespring bridge.  
**Status:** New **conditional mathematical construction and exact counterchecks**; not a deduction of the physical electromagnetic coupling or its QED beta function from the paper's axiom.

## Abstract

We impose locality to the two joining qutrits, preservation of the original qutrit Bell projection as an **exact Kraus branch**, complete positivity and trace preservation, no-signaling under a local instrument, and covariance under jointly transported reference frames.

Starting from the *particular* 31,104-dimensional operator-Hilbert-space embedding used in v11 and the unital CP map built in v14, we localize the resulting five-outcome physical POVM by a normalized partial trace over the two read qubits and the unread qutrit–qubit pair. All four nonzero effects reduce **exactly** to elements of the two-qutrit swap algebra, `span{I9,S}`; their sum is `(172/31104)I9`. They commute with the qutrit Bell effect `P3`. Consequently, a **physical, local, six-outcome CPTP instrument** exists with selected Kraus operator `K_Bell=P3` and the record-score outcomes assigned only after the complementary branch.

For the stand-alone localized record POVM, the paper's number `alpha_cand` occurs *exactly* on the maximally mixed nine-dimensional qutrit state, or on any state whose swap expectation is `1/3`. But the specified unjoined input has swap expectation `1/2`, producing `0.007338539835821646...`, and the selected Bell state has swap expectation `1`, producing `0.007462101652299803...`. Thus the source's physical pure states do **not** produce `alpha_cand` with this choice of local bridge. The six-outcome Bell-first instrument gives `0.006094856227105012...` on the unjoined input if its Bell event is assigned score zero. That extra score convention is not supplied by the source, and the result also disagrees with the paper's candidate.

Most importantly, the record effects are covariant when transported with the reference frame, but are not **fixed invariant** under all correlated `U⊗conj(U)` matchings that preserve `P3`. A simple exact witness `U=diag(i,-i,1)` proves this. A family of equally local, CP, covariant completions changes the score continuously without changing the selected Bell outcome. The result is therefore an improved operational bridge **and** a precise demonstration that the axiom has not uniquely selected the physical measurement or its score.

## 1. Spaces, source assumptions, and what is additional

For the original three qutrit–qubit links, retain the tensor-factor labels `Q1,q1,Q2,q2,Q3,f`. The standard external operational extension uses

\[
\mathcal H_{\rm phys}
=(\mathbb C^3\otimes\mathbb C^2)^{\otimes3},
\quad\dim\mathcal H_{\rm phys}=216.
\]

Alice's prospective local join is on `A=Q1Q2`, a 9-dimensional subsystem. The other two read qubits are `B=q1q2` (dimension 4), and the unread pair is `C=Q3f` (dimension 6). The source's readout space is instead

\[
\mathcal R=\operatorname{End}(\mathbb C^6)^{\otimes2}\otimes
[\mathfrak{su}(3)_{\rm Herm}\otimes\mathfrak{su}(2)_{\rm Herm}]_{\mathbb C},
\quad\dim\mathcal R=36^2\times24=31104.
\]

The original paper defines a **normalized trace on this operator space** and *conjecturally* identifies it with `alpha(0)`; it does not define a physical POVM whose Born probabilities carry the four weights. The paper also excludes the complementary outcome to its Bell projection. The present constructions add a conventional CPTP measurement calculus, a representative basis embedding and a normalized conditional expectation. These additions must never be called consequences of its sole axiom.

The record operator `M` has spectral weights `0,1,sqrt(2),sqrt(3),sqrt(6)` with multiplicities `30932,79,48,42,3` respectively. Hence

\[
\alpha_{\rm cand}=\frac{79+48\sqrt2+42\sqrt3+3\sqrt6}{31104}
=0.007297352563662260\ldots
\]

## 2. From the full-space POVM to a local measurement

The v14 Stinespring map is constructed from a Hilbert–Schmidt orthonormal operator frame `{A_a}` representing the source-motivated `R⊂End(H_phys)` and its exact tight-frame identity

\[
\sum_{a=1}^{31104}A_a^\dagger A_a=144I_{216}.
\]

For each record sector `w` (and the zero sector), v14 gives a physical POVM effect `E_w^global=Phi(P_w)` on the full 216-dimensional space. Its nonzero effects factor as

\[
E_w^{\rm global}=\frac{S_w\otimes C_*^\dagger C_*}{144},
\quad S_w=\sum_{r\in w}R_r^\dagger R_r,
\quad \operatorname{tr}(C_*^\dagger C_*)=1.
\]

Define a **localizing quantum channel in the Heisenberg picture** by normalized partial trace over `B` and `C`:

\[
\mathcal L(X)=\frac{\operatorname{Tr}_{BC}X}{24}
\quad\text{for }X\in\operatorname{End}(A\otimes B\otimes C).
\]

It is unital, completely positive and preserves identity. Applied to the v14 effects, it produces

\[
E_w^A=(\mathcal L\circ\Phi)(P_w)
=\frac{\operatorname{Tr}_{q_1q_2}S_w}{4\cdot6\cdot144}
=\frac{\operatorname{Tr}_{q_1q_2}S_w}{3456}.
\]

The `q1q2` partial trace of each sector was computed independently from **all 172 explicit record eigenoperators** using NumPy, then matched against simple exact rational `9×9` matrices in SymPy. Maximum entrywise difference: `1.52e-18`.

Let `S|a,b>=|b,a>` denote the swap on `C3⊗C3`. The exact answer is

\[
\boxed{\begin{aligned}
E_1&=\frac{55I_9-7S}{20736},\\
E_{\sqrt2}&=\frac{I_9}{648},\\
E_{\sqrt3}&=\frac{77I_9+21S}{62208},\\
E_{\sqrt6}&=\frac{I_9}{10368}.
\end{aligned}}
\]

The swap has eigenvalue `+1` on the 6-dimensional symmetric sector and `-1` on the 3-dimensional antisymmetric sector. The effects have these exact eigenvalues:

| Effect | Symmetric (6×) | Antisymmetric (3×) |
|---|---:|---:|
| `E_1` | `1/432` | `31/10368` |
| `E_sqrt2` | `1/648` | `1/648` |
| `E_sqrt3` | `49/31104` | `7/7776` |
| `E_sqrt6` | `1/10368` | `1/10368` |

Every eigenvalue is nonnegative. Their remarkable sum is **scalar**:

\[
\boxed{
E_1+E_{\sqrt2}+E_{\sqrt3}+E_{\sqrt6}
=\frac{172}{31104}I_9.}
\]

The cancellation in the swap term is exact: `-7/20736+21/62208=0`. Thus the complementary effect is also scalar,

\[
E_0=\frac{30932}{31104}I_9.
\]

The total chance of a nonzero record event is therefore state-independent, but its distribution across the weighted outcomes depends on the swap expectation.

## 3. Exact physical-state evaluation

For any two-qutrit state `rho_A` let

\[
s=\operatorname{tr}(\rho_A S),\qquad -1\le s\le 1.
\]

Then

\[
\begin{aligned}
p_1(s)&=\frac{55-7s}{20736},\\
p_{\sqrt2}(s)&=\frac1{648},\\
p_{\sqrt3}(s)&=\frac{77+21s}{62208},\\
p_{\sqrt6}(s)&=\frac1{10368},\\
p_0(s)&=\frac{30932}{31104}.
\end{aligned}
\]

A weighted-score convention gives

\[
\boxed{
f(s)=\frac{55-7s}{20736}
+\frac{\sqrt2}{648}
+\frac{(77+21s)\sqrt3}{62208}
+\frac{\sqrt6}{10368}.}
\]

This score is **an expectation of numeric labels**, not the probability of a distinguished measurement outcome, and need not equal the physical electromagnetic coupling.

For the specified unjoined input, tracing out both read qubits gives

\[
\rho_A^{\rm initial}=\rho_Q\otimes\rho_Q,
\quad\rho_Q=\tfrac12\operatorname{diag}(1,1,0).
\]

The swap identity gives `s=tr(rho_Q^2)=1/2`. The selected Bell state is `rho_A^Bell=P3=|omega3><omega3|` and has `s=1`. The maximally mixed state `I9/9` has `s=1/3`.

| State | Exact `s` | Score `f(s)` |
|---|---:|---:|
| `I9/9` | `1/3` | `0.007297352563662260` |
| Original unjoined input | `1/2` | `0.007338539835821646` |
| Selected `P3` branch | `1` | `0.007462101652299803` |

The deviation for the input state is **exactly**

\[
\boxed{f(1/2)-\alpha_{\rm cand}
=\frac{7(\sqrt3-1)}{124416}
=0.00004118727215938578\ldots}
\]

or about `+0.5644%` relative. Matching `alpha_cand` requires `s=1/3`. Even a family of other valid three-pair inputs can produce this: if two rank-two qutrit support projections overlap along one direction and have principal-angle cosine satisfying `cos(theta)^2=1/3`, then `s=[1+cos(theta)^2]/4=1/3`. This demonstrates sensitivity to an unselected support relation, not a new parameter-free prediction.

## 4. A CPTP *local* instrument that preserves the literal `P3` branch

The Bell join projector

\[
P=P_3=|\omega_3\rangle\langle\omega_3|,
\quad F=I_9-P,
\]

satisfies `S|omega3>=|omega3>`. Thus **every `E_w` commutes with `P`**, and the two measurements are compatible as POVMs.

Define the six Kraus operators on `A`:

\[
\boxed{K_{\rm Bell}=P,\qquad
K_w=\sqrt{E_w}\,F\quad(w=0,1,\sqrt2,\sqrt3,\sqrt6).}
\]

Each acts as the identity on `B,C`. Their completeness is exact:

\[
\sum_jK_j^\dagger K_j
=P+F\Bigl(\sum_wE_w\Bigr)F
=P+F=I_9.
\]

So the total instrument is CPTP and local. On the original three-pair input,

\[
\boxed{\Pr({\rm Bell})=\operatorname{tr}(\rho_A^{\rm initial}P)=\tfrac16}
\]

and its normalized selected branch is precisely the source's joined pair state. Unconditional no-signaling follows from trace preservation and was verified explicitly on Bob's two qubits to `2.78e-17` max entrywise error.

If the selected Bell outcome is assigned **score zero**, and the record weights apply only on the complementary branch, the weighted expectation becomes

\[
\begin{aligned}
 f_{\rm seq}(\rho_A)
 &=\sum_{w>0}w\operatorname{tr}(\rho_A F E_wF),\\
 f_{\rm seq}(\rho_A^{\rm initial})
 &=f(1/2)-\frac16 f(1)\\
 &=\boxed{0.006094856227105012\ldots}.
\end{aligned}
\]

This **does not** reproduce the original candidate. Assigning the Bell outcome another score would alter the expectation; the source supplies no physical convention to select that score. In particular, choosing it to force agreement would be a calibration, not a derivation.

## 5. A sharper gauge test: transported covariance versus fixed invariance

There are two different requirements:

1. **Covariant families:** If `rho -> W rho W†`, move the full instrument along with the frame, `K_j -> W K_j W†`. All probabilities are unchanged, by cyclicity of the trace. Thirty-two numerical random correlated `SU(3)` tests passed, with max score-branch probability residual below `2.9e-15`. This is a coherent extension but uses an external frame-transport prescription.
2. **Fixed physical effects:** Keep `P` and each `E_w` fixed while invoking a matching transformation that is supposed to act as a gauge relabeling. This is stronger and **fails** for the localized record measurement even under a transformation that leaves the Bell join itself invariant.

Here is an exact counterexample. Set

\[
U=\operatorname{diag}(i,-i,1)\in SU(3),\qquad W=U\otimes\overline U.
\]

Then `W P W†=P` exactly, the familiar maximally entangled-pair matching invariance. But

\[
\boxed{\|WSW^\dagger-S\|_{\rm HS}^2=16.}
\]

Because `E_1` has swap coefficient `-7/20736` and `E_sqrt3` has coefficient `+7/20736`, it follows that

\[
\boxed{
\|WE_1W^\dagger-E_1\|_{\rm HS}^2
=\|WE_{\sqrt3}W^\dagger-E_{\sqrt3}\|_{\rm HS}^2
=\frac{784}{20736^2}>0.
}
\]

The `E_sqrt2` and `E_sqrt6` effects, proportional to identity, are individually invariant. This is an exact algebraic witness that the **particular localized record effects do not descend as fixed effects** under every such correlated matching. A transported family remains covariant; proving that the source alone supplies the transport is outstanding.

Note that `S` commutes with *diagonal* `U⊗U`, while the entangled join `P` is fixed by *correlated conjugate* `U⊗conj(U)`. Confusing these two different stabilizer actions would hide the obstruction.

## 6. Complete-positivity and gauge covariance do not enforce uniqueness

Let the depolarizing record effects be

\[
E_w^{\rm tr}=\frac{n_w}{31104}I_9,
\qquad (n_1,n_{\sqrt2},n_{\sqrt3},n_{\sqrt6})=(79,48,42,3).
\]

For any `0<=t<=1`,

\[
E_w^{(t)}=tE_w+(1-t)E_w^{\rm tr}
\]

are positive, local, covariantly transportable POVM effects. Their sum remains `(172/31104)I9`, and they commute with `P`. Hence they all admit the same exact physical `K_Bell=P` completion above (using the appropriate positive square roots).

Their *stand-alone* score on the original input is

\[
\boxed{f_t(\rho_A^{\rm initial})=
\alpha_{\rm cand}+t\frac{7(\sqrt3-1)}{124416}.}
\]

For `t=0`, the exact target is obtained only because the trace map explicitly erases physical state dependence. It is not a derivation. For `t=1`, one recovers the operator-frame localization and its discrepant source-state score.

| `t` | Stand-alone score on source input | Bell-first score, Bell label 0 |
|---:|---:|---:|
| `0` | `0.007297352564` | `0.006081127136` |
| `1/2` | `0.007317946200` | `0.006087991682` |
| `1` | `0.007338539836` | `0.006094856227` |

Even for a **fixed** full correlated-matching group, the trace-map member `t=0` is invariant, whereas the other members need transport. The original axiom provides no selection principle fixing `t` or the outcome-label rule.

## 7. Independent checks and scope limitations

`chain_beta_extension_v15.py` performs:

- A direct build of 172 record directions using v11/v14 code, their four weighted sectors, and normalized partial traces from all sectors. It does **not** input the rational `E_w` coefficients into that computation.
- A separate SymPy exact-matrix verification of those rational formulas, sector eigenvalues, trace multiplicities, positivity, and the scalar sum identity.
- Exact evaluation of the initial, Bell-selected and maximally mixed scores, plus a rotated-support state that accidentally gives the target value.
- Exact construction/normalization of the six-outcome Bell-first instrument, and a direct numerical Born-rule and no-signaling check.
- An exact `U=diag(i,-i,1)` gauge witness; 32 randomized frame-transport covariance trials.

The direct NumPy and exact formulas agree within `1.52e-18` entrywise. Kraus completeness is better than `1.2e-16`, and the no-signaling residual is below `2.8e-17` in the aligned input. Randomly transported probability residuals are below `2.9e-15`.

**Limitations:**

1. The v11 representation of the 172 record directions, the v14 CP map, and the local conditional expectation are mathematically explicit but **not uniquely prescribed** by Wissner-Gross's axiom. The source's Lean implementation certifies record enumeration and the weighted count but explicitly does **not** construct the 31,104-dimensional readout space. Replacing our chosen embedding could alter this result.
2. A conventional local CPTP instrument includes the failure branch excluded by the original axiom. It is an **operational extension**.
3. The scalar readout is a weighted mean of numerical outcomes, **not** itself the success probability of a physical photon measurement.
4. The gauge-covariant instrument is a **family transported between frames**; fixed measurement invariance fails in the explicit witness. A canonical source-derived relational connection remains missing.
5. No spacetime, momentum transfer, propagator, fermionic dispersion, loop integration, or renormalization scale has been derived. Thus no independent QED beta-function coefficient has been obtained.

## 8. Next mathematically decisive problem

A useful v16 should tackle **the intersection of physical locality, fixed gauge descent, and record-weight uniqueness**, not add another unconstrained dynamical model. In particular:

1. Characterize covariant CP maps from `R` to `End(Q1Q2)` that intertwine the actual `U⊗conj(U)` matching representation and admit the exact Bell branch `P3`.
2. Impose source-motivated channel compatibilities and study their Choi operators. Identify how many degrees of freedom remain and whether any principle selects one without injecting `alpha_cand`.
3. Evaluate every admissible map on the **actual three-pair state**, rather than only `I9/9`, and measure how their score depends on the allowed link support orientations.
4. Decide whether the source's normalized trace can be operationally identified with an electromagnetic coupling at all. Even success here would leave the separate physical QED beta-function problem open.

**Bottom line:** v15 succeeds in constructing an exactly local, CPTP, no-signaling, Bell-preserving instrument tied to the proposed record weights. It simultaneously demonstrates that its score is neither uniquely fixed by the axiom nor equal to the paper's candidate on the physical source input. The missing step is a uniquely derived gauge-compatible physical observation rule, and subsequently dynamical running with energy.