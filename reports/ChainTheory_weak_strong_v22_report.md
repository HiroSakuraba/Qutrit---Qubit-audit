# Chain Theory v22: Can its qutrit–qubit construction derive the weak and strong couplings?

**Critical source-based and independent theoretical audit — 27 October 2026**

Source: A. D. Wissner-Gross, *A Parameter-Free Candidate for the Fine-Structure Constant from Entangled Qudits* (23 October 2026), `ChainTheoryAlpha.pdf`. Numerical benchmark sources: Particle Data Group, *2025 Review of Particle Physics*, sections 9 (QCD) and 10 (electroweak), versions reviewed August/November 2025 (electroweak document stamped April 2026). Linked below. All additional quantum-field-theory assumptions are identified as such.

## Executive verdict

The original construction correctly exhibits the Lie algebra `su(3) ⊕ su(2) ⊕ u(1)` and the representation `(3,2)_{1/6}` of the Standard Model. Those are nontrivial, legitimate algebraic correspondences. **It neither calculates nor fixes the physical weak coupling `g_2`, the QCD coupling `g_3`, or their evolution with energy.** More strongly, its **central-charge selection mechanism is intrinsically abelian**; the nonabelian simple ideals `su(2)` and `su(3)` have trivial centers, so the exact same charge-counting prescription does not transfer directly to the weak and strong interactions.

As a proposal *for the values of the two couplings*, this construction presently **fails to supply a derivation**. That is not a proof that a future dynamical extension is impossible or that the 172-record arithmetic is incorrect.

## 1. What the source itself states

- Letter, Eq. (3) and following: qutrit/qubit link has the identification algebra `su(3) ⊕ su(2) ⊕ u(1)`. Its one-dimensional center is spanned by `Q = i diag(3,3;−2,−2,−2)`, proportional to hypercharge `6iY`, which acts as a scalar `i` on the pair ket.
- Supplemental S3: the link carries the `(3,2)_{1/6}` representation, whose quantum numbers coincide with one left-handed quark doublet. The paper **does not establish that this link is a left-handed dynamical fermion**.
- Letter, Eqs. (4)–(7), Supplemental S7/S8: a unique chosen readout produces `alpha_cand = (79+48 sqrt(2)+42 sqrt(3)+3 sqrt(6))/31104 = 0.007297352563662262…`.
- Letter Conjecture 1 states that the link is an elementary photon link; only under this conjecture is the record count called `alpha(0)`.
- The letter explicitly says that translating its **hypercharge** into **electric charge**, and deriving spacetime and helicity, lie beyond the paper.
- Supplemental S8, Fig. S2 explicitly states that the model has **no charged loop, vertex factor, loop momentum, expansion in a coupling, momentum-dependent propagator or derived massless pole**. It describes a *correspondence of roles*, not equality of Feynman diagrams.
- No formula for `g_2`, `g_3`, `alpha_2`, `alpha_3`, `G_F`, the weak-mixing angle, QCD confinement, or corresponding beta functions is derived in the inspected source. The data below should not be called the paper's predictions.

## 2. Rigorous structural obstruction: there are three free couplings

For a connected Lie group whose Lie algebra is `su(3) ⊕ su(2) ⊕ u(1)`, each simple ideal has its own invariant nondegenerate quadratic form, and the abelian ideal has an independent quadratic form. Cross-terms between distinct ideals vanish by invariance. The most general positive gauge kinetic form is consequently parameterized by **three independent positive coefficients**, schematically

`L_gauge = -F_3^2/(4 g_3^2) - F_2^2/(4 g_2^2) - F_Y^2/(4 g_Y^2)`.

The factors `g_3`, `g_2`, `g_Y` can be varied without altering the gauge group or its representations. Quotienting the global group by its finite `Z_6` center does not eliminate this freedom. Merely obtaining the *correct gauge group* therefore **cannot** determine the physical coupling strengths. This is an exact algebraic point, conditional only on a conventional local Yang–Mills interpretation, which the source has not yet constructed.

### The center is different for the weak and strong factors

`Z(su(3)) = Z(su(2)) = {0}`, whereas `Z(s(u(3) ⊕ u(2))) ≅ u(1)`.

The paper's central charge `Q` selects that last, abelian direction. There is no analogous nonzero central `SU(2)` or `SU(3)` generator whose conserved scalar charge could simply be substituted into its `alpha` trace formula. Nonabelian gauge couplings would require a different observable, such as a properly derived nonabelian two-point function or physical Wilson-loop/string-tension construction.

## 3. Weak interaction audit

The source supplies a fundamental `SU(2)` representation, but it does not select **left-handed** weak-isospin couplings rather than a vectorlike or right-handed assignment. The pair ket's diagonal stabilizer is not a spontaneous-symmetry-breaking mechanism; its mathematical stabilizer does not generate Higgs physics, massive W and Z bosons, neutrinos, or weak charged-current decays.

At the `Z` mass, using the PDG 2025/2026 `MS-bar` weak mixing angle `sin^2(theta_W)(M_Z)=0.23122` and a representative five-flavour `MS-bar` electromagnetic input `1/alpha_em(M_Z)=127.930`, the standard electroweak matching formula gives

`alpha_2(M_Z) = alpha_em(M_Z)/sin^2(theta_W)(M_Z) = 0.0338066551281`, `g_2(M_Z)=sqrt(4pi alpha_2)=0.651787509524`.

**Both external inputs and the electroweak mixing relation were taken from the established Standard Model.** Nothing in the qudit record count derives these numbers or fixes the angle. The PDG explicitly distinguishes `MS-bar`, effective, and on-shell weak angles; conflating them would produce an artificial mismatch. The low-energy Fermi constant `G_F`, with units of inverse mass squared, is *not itself the same dimensionless quantity* as `alpha_2`.

## 4. Strong interaction audit

The source supplies the `SU(3)` algebra, hence its eight generators and their commutators, but does not supply color gluon fields, a gauge connection, gauge-field dynamics, quark masses, a renormalization scale, or nonperturbative confinement.

PDG 2025: `alpha_s(M_Z)=0.1180 ± 0.0009`, so `g_3(M_Z)=sqrt(4pi alpha_s)=1.21771578478`. This is a measured input **at 91 GeV**, not a low-energy constant analogous to `alpha_em(0)`.

One-loop QCD running, **after adding a conventional QCD Lagrangian and six quark flavors at high energies**, is controlled by `b_3=-7` (the coefficient convention `d alpha_3^-1 / d ln mu = -b_3/(2pi)`); a pure gauge `SU(3)` theory without quarks instead has `b_3=-11`. Thus knowing `SU(3)` alone does not determine the physical QCD beta coefficient, because it depends on the dynamical matter spectrum. The magnitude at a reference scale and the confinement scale remain separate input information unless new physics fixes them.

## 5. Numerical three-coupling comparison at a shared scale

All values below at `M_Z≈91.1880 GeV`. The electroweak numbers use the above `MS-bar` inputs at leading matching order, with small thresholds and higher-order matching effects omitted. `alpha_1=(5/3)alpha_Y` is **GUT-normalized hypercharge**, not `alpha_em`.

| Coupling | Gauge role | Reference number | How obtained |
|---|---|---:|---|
| `alpha_em(M_Z)` | photon (electroweak mixture) | 0.007816774799 | External electroweak input |
| `alpha_1(M_Z)` | normalized hypercharge | 0.016946275915 | Standard Model mixing and `5/3` normalization |
| `alpha_2(M_Z)` | weak isospin | 0.033806655128 | Standard Model mixing, external `sin²θ_W` |
| `alpha_3(M_Z)` | strong color | 0.118000000000 | PDG QCD world average |
| `alpha_cand(0)` | proposed *on-shell* photon candidate | 0.007297352564 | Source's finite record count, **not at the same scale** |

The ratios `alpha_2/alpha_em≈4.32489` and `alpha_3/alpha_em≈15.09574` are data-derived. The counts of nonabelian generators (`3` and `8`) do not imply those numbers. Even an explicitly unsupported guess `alpha_weak = 3 alpha_cand`, `alpha_strong=8 alpha_cand` would give `0.02189` and `0.05838`, far from the benchmark **and with inconsistent scales**, so this is a failed guess, not a physical falsification of the source.

## 6. The SU(5) normalization test: a real mathematical link, not an automatic prediction

The original `S(U(3)×U(2))` is indeed a standard subgroup of `SU(5)`. With `Y=diag(1/2,1/2,-1/3,-1/3,-1/3)` and `T3=diag(1/2,-1/2,0,0,0)`, the exact traces are

`Tr(Y²)=5/6`, `Tr(T3²)=1/2`.

Canonical `SU(5)` generator normalization therefore yields `T_Y=sqrt(3/5)Y` and `g_1=sqrt(5/3)g_Y`. **If** one imposes the additional hypothesis of a common `SU(5)` gauge coupling at a unification scale, one derives `sin²(theta_W)=3/8` **at that scale**. This is not a prediction at the Z mass, and simply identifying the subgroup does not impose a common `g_5`.

The usual one-loop Standard Model RG equations have `b=(41/10,-19/6,-7)` for `(alpha_1,alpha_2,alpha_3)`. Using the external Z-scale values, *pairwise* coupling crossings occur approximately at:

- `alpha_1=alpha_2`: `1.03×10^13 GeV`;
- `alpha_1=alpha_3`: `2.42×10^14 GeV`;
- `alpha_2=alpha_3`: `9.63×10^16 GeV`.

They do not meet at a single point in this minimal one-loop extrapolation. This illustrates why a *further* dynamical unified theory and matter spectrum would be needed. It is **not a refutation of the qudit paper**, which has no high-scale unification claim; nor of all unified extensions, which can add thresholds and matter.

## 7. Conditional consistency test: chiral gauge anomalies

The paper observes that its link carries the **representation** `(3,2)_{1/6}`, exactly that of the left-handed quark doublet. Suppose one attempts to declare it the *complete* chiral fermion content. Then its gauge anomalies are nonzero:

| Chiral-anomaly contribution of a single LH `(3,2)_{1/6}` | Value |
|---|---:|
| `SU(3)^3` in fundamental-unit normalization | `2` |
| `SU(3)^2 U(1)_Y` | `1/6` |
| `SU(2)^2 U(1)_Y` | `1/4` |
| `U(1)_Y^3` | `1/36` |
| gravity–`U(1)_Y` | `1` |
| `SU(2)` global (Witten) anomaly | present: **3** LH doublets, an odd number |

An anomaly-free one-generation Standard Model adds `(bar 3,1)_{-2/3}`, `(bar 3,1)_{1/3}`, `(1,2)_{-1/2}`, `(1,1)_{+1}` (all as LH Weyl fields). Then the five listed local anomalies vanish, and there are **four** `SU(2)` doublets, so the global obstruction vanishes.

**Caveat:** This does *not* demonstrate that Chain Theory has gauge anomalies. The original source identifies one quantum-information representation and conjectures a **photon**, not a complete spectrum of left-handed Weyl fermions. This is an actionable consistency condition for any attempted chiral Standard Model extension. The presence of the quark-doublet quantum numbers by itself cannot be advertised as a full reconstruction of matter content.

## 8. Verdict by claim

| Claim | Assessment |
|---|---|
| The link admits the correct Standard Model Lie algebra | **Mathematically supported** |
| The link has a genuine `(3,2)_{1/6}` representation | **Mathematically supported** |
| Hypercharge-normalization is compatible with the SU(5) subgroup | **Supported, conditional on convention** |
| The source predicts `g_2` or `alpha_2` | **No** |
| The source predicts `g_3` or `alpha_s` | **No** |
| Nonabelian coupling follows from 3/8 generator multiplicities | **No** |
| Weak chirality, electroweak mixing, W/Z masses | **Not supplied** |
| QCD beta coefficient, gluon propagator, confinement | **Not supplied** |
| The 172-record formula is invalidated by these tests | **No** |
| The theory has explained weak and strong interactions | **Not established** |

The source's mathematical symmetry identification does **not** fall apart under scrutiny. Its interpretation as a mechanism that explains all three force strengths currently does. The numerically successful electromagnetic readout remains a conjectural physical identification, and extending it to the two nonabelian sectors lacks a specified readout, a physical state/measurement rule, and dynamics.

## 9. Highest-value next tests

1. **Nonabelian observable definition:** From the original allowed distillations only, define an operational `SU(2)`/`SU(3)` gauge-covariant correlator. Show it reduces to a gauge-invariant measurable quantity (e.g. a suitable Wilson loop or physical scattering amplitude), not an arbitrary count.
2. **Matter and chirality:** Derive the complete anomaly-free chiral multiplet content, or prove that the link should *not* be interpreted as matter and identify what it represents instead.
3. **Physical scale and beta coefficients:** Derive a momentum/scale transformation and both nonabelian beta coefficients (including the matter contributions) *without* supplying the Standard Model field content and dynamics as axioms.
4. **Predictions before comparisons:** Publish at least two independent, numerical predictions with fixed schemes/scales and no choices tuned to PDG values: e.g. `sin²θ_W(M_Z)` and `alpha_s(M_Z)` with uncertainties, or one robust ratio at a common scale.

**Stop criterion:** If the only way to obtain numbers close to weak/strong couplings is to choose a separate normalization, matching scale, state, or record weighting for each sector, the extension is numerological rather than predictive.

## Primary references

- A. D. Wissner-Gross, *A Parameter-Free Candidate for the Fine-Structure Constant from Entangled Qudits* (Oct 23 2026), PDF attached in this conversation; sections cited above.
- Particle Data Group, *Quantum Chromodynamics*, 2025, https://pdg.lbl.gov/2025/reviews/rpp2025-rev-qcd.pdf (world average `alpha_s=0.1180±0.0009`, section 9.4).
- Particle Data Group, *Electroweak Model and Constraints on New Physics*, 2025 updated Apr 2026, https://pdg.lbl.gov/2025/reviews/rpp2025-rev-standard-model.pdf (sections 10.2, 10.3, 10.4; `MS-bar sin²(theta_W)=0.23122±0.00006`; `MS-bar alpha^-1(M_Z)=127.930±0.008` five-flavor reference).

### Reproducibility

Run `python chain_weak_strong_audit_v22.py`. The script uses SymPy for exact generator traces and chiral anomaly sums, and Python arithmetic for benchmark couplings and illustrative one-loop RG crossings. The companion JSON records provenance, assumptions, and results. Nothing in the script treats empirical weak/strong coupling inputs as outputs of Chain Theory.