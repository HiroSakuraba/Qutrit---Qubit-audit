# Chain Theory → QED beta program, v21
## Source-level audit of read order, matching, frequencies, and sequential distillation

**Date:** 25 October 2026  
**Primary source:** A. D. Wissner-Gross, *A Parameter-Free Candidate for the Fine-Structure Constant from Entangled Qudits*, 23 October 2026, `ChainTheoryAlpha.pdf`.  
**Continuity:** v2–v20, especially v19 (support-overlap transfer) and v20 (freely chosen operational completions).  
**Status:** Exact algebraic results, a source-faithful double-swap calculation, and a conditional counterfamily. **Not** a no-go theorem for all possible added physics; **not** a derivation of a clock, propagator, or QED beta function.

## Abstract

We test whether the original matching and distillation axioms select an **event order**, **statistical weights for different acts**, or **relative operation frequencies**. The source says in Supplement S7 that the two central read insertions have **no time order**; Supplement S1 distinguishes the usual laboratory Born probability of an entanglement-swap outcome from the model's single selected branch, which has **no failure outcome**; Corollary S5 restricts reading to **one copy** of a link. The exchanges of reads and links, `X` and `Y`, are **involutive symmetries used to classify records**, not operators declared to occur randomly in time. The qubit anticommutator records attach the **sign of an algebraic exchange** to an operator. These are not intrinsic event timestamps or transition rates.

From the exact 172-record signature table, we reconstruct `X,Y`, prove the four character-projector ranks `(79,48,42,3)`, and derive the signed character traces `(82,70,-8)`. The weighted insertion remains the known source quantity. If one **adds** a quantum stochastic dynamics on this *auxiliary* record Hilbert space, there are independent nonnegative rates `γ_X,γ_Y` that preserve all source-sector weights and symmetries. We explicitly classify its `172²=29584` Liouville modes into multiplicities `(10318,6924,7836,4506)` with distinct freely choosable damping rates. This is a witness to underdetermination **even after** adopting a specific kind of dynamics, not a proposal that the source records themselves are physical states.

Independently, two Bell projections on distinct intermediate qudit copies in a chain of three maximally entangled links **commute**. Both histories give exactly `d⁻² |ω_d>_{16}|ω_d>_{23}|ω_d>_{45}` as an unnormalized pure state; the ordinary laboratory Born probability is `d⁻⁴`. The normalized source-style selected output is independent of the order. This proves that a nontrivial two-operation example contains no operational chronology in the *final selected state*. It does not rule out clocks or records added to a larger theory.

## 1. Precise source audit (what is supplied)

Primary-source references refer to the actual PDF's printed page numbers, not to claims from previous v2–v20 reports.

| Source | Exact role | What it does **not** say |
|---|---|---|
| Supplement S1, pp. 8–10: Corollaries S2–S5 | Distillation is the selected maximally entangled projection; a read operates on the one copy of its link; for anticommting qubit reads, interchanging the algebraic product supplies a sign. | It does not introduce a probability distribution over which two junctions are acted on next. |
| Supplement S1, p. 10: “The Born square” and “Distillation and postselection” | The paper deliberately uses Born squares as readout **weights**, not outcome frequencies within its ontology. Ordinary Bell swapping has a lab probability `1/d²`; the model keeps only the selected projection. | It does not equate these laboratory probabilities with internal time-step transition probabilities. |
| Supplement S7, p. 16: “The joint state” | The two central reads are insertions on the **joined state** and “the model has no time order.” The selected qutrit Bell join has ordinary lab probability `1/6`. | It does not say that read order is a causal evolution parameter. |
| Supplement S7, pp. 16–17: Corollary S6, Theorem S7 | `X` exchanges reads, `Y` exchanges links; they commute on the 172 selected records and yield exchange eigenvalue signs. The multiplicities are `79,48,42,3`. | It does not make `X,Y` stochastic transitions or fix their rates. |
| Supplement S11, pp. 22–23 | Lean checks the finite combinatorics and the computed quantity's arithmetic. | It explicitly does not build the complete physical readout measurement space or formalize the Photon Conjecture. |

**Important terminological separation.** The paper's **matchings** include local basis changes of a pair, of the form `U⊗Ū`, which fix a pair ket. Its `X,Y` are the **exchange-of-reads and exchange-of-links** operations used in the record weight law. We do not identify them with each other or with a physical time-evolution generator.

## 2. Source-defined record exchange algebra

Let the 172 mutually orthogonal source record directions define a real or complex auxiliary coordinate Hilbert space `K`, solely for bookkeeping. On this space, define `X,Y` by their record eigenvalues. From Theorem S7:

| Weight `sqrt(k)` | `X` | `Y` | Rank |
|---|---:|---:|---:|
| `1` | +1 | +1 | 79 |
| `sqrt(2)` | +1 | −1 | 48 |
| `sqrt(3)` | −1 | +1 | 42 |
| `sqrt(6)` | −1 | −1 | 3 |

The eigenprojector onto character `(x,y)∈{±1}²` is

\[
\boxed{P_{xy}=\tfrac14(I+xX)(I+yY).}\tag{21.1}
\]

The defining relations `X²=Y²=I`, `XY=YX`, make all four `P_xy` mutually orthogonal projectors with ranks listed above. Their traces imply

\[
\begin{split}
\operatorname{Tr} I&=172,\\
\operatorname{Tr} X&=79+48-42-3=82,\\
\operatorname{Tr} Y&=79-48+42-3=70,\\
\operatorname{Tr}(XY)&=79-48-42+3=-8.
\end{split}\tag{21.2}
\]

The paper's squared-weight operator on these record directions can be expressed as `K=(2I-X)(3I-Y)/2`; hence the positive weighted insertion restricted to them is

\[
\boxed{M=K^{1/2}=\sum_{x,y}3^{(1-x)/4}2^{(1-y)/4}P_{xy}.}\tag{21.3}
\]

Its trace is exactly

\[
\operatorname{Tr}M=79+48\sqrt2+42\sqrt3+3\sqrt6,
\qquad \alpha_{\mathrm{cand}}=\operatorname{Tr}M/31104.
\]

This expression is **a spectral classification of observables**. There is no term in Eq. (21.1)–(21.3) with units of inverse time, no forward transition kernel, and no positive probability distribution on histories. Both `X` and `Y` are involutions, so their labels cannot themselves distinguish an orientation of time.

### Does “the order of two reads is its sign” supply chronology?

Corollary S4 uses two anticommuting qubit observables `A,B` and the fact that `AB=−BA`. For any Pauli pair, `AB` and `BA` differ only by a global factor `−1`. As ordinary single-outcome Kraus operators, they therefore produce **the same physical channel** on any density operator:

\[
(AB)\rho(AB)^\dagger=(BA)\rho(BA)^\dagger.\tag{21.4}
\]

The distinction is a valid **algebraic record sign**, but a standalone unconditioned transformation cannot reveal which product order was applied. Coherent interference between two controlled orders *can* reveal a relative phase, but a control system, coherent alternative histories, and their operational interpretation would be extra structure. Equation (21.4) does not imply that commutators are never observable.

The symbolic test covers all three unordered Pauli pairs, with exact entries.

## 3. Conditional rate-family witness constructed from the exact source exchanges

To test whether the two involutions themselves could *select* a stochastic evolution, give the auxiliary 172-dimensional record space an **additional**, conventional quantum-channel interpretation. For arbitrary independent `γ_X,γ_Y≥0`, set

\[
\boxed{
\mathcal L_{\gamma_X,\gamma_Y}(\rho)
=\gamma_X(X\rho X-\rho)+\gamma_Y(Y\rho Y-\rho).
}\tag{21.5}
\]

Each term is a valid Gorini–Kossakowski–Sudarshan–Lindblad generator with jump operator `sqrt(γ_X) X` or `sqrt(γ_Y) Y`. Consequently `e^{tL}` is completely positive, trace preserving, and unital for `t≥0`. Since `X,Y` commute, it is covariant under both source exchange actions. It fixes each of the four `P_xy` and all weights in `M`, independently of its rates:

\[
\mathcal L_{\gamma_X,\gamma_Y}(M)=0.\tag{21.6}
\]

For a matrix unit `|r><s|`, define `ξ=x_rx_s` and `η=y_ry_s`. It is an eigenoperator of the superoperator with eigenvalue

\[
\boxed{\lambda_{rs}=\gamma_X(\xi-1)+\gamma_Y(\eta-1).}\tag{21.7}
\]

These `172²=29584` complex operator directions decompose as follows:

| Relation between record characters | Superoperator dimension | Eigenvalue |
|---|---:|---|
| Neither X nor Y parity differs | 10,318 | `0` |
| X parity differs only | 6,924 | `−2γ_X` |
| Y parity differs only | 7,836 | `−2γ_Y` |
| Both parities differ | 4,506 | `−2(γ_X+γ_Y)` |
| **Total** | **29,584** | |

The counts follow from ordered pairs of character blocks. For example,

\[
\begin{split}
D_{XY}&=2(79\cdot3+48\cdot42)=4506,\\
D_X&=2(79\cdot42+48\cdot3)=6924,\\
D_Y&=2(79\cdot48+42\cdot3)=7836,\\
D_0&=79^2+48^2+42^2+3^2=10318.
\end{split}\tag{21.8}
\]

The numerical suite explicitly checks positivity, completeness of the four-unitary convex-mixture representation, trace preservation, semigroup composition, and that `(γ_X,γ_Y)=(0.2,0.8)` differs from `(0.8,0.2)` while conserving `M`.

**Interpretive limit:** Eq. (21.5) is not a source-allowed repeated read. The source's `172` directions belong to an **operator-space record algebra**, not a physical 172-level system; it also explicitly rules out repeated acts on the same copy. The theorem proves only that the symmetry data, even granted an operational Hilbert-space interpretation, do **not choose the dimensionless ratio** `γ_X/γ_Y`. In particular, the original `sqrt(k)` record weights are **not** these `γ` rates.

### The canonical twirl still does not provide a clock

The uniform group twirl of `Z_2×Z_2`,

\[
\mathcal P(\rho)=\tfrac14\sum_{a,b\in\{0,1\}}X^aY^b\rho Y^bX^a,
\]

is an idempotent CP trace-preserving projection onto exchange-invariant operators. It has no intrinsic duration or nontrivial decay-rate ratio. Its mathematical uniform average reflects averaging over a **finite symmetry group**, not the assertion that random physical exchanges happen at equal frequencies or that the twirl is a law of motion.

## 4. Source-compatible selected Bell chain: two joins with no detectable order

Now work entirely with the source's allowed maximally entangled `d×d` pairs (Supplement S1, Eq. S1), and the **selected** branch of the same rank-one Bell projection.

Start with three links in six distinct qudit-copy slots:

\[
|\Psi_0\rangle=|\omega_d\rangle_{12}|\omega_d\rangle_{34}|\omega_d\rangle_{56},
\quad |\omega_d\rangle=\tfrac1{\sqrt d}\sum_{a=0}^{d-1}|aa\rangle.
\]

Let `J_1=P_d^{(23)}` and `J_2=P_d^{(45)}` be Bell projections at the two junctions. They act on **disjoint tensor slots**, so

\[
\boxed{[J_1,J_2]=0.}\tag{21.9}
\]

Applying both projections and evaluating the tensor contraction gives

\[
\boxed{
J_2J_1|\Psi_0\rangle=J_1J_2|\Psi_0\rangle
=\frac1{d^2}|\omega_d\rangle_{16}|\omega_d\rangle_{23}|\omega_d\rangle_{45}.
}\tag{21.10}
\]

This uses the one-swap result of Supplement S1 Eq. (S1) twice. The right-hand side is a normalized tensor product multiplied by `d⁻²`. Therefore its **ordinary laboratory Born probability** is

\[
\boxed{p_{\mathrm{lab}}=d^{-4}.}\tag{21.11}
\]

The source's postselection instead keeps the **normalized single selected configuration**, independent of any failure outcome:

\[
\boxed{
|\Psi_{\mathrm{selected}}\rangle
=|\omega_d\rangle_{16}|\omega_d\rangle_{23}|\omega_d\rangle_{45}.
}\tag{21.12}
\]

The two descriptions “join 23 first” and “join 45 first” are distinct externally labelled sequences, but the selected output ket is identical. No observable of that final state can identify their order. The calculations were checked with explicit state tensors for `d=2,3,4`, comparing the two orders directly and the formula `d⁻²`.

The argument does not exclude an *external* control apparatus recording the order. Nor does it prove that all possible compositions commute: projections sharing the same tensor slot can fail to commute, and the source's shared-qudit-as-virtual-copies definition must be respected. It demonstrates a concrete case in which the original distillation data are **insufficient to infer** an event order.

### Why the Born number does not define a source-internal waiting time

The laboratory probability `d⁻⁴` applies to a standard experimental **sequence of postselections** (where failures exist). In the source's ontology those failures are excluded and the selected result is the one admitted state. Neither a wait-time distribution nor a repetition rate is specified. Even if one voluntarily posits two event orders with probabilities `q` and `1−q`, the final selected state is identical for **every** `q∈[0,1]`; this extra scheduling parameter is invisible to the source's selected result.

## 5. Relationship to earlier results and limits of inference

- **v19:** An independently chosen support-overlap path admits nontrivial holonomy, but microscopic fluctuations change its continuum limit. That assumed a path order and shared ambient support comparison.
- **v20:** A source-compatible selected qutrit Bell branch admits continuum-many distinct conventional CPTP **failure completions**, parametrized by `s`. These failure states are outside the source's ontology.
- **v21:** Even staying at the source's exact exchange-character algebra and selecting only Bell-join outcomes, we find no stochastic rate or event order. The qubit “order sign” does not supply a clock, and exact commuting join histories are indistinguishable in the selected final state.

These results support a **conditional underdetermination statement**: the definitions tested so far contain no operational clock, propagator, or relative event rates. They do not prove that a future dynamical extension is impossible, nor a contradiction of the source's arithmetic.

## 6. What a next research advance would require

A productive new postulate should be formulated so its *independent content can be tested*, rather than chosen to yield an already-known QED coefficient. Candidate minimal ingredients are:

1. **Causal/event structure:** a rule fixing which joins are admissible together and which, if any, must precede others, with operationally identifiable history records.
2. **Physical normalization:** a consistent one-copy probability instrument or alternative semantics for selected-only outcomes, giving measurable rates without signaling.
3. **Dynamical scale:** an internally fixed conversion from counted events, graph geometry, or spectral modes to dimensionful physical time and energy; ideally with no free dimensionless rate ratios.
4. **Matter and gauge response:** a controlled scaling limit with propagating charged fermions and a conserved current, from which the photon polarization and renormalization could be derived rather than imported.

The sharply testable v22 question is whether the full, source-admissible **rewrite system on pair matchings** has a canonical causal partial order or confluent normal form, even after considering overlapping possible joins and virtual-copy sharing. If there is no canonical order, document an explicit critical-pair witness. If all such rewrites are confluent, distinguish mathematical confluence from a physical clock.

## 7. Independent verification

Script `chain_beta_extension_v21.py` constructs the exact record-character multiplicities, validates their spectral projectors and trace formulas, diagonalizes the induced Liouville-mode counts by combinatorial masks, tests random-unitary semigroups, symbolically checks every Pauli-pair order identity, and applies two literal Bell projections to six-qudit tensors of local dimensions 2, 3, and 4.

All new v21 assertions pass. The generated JSON includes each assertion, quantitative residuals, exact character tables, and ordinary Bell probabilities. Prior version regressions are recorded separately and are not prerequisites for the exact source-level theorem. No GPU, QED action, lattice fermion, physical time parameter, or fitting to the measured alpha was used.

**Bottom line:** the source produces an exact *order-sign algebra*, not a derivation of *chronological dynamics*. The distinction is now exhibited at the record-algebra level and by a two-join physical selected-state example. The first QED beta-function coefficient remains unestablished from the original axiom.