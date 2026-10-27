# Qutrit–Qubit / Chain Theory Audit

Independent mathematical, operational, and phenomenological checks of A. D. Wissner-Gross, *A Parameter-Free Candidate for the Fine-Structure Constant from Entangled Qudits* (23 October 2026). This is an **audit and reproducibility repository**, not an endorsement of a derived fundamental physical theory.

> **Current conclusion (v22, 27 October 2026):** The proposed 172-record arithmetic is mathematically reproducible and the qutrit–qubit symmetry construction has a genuine relationship to the Standard Model gauge algebra. An operational identification of its weighted count with electromagnetic coupling, a source-derived physical dynamics, an independent QED beta function, and predictions for the weak and strong couplings remain **unestablished**.

## Headline claim and status

The source proposes the following **candidate**, not an independently established physical coupling:

\[
\alpha_{\mathrm{candidate}}
=\frac{79+48\sqrt{2}+42\sqrt{3}+3\sqrt{6}}{31104}
\approx 0.00729735256366226,
\]

with inverse approximately **137.035999189566**.

The four sector multiplicities (79, 48, 42, 3) and the readout dimension \(31104=24\times36\times36\) come from the source's record enumeration. **The arithmetic does not establish that the result is electromagnetism.** The photon identification is conjectural. The weak/strong audit (v22) does not derive \(\alpha_2\) or \(\alpha_s\).

### Audit outcomes

| Theme | Main outcome | Versions |
| --- | --- | --- |
| QED running and lattice extensions | Ordinary QED beta coefficient can be obtained after importing additional continuum fermion/gauge dynamics; not derived from qudit axioms | v2–v4 |
| Postselection and readout frames | Bell-join branch has quantum-instrument and frame-identification ambiguities | v5–v10 |
| Full readout and physical measurement | 172 records can be embedded in the 31,104-dimensional operator space; CP/gauge-covariant physical measurement maps are not uniquely selected | v11–v17 |
| Geometry, transfer, and time | Virtual-copy entropies reconstruct link graph; support-overlap transfer and continuum limit depend on extra support-path regularity; original rules do not fix a physical clock | v18–v21 |
| Weak and strong interactions | Correct gauge Lie algebra and a quark-doublet representation analogy; nonabelian coupling magnitudes, running, full chiral spectrum, confinement not derived | v22 |

These are **conditional model analyses** and source-based checks, not validation against independent experimental predictions.

## Repository layout

- [`reports/`](reports/) — versioned mathematical and physics audits v2–v22.
- [`scripts/`](scripts/) — standalone Python verification programs v2–v22. Some later scripts import earlier scripts; keep the directory together.
- [`results/`](results/) — captured program outputs from the original audit runs. Re-running programs may write fresh outputs inside `scripts/`.
- [`notes/`](notes/) — short research summaries and limitations.
- [`sources/`](sources/) — primary-source citation and provenance. **The original PDF is not included in the repository at this time**; compare the reports with the source before citing them as primary results.

No historical v2–v21 bundle ZIPs are checked in: they largely duplicate the tracked text files and paper. Preserve history through Git instead.

## Reproduce

Python 3 with NumPy and SymPy is required for the later programs; some earlier diagnostics may also require SciPy. Install with:

```bash
python -m pip install -r requirements.txt
cd scripts
python chain_beta_extension_v21.py
python chain_weak_strong_audit_v22.py
```

To run the available research sequence from inside `scripts/`:

```bash
for v in $(seq 2 21); do
  echo "=== v$v ==="
  python "chain_beta_extension_v$v.py" || exit 1
done
python chain_weak_strong_audit_v22.py
```

This may take time and re-create result files next to the scripts. Historical results are in `results/`. Scripts are research diagnostics rather than a maintained production package; inspect the corresponding report for the assumptions and limitations of each test.

## Scientific interpretation and red lines

1. **Numerical agreement is not a physical derivation.** A weighted trace equal to a known constant needs a uniquely selected physical measurement and coupling definition.
2. **Symmetry is not dynamics.** (mathfrak{su}(3)oplusmathfrak{su}(2)oplusmathfrak{u}(1)) does not fix (g_3,g_2,g_Y), charged excitations, spacetime, or the beta functions.
3. **Conditional extensions are identified.** Physical CPTP completions, model Hamiltonians, continuum path rules, and Standard Model RG equations introduced for diagnostic purposes are not consequences of the original axiom.
4. **Weak/strong comparisons use externally measured inputs.** The results do not retroactively predict the PDG values.
5. **No independent verification of the first QED beta coefficient.** The coefficient (2/(3pi)) for one unit-charged Dirac fermion follows only in models importing corresponding fermion/gauge dynamics.

## Provenance

- Primary source: A. D. Wissner-Gross (2026), *A Parameter-Free Candidate for the Fine-Structure Constant from Entangled Qudits*. The author/paper is not affiliated with this audit repository by implication.
- The paper's supplementary material describes the 172 records, readout dimension and a Lean 4 checking artifact; its proof verifies an arithmetic identity, not a photon propagator.
- Reference experimental inputs and the provenance of any non-source hypotheses are described in [`reports/ChainTheory_weak_strong_v22_report.md`](reports/ChainTheory_weak_strong_v22_report.md).

## Next steps

Prioritize (a) a physical source-derived observable that uniquely selects the 172-record score, (b) scale-dependent propagation from permitted distillations without an imported clock, and (c) a nonabelian response and running-coupling prediction that does not insert the measured values by hand. A negative theorem about what the original axioms cannot determine is an equally valid research result.

## Integrity and limitations

Versioned reports and scripts are preserved as generated research records, including their own stated qualifications. Some earlier model conclusions are refined or corrected by subsequent investigations. Read the newer audits (especially v14–v22) before relying on an earlier intermediate claim. Independent peer review and verification against upstream source definitions remain desirable.
