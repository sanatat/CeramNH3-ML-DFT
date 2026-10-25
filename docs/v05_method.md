# v0.5 Research method and limitations

- Uses user-downloaded MamunHighT2019 data, without bundling it. The script reads the existing `outputs/feni_adsorption.csv` generated in v0.4.
- GroupKFold splits by `material`, preventing adsorption configurations of the same composition from appearing on both sides. An optional `--holdout-feni` excludes FeNi records **without evaluating them**; prior iterations already inspected them multiple times, so they are not a fresh holdout.
- Descriptors from ASE: periodic atomic number, atomic mass, covalent radius; weighted mean, weighted std, min, max, number of elements. Elemental fractions retained. No electronegativity or d-band values are fabricated.
- Local features require separate matched structure files and exact adsorption atom indices; JSON `raw` structure hashes/refs alone are NOT geometries. The local descriptors use a simple distance cutoff and are exploratory, not validated electronic features.
- Validate published adsorption reference-energy convention and DFT functional before merging additional datasets; predicted adsorption energies do not establish activity, kinetics, or operando stability for ceramic ammonia cells.
- Comparison on same 5 grouped folds supports feature ablation; for a publication-grade study implement nested group CV, uncertainty quantification, independent structures and site-aware validation, verify reference definitions and computational settings.
- With a few dozen FeNi records from only 3 unique compositions, there is no reliable independent FeNi benchmark; do not claim FeNi catalyst optimization.
