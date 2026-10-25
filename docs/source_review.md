# Verified literature assessment (2026-10-23)

## Baselines
- Chen et al., *Catalysis Science & Technology* (2023), DOI 10.1039/D3CY00462G: DFT + microkinetic ordering Fe(111) > Fe(211) > Fe(210) > Fe(100) > Fe(110). These are model-derived activities, not directly measured ceramic-cell Faradaic efficiencies.
- Pan, Hansen & Vegge, *J. Mater. Chem. A* (2020), DOI 10.1039/D0TA08313E: RPBE calculations on vanadium oxynitrides. Reported values from abstract: onset potential magnitude 0.4 V-RHE (12.5% O), approx. 0.3 eV nitrogen adsorption/dissociation barrier (31.25% O), onset 0.6 V-RHE (31.25% O). **Do not label a barrier as adsorption energy**.
- Mamun et al., *Scientific Data* (2019), DOI 10.1038/s41597-019-0080-z, Catalysis-Hub publication MamunHighT2019: alloy-surface adsorption dataset with Quantum ESPRESSO source data. Check whether actual Fe/Ni surfaces, facets and N/NH/H adsorbates exist before training.
- Fe(110) adsorption has already been studied: *Surface Science* (2018), DOI 10.1016/j.susc.2017.09.002.
- Important novelty warning: a 2026 *Langmuir* article titled "AI-Driven Coverage-Dependent Kinetics for NH3 Synthesis on Fe(110)" already reports DFT + GNN screening and coverage-dependent kinetics. Do not claim an Fe(110) AI screening workflow is novel without distinguishing it.

## Curation and validity rules
1. Never mix direct N adsorption, N2 dissociation barriers, overall NH3 cell production rates, electrochemical onset potentials, and catalyst activity rankings as the same ML target.
2. A single database dataset may contain thousands of datapoints but **not thousands of independent publications**. Group validation by publication AND composition/slab/adsorption family; leave-one-paper-out is impossible for only one paper. Use composition-holdout for one-source datasets, with honest limitations.
3. Match XC, reference state (e.g. 1/2 N2 vs atomic N), facet, coverage, geometry, spin and computational settings. PBE bulk calculations cannot be mixed blindly with RPBE surface data.
4. The validated adsorption template remains empty until each row is traced to source equation, functional, system and reference.
5. Use Fe/Ni adsorption database as *pretraining / methodology*, not evidence of VN/Fe ceramic membrane performance.
