# CeramNH3-ML-DFT

**Machine-learning-guided screening of potential ammonia-synthesis catalysts, supported by adsorption-energy analysis and a Quantum ESPRESSO DFT case study of N₂ on Fe-based Cu₃Fe₉.**

The long-term objective is to identify and evaluate promising materials for **N₂-to-NH₃ conversion**. This repository documents an **early-stage research workflow**, not a completed ammonia-synthesis reaction mechanism or experimentally validated catalyst. It links adsorption data and interpretable computational screening to an initial first-principles geometry test on one model-selected candidate.

> **Research status:** preliminary computational results, **not** a verified ammonia-synthesis catalyst discovery or proof of nitrogen-reduction activity.

## Project goal and current scope

The intended research pipeline is **adsorption-data curation → model evaluation → heuristic catalyst nomination → targeted DFT checks → future reaction-pathway and experimental validation**. Only the first four parts are represented here, and the available DFT check assesses an optimized N₂ adsorption geometry—not ammonia formation, reaction kinetics, or catalyst ranking.

## What this repository contains

- Python pipelines for adsorption-record extraction, descriptor analysis, material-grouped cross-validation, relaxed-geometry site features, and retrospective model-based nominations.
- Real tabular outputs exported from the project (not illustrative ML scores).
- Quantum ESPRESSO inputs and terminal output for the **clean Cu₃Fe₉ slab SCF** and **converged end-on N₂-on-Fe BFGS optimization**, together with final XYZ and a small result report.
- A coverage audit identifying missing NRR intermediates and unrelaxed starting geometries for further calculations.

## ML screening for ammonia-synthesis research

A curated subset derived from the **MamunHighT2019 adsorption dataset** contains **1,309 Fe- or Ni-containing reaction records spanning 184 material compositions**, covering H*, N*, and NH* reference adsorption/reaction energies. FeNi compositions were excluded from the 5-fold material-grouped model development (1,285 rows; 181 materials) and treated as a separate exploratory holdout.

| Model / representation | Grouped CV MAE (eV; mean of folds) |
|:--|--:|
| Random forest; composition and elemental physics | 0.4894 |
| SVR; composition and elemental physics | 0.5053 |
| Random forest; **relaxed adsorption-site geometry added** | **0.2706** |
| SVR; relaxed adsorption-site geometry added | 0.3414 |

The geometry-informed random forest has an approximate **out-of-fold R² = 0.8089** on the material-grouped development split; see `outputs/local_comparison/cv_predictions.csv` and `fold_results.csv`. Reported MAE is in eV under the source dataset's energy-reference conventions. **The geometry model uses relaxed adsorbate configurations derived from DFT: these inputs are not known before an adsorption geometry calculation. Consequently, the 0.2706 eV result is retrospective and must not be portrayed as prospective, DFT-free catalyst screening.** Cross-validation is grouped by material, not an independent literature/external experimental validation.

### Retrospective nominees

The relaxed-site, N*/H* adsorption-energy proxy produced the following examples: **Al₃Fe₉, Cu₃Fe₉, Fe₉Ga₃, Fe₃Sn₉**. These are **heuristic Pareto nominees**, not predicted turnover rates, electrochemical Faradaic efficiency, reaction barriers, or ammonia yield. See `outputs/local_nominations/pareto_nominees.csv` and `outputs/nominee_validation/` for additional post hoc checks. Note that source `ref_ads_eng` values can be reaction-referenced (including dissociative conventions); do not directly label every datum a molecular N₂ adsorption energy.

## DFT case study for a candidate ammonia catalyst: Cu₃Fe₉ + N₂

**Code:** Quantum ESPRESSO 7.5 (`pw.x`), spin-polarized PBE-family pseudopotentials, slab model with 9 Fe + 3 Cu atoms; two nitrogen atoms adsorbed at an Fe site. Calculation parameters below describe the supplied run, **not a demonstrated basis-set/k-point/pseudopotential convergence benchmark**.

| Quantity | Recorded value |
|:--|--:|
| Clean slab SCF total energy | −1168.47185003 Ry |
| Converged N₂/Fe/Cu₃Fe₉ slab total energy | **−1210.00716143 Ry** |
| End-on N–N distance | **1.125855 Å** |
| Fe–N distance | **1.896271 Å** |
| Final reported total force | 0.001495 Ry/Bohr |
| Final reported BFGS gradient error | 0.000750 Ry/Bohr |
| BFGS message | `bfgs converged in 5 scf cycles and 4 bfgs steps` |

The reported relaxation was reached through controlled stop/restart segments following an interrupted earlier run. The final artifact is `outputs/qe_Cu3Fe9_N2_Fe/n2_restart4.out`; the final coordinates are in `Cu3Fe9_N2_Fe_final.xyz`. The clean-slab SCF was **not** a clean-slab ionic relaxation. The N–N bond shows **no appreciable elongation relative to the 1.13 Å starting guess** in the tested end-on adsorption geometry. This does not exclude alternative adsorption sites or reaction routes.

**An N₂ adsorption energy is NOT calculated in this repository.** It requires a consistently referenced isolated N₂ calculation, a well-controlled clean slab, and an appropriate comparison of geometries/pseudopotentials and energy conventions:

`E_ads(N₂) = E(slab+N₂) - E(slab) - E(N₂ gas)`.

Raw total energies cannot be compared directly with the MamunHighT2019 reference energies. Formation of NNH*, other NRR intermediates, transition barriers, selectivity against hydrogen evolution, and experimental NH₃ production remain **unvalidated**. The NRR source audit in `outputs/nrr_v10/` shows missing N₂*, NNH*, and several downstream species from the parent adsorption subset.

## Reproduce and inspect

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

# Recompute grouped CV on the included small derived tables; may take CPU time
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python scripts/compare_local_models.py

# Inspect included saved results, without requiring Quantum ESPRESSO
python scripts/report_results.py
```

The scripts expect to be run from the repository root. QE inputs reference pseudopotential files under a **local `pseudo/` directory** (`Fe.jry.pbe.UPF`, `Cu.pbe-dn-kjpaw_psl.0.2.UPF`, `N-PBE.upf`), which are **not distributed** here. Some upstream extraction scripts require the much larger original dataset and structures (`raw_data/MamunHighT2019_adsorption.json`), also omitted from this archive. Read `docs/DATA_PROVENANCE.md` before rerunning or republishing data.

The final QE run's `n2_restart4.in` is a **historical restart input** and requires checkpoint files to run unchanged. `final_geometry_scf_template.in` is a **separate from-scratch SCF verification template**, not the input used to generate the reported relaxation energy; it still requires compatible pseudopotentials.

## Structure

- `scripts/`: executable analyses and original research scripts.
- `data/`: source-review notes, literature and example templates. Some rows are templates; **do not mistake them for verified published measurements**.
- `outputs/local_comparison/`: out-of-group prediction and fold results.
- `outputs/local_nominations/`, `outputs/nominee_validation/`: exploratory nominees and audits.
- `outputs/qe_Cu3Fe9_clean/` and `outputs/qe_Cu3Fe9_N2_Fe/`: QE evidence and geometry.
- `outputs/nrr_v10/`: intermediate-coverage audit.
- `docs/`: methods, limitations and provenance.

## Current conclusions and scientific limitations

The main deliverable is a traceable **candidate-screening and DFT-feasibility workflow for ammonia-synthesis research**. Successful numerical optimization of N₂ on the tested Cu₃Fe₉ site demonstrates a stable local geometry under the reported computational setup; it does **not** establish favorable molecular N₂ adsorption energy, catalytic turnover, viable NH₃ reaction steps, or improvement over pure Fe.

## Responsible interpretation

This is an **educational/research portfolio prototype** with traceable calculations and explicit uncertainties. No measured catalytic activity, electrochemical efficiency, N₂ bond dissociation barrier, NH₃ production rate, or independent experimental confirmation is claimed.
