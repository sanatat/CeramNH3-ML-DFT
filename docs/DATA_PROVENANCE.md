# Data provenance, reproducibility and limits

The 2026-10-25 project export provided the scripts and CSV outputs included here. Raw `MamunHighT2019_adsorption.json` (>90 MB in the original workspace) and original atomistic structure cache were intentionally not included; the derived tables cannot substitute for the upstream raw archive in a full audit.

- `outputs/feni_adsorption.csv`: derived adsorption reaction records with material and adsorbate labels; inspect the original records, source conventions and licenses before external redistribution of additional data.
- `outputs/local_sites.csv`: features extracted from relaxed DFT adsorption-site structures. They cannot be assumed available for an uncalculated candidate surface.
- `outputs/local_comparison/`: material-grouped 5-fold CV outputs, not external experimental validation.
- `outputs/local_nominations/`: retrospective prediction and multiobjective heuristic proxies, not catalytic rates.
- `outputs/qe_Cu3Fe9_clean/clean_scf.out`: converged single-point SCF of the as-deposited clean slab, not optimized bare slab.
- `outputs/qe_Cu3Fe9_N2_Fe/n2_restart4.out`: last segment of a multi-run relaxation using saved checkpoint states; the full initial and intermediate QE .out segments are not provided in this export.
- QE pseudo UPF files, checkpoints, and large caches are absent; a from-scratch recreation needs pseudopotential sourcing and the original slab/setup and should test numerical convergence.
- `data/` contains a combination of template/example and literature-review material and is not universally a curated empirical dataset.

No N₂ gas reference or adsorption free energies are provided; no kinetics, HER/NRR selectivity or ammonia synthesis performance conclusion is warranted.
