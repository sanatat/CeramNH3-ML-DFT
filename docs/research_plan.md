# Practical research plan

## Problem
Existing experimental VN–Fe ammonia work and FeNi/FeNi3 oxide-supported particle research motivate testing whether interface motifs transfer to nitrogen reduction. This is a hypothesis, not established transfer.

## Data-gathering order
1. Check supporting information for Pan et al. and related VN/VON nitrogen vacancy free energies and elementary steps. Capture surface, vacancy, coverage, functional, and reference states.
2. Query Catalysis-Hub for adsorbate-energy records with explicit slab formula, facet, site, method and reaction-energy definition. Keep provenance identifiers.
3. Inspect the OC20 bulk and adsorbate mappings before downloading anything large. Use N/N2/NH/H entries only with proper chemistry checks. OC20 may not be directly representative of high-temperature protonic ceramic cells.
4. Keep experimental NH3 rate and Faradaic efficiency as a separate table with pressure, temperature, gas composition, current, analytical blank and isotope controls where available.

## Candidate mechanistic tests
- N vacancy stability and refill on VN/VON; Fe–VN interface versus unmodified VN.
- H versus N-species adsorption and selectivity as *screening descriptors*, not as definitive pathways.
- Substitution/site effects for FeNi vs FeNi3 supported on relevant oxides, only after feasible slab structures and verified data.

## Validation
Hold out whole publications, not individual geometries from same article; beware related systems repeated across multiple publications. Keep RPBE and PBE separate. Baseline against global mean and composition-only models. Report uncertainty and out-of-domain cases.

## Computing
Preserve prior QE results but discontinue expensive 4-atom Fe relaxation as an immediate ML data source if impractical. For future target calculations prefer benchmarked smaller problems and consistent pseudopotentials before expensive slabs.
