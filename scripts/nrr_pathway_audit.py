#!/usr/bin/env python3
"""Audit coverage of NRR intermediates for the pre-registered ML nominees.

Does not transform adsorption reaction energies into CHE free energies, does not
infer transition states, and does not rank catalytic ammonia production.
"""
import argparse, json, re
from collections import Counter
from pathlib import Path
import pandas as pd

NOMINEES=['Al3Fe9','Cu3Fe9','Fe9Ga3','Fe3Sn9']
# Representative distal/alternating NRR intermediates, plus dissociative-route
# species. A label here does not imply an experimentally established mechanism.
SPECIES=['H*','N*','NH*','NH2*','NH3*','N2*','NNH*','NNH2*','NHNH*','NHNH2*','NH2NH2*']

def canonical_species(raw):
    # Handle trailing site identifiers, e.g. NH*_2
    return re.sub(r'_\d+$','',raw.strip()).replace(' ', '')

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--json',default='raw_data/MamunHighT2019_adsorption.json')
    ap.add_argument('--validation',default='outputs/nominee_validation/nominee_adsorbate_comparison.csv')
    ap.add_argument('--outdir',default='outputs/nrr_v10')
    args=ap.parse_args()
    path=Path(args.json)
    if not path.is_file():raise SystemExit(f'Missing {path}. Copy existing v0.5 raw_data JSON.')
    out=Path(args.outdir);out.mkdir(parents=True,exist_ok=True)
    with path.open(encoding='utf8') as f: data=json.load(f)
    rows=[]; total=Counter(); types=Counter()
    for k,v in data.items():
        if k.startswith('_') or not isinstance(v,dict) or '->' not in k:continue
        material=k.split('_',1)[0]
        species=canonical_species(k.split('->',1)[1])
        total[(material,species)]+=1;types[species]+=1
        if material in NOMINEES:
            rows.append({'material':material,'species':species,'record_id':k,'reference_reaction_energy_eV':v.get('ref_ads_eng'), 'reaction_label':k.split('_',1)[1]})
    pd.DataFrame(rows,columns=['material','species','record_id','reference_reaction_energy_eV','reaction_label']).to_csv(out/'raw_nominee_reactions.csv',index=False)
    coverage=[]
    for material in NOMINEES:
        for species in SPECIES:
            n=total[(material,species)]
            coverage.append({'material':material,'species':species,'n_reaction_records':n,'available_in_parent_dataset':bool(n), 'role':'adsorption reference only; not free energy of NRR step'})
    cv=pd.DataFrame(coverage);cv.to_csv(out/'intermediate_coverage.csv',index=False)
    pd.DataFrame([{'species':a,'count_all_materials':n} for a,n in types.most_common()]).to_csv(out/'parent_dataset_species.csv',index=False)
    wide=cv.pivot(index='material',columns='species',values='n_reaction_records').reindex(index=NOMINEES,columns=SPECIES).fillna(0).astype(int)
    print('=== PARENT-DATASET INTERMEDIATE COVERAGE (number of reaction entries) ===')
    print(wide.to_string())
    print('\nLegend: 0 = no such label in this JSON; counts do NOT supply CHE reaction free energies.')
    print('Known gaps: N2 adsorption, NNH hydrogenation, reaction barriers, Gibbs corrections, potential/pH, and HER energetics are required for NRR selectivity.')
    print('No activity ranking or limiting potential computed. Files in',out)

if __name__=='__main__':main()
