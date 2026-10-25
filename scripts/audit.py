#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
root=Path(__file__).resolve().parents[1]
papers=pd.read_csv(root/'data/papers.csv').fillna('')
obs=pd.read_csv(root/'data/literature_observations.csv').fillna('')
ads=pd.read_csv(root/'data/adsorption_template.csv').fillna('')
issues=[]
for name,df in [('papers',papers),('literature_observations',obs),('adsorption_template',ads)]:
    print(f'{name}: {len(df)} rows')
    if df.duplicated().any(): issues.append(f'duplicate rows in {name}')
for _,r in obs.iterrows():
    if r.paper_id not in set(papers.paper_id): issues.append(f'unknown paper for {r.record_id}')
    if not r.source_location: issues.append(f'no provenance location for {r.record_id}')
    if not r.unit: issues.append(f'no unit for {r.record_id}')
for _,r in ads.iterrows():
    if r.paper_id not in set(papers.paper_id): issues.append(f'unknown paper for {r.record_id}')
    for key in ['xc_functional','energy_definition','source_location','facet','adsorbate']:
        if not str(r[key]).strip(): issues.append(f'{r.record_id}: missing {key}')
print('Audit:', 'PASS' if not issues else 'ISSUES')
for problem in issues:print('-',problem)
raise SystemExit(bool(issues))
