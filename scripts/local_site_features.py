#!/usr/bin/env python3
"""Explicit structure-linked local-site descriptors from ASE-readable structures.
Requires manifest with record_id,structure_file,adsorbate_indices (0-based ';'-delimited),
and material; full path must exist. No geometry inference from composition alone.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from ase.io import read
from ase.data import atomic_numbers

p=argparse.ArgumentParser();p.add_argument('--manifest',required=True);p.add_argument('--output',default='outputs/local_site_features.csv');p.add_argument('--cutoff',type=float,default=3.0);a=p.parse_args()
manifest=Path(a.manifest); df=pd.read_csv(manifest,dtype=str).fillna('')
required={'record_id','structure_file','adsorbate_indices','material'}
if not required.issubset(df.columns): raise SystemExit('Manifest missing: '+', '.join(required-set(df.columns)))
rows=[]
for r in df.to_dict('records'):
 if not r['adsorbate_indices'].strip(): raise ValueError('adsorbate_indices required: '+r['record_id'])
 file=Path(r['structure_file']); file=file if file.is_absolute() else manifest.parent/file
 atoms=read(str(file)); idx=[int(i) for i in r['adsorbate_indices'].split(';')]
 if any(i<0 or i>=len(atoms) for i in idx): raise ValueError('Invalid indices '+r['record_id'])
 metal_indices=[i for i,s in enumerate(atoms.get_chemical_symbols()) if s not in {'H','N','O','C','S'}]
 if not metal_indices: raise ValueError('No metal atoms '+r['record_id'])
 touched=[]
 for i in idx:
  for j in metal_indices:
   d=atoms.get_distance(i,j,mic=True)
   if d<=a.cutoff: touched.append((j,d))
 nearest=min((atoms.get_distance(i,j,mic=True),j) for i in idx for j in metal_indices)
 touched_unique=sorted(set(j for j,_ in touched))
 rows.append({'record_id':r['record_id'],'material':r['material'],'metal_neighbors_cutoff':len(touched_unique),'nearest_metal_distance_A':nearest[0], 'nearest_metal_Z':atomic_numbers[atoms[nearest[1]].symbol], 'neighbor_mean_Z':np.mean([atomic_numbers[atoms[j].symbol] for j in touched_unique]) if touched_unique else np.nan,'cutoff_A':a.cutoff,'structure_file':str(file)})
out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True);pd.DataFrame(rows).to_csv(out,index=False)
print(f'Saved {len(rows)} structure-linked site records to {out}')
print('WARNING: geometries must be matched to exact adsorption records; descriptors do not imply catalytic activity.')
