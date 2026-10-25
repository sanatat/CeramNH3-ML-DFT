#!/usr/bin/env python3
"""Extract geometry-backed local adsorption site descriptors from CatBench JSON.
No DFT inference; skips records whose structures cannot be resolved.
"""
import argparse, io, json, re, math
from collections import Counter
from pathlib import Path
import numpy as np
import pandas as pd
from ase.io import read
from ase.data import covalent_radii, atomic_numbers

TARGETS={'H*','N*','NH*'}

def load_atoms(struct):
    if isinstance(struct,dict) and 'atoms_json' in struct:
        struct=struct['atoms_json']
    if isinstance(struct, str):
        return read(io.StringIO(struct), format='json')
    if isinstance(struct,dict):
        # ASE's JSON database may be a serialized string OR nested ASE encoding
        return read(io.StringIO(json.dumps(struct)), format='json')
    raise ValueError(f'Unsupported structure representation: {type(struct).__name__}')

def get_ads_slab(record, structures):
    raw=record.get('raw',{})
    candidates=[(k,v) for k,v in raw.items() if isinstance(v,dict) and float(v.get('stoi',0))>0]
    if len(candidates)!=1:
        raise ValueError(f'Expected one positive-stoichiometry product structure, found {len(candidates)}')
    key, value=candidates[0]
    obj=value.get('atoms_json')
    if obj is None:
        ref=value.get('ref')
        if ref is None or ref not in structures:
            raise KeyError('Missing structure ref in _structures')
        obj=structures[ref]
    return key,load_atoms(obj)

def fmt_material(name):
    counts=re.findall(r'([A-Z][a-z]?)(\d+)',name)
    if not counts: return None
    total=sum(int(n) for _,n in counts)
    return {el:int(n)/total for el,n in counts}

def extract_features(atoms, ads_idx, name, adsorbate):
    if any(i<0 or i>=len(atoms) for i in ads_idx):
        raise ValueError('Invalid adsorbate atom index')
    ads_indices=set(ads_idx)
    metal_idxs=[i for i,a in enumerate(atoms) if i not in ads_indices and a.symbol not in ('H','N','O','C','S')]
    if not metal_idxs: raise ValueError('No metal atoms identified')
    ads_atoms=[atoms[i] for i in ads_idx]
    # Chemically relevant anchor: N for N* and NH*, H for H*
    anchor_symbol='N' if adsorbate in ('N*','NH*') else 'H'
    anchors=[i for i in ads_idx if atoms[i].symbol==anchor_symbol]
    if len(anchors)!=1:
        raise ValueError(f'Cannot identify unique {anchor_symbol} anchor')
    idx=anchors[0]
    distances=atoms.get_distances(idx, metal_idxs, mic=True)
    ordered=sorted(zip(metal_idxs,distances),key=lambda z:z[1])
    nearest=float(ordered[0][1]); near=[(i,float(d)) for i,d in ordered if d<=nearest+0.45 and d<3.5]
    counts=Counter(atoms[i].symbol for i,_ in near)
    # Local first-shell definition is heuristic and MUST be sensitivity-tested.
    return {
      'ads_anchor':anchor_symbol,
      'site_n_metal':len(near),
      'site_n_Fe':counts['Fe'], 'site_n_Ni':counts['Ni'],
      'nearest_metal_A':nearest,
      'mean_site_distance_A':float(np.mean([d for _,d in near])),
      'site_Fe_fraction':counts['Fe']/len(near) if near else np.nan,
      'site_Ni_fraction':counts['Ni']/len(near) if near else np.nan,
      'site_type_proxy': 'top' if len(near)==1 else 'bridge' if len(near)==2 else 'multifold',
      'ads_structure_natoms':len(atoms),
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--input',default='raw_data/MamunHighT2019_adsorption.json')
    ap.add_argument('--output',default='outputs/local_sites.csv')
    ap.add_argument('--limit',type=int,default=0,help='0 = all matching records; positive for fast smoke test')
    a=ap.parse_args()
    p=Path(a.input)
    with p.open() as f: data=json.load(f)
    structures=data.get('_structures',{})
    print('Top-level entries:',len(data),'deduplicated geometries:',len(structures))
    if not structures:
        print('WARNING: _structures not present. Only inline atoms_json entries can be resolved.')
    rows=[]; errors=Counter(); examples={}; attempted=0
    for key, record in data.items():
        if key.startswith('_') or not isinstance(record,dict): continue
        material=key.split('_',1)[0]
        composition=fmt_material(material)
        if composition is None or not ({'Fe','Ni'}&composition.keys()):continue
        ads=re.sub(r'_\d+$','',key.split('->',1)[1].strip()) if '->' in key else ''
        if ads not in TARGETS:continue
        val=record.get('ref_ads_eng')
        if not isinstance(val,(float,int)) or not math.isfinite(val):continue
        attempted+=1
        try:
            prod,atoms=get_ads_slab(record,structures)
            idx=record.get('adsorbate_indices')
            if not isinstance(idx,list) or not idx: raise ValueError('Missing adsorbate_indices')
            props=extract_features(atoms,idx,material,ads)
            rows.append({'record_id':key,'material':material,'adsorbate':ads,'energy_eV':val,
              'Fe_fraction':composition.get('Fe',0),'Ni_fraction':composition.get('Ni',0),
              'contains_FeNi':('Fe'in composition and 'Ni'in composition),
              'raw_product_key':prod,**props})
        except Exception as exc:
            tag=f'{type(exc).__name__}: {str(exc)[:110]}'
            errors[tag]+=1
            if tag not in examples: examples[tag]=key
        if a.limit and attempted>=a.limit:break
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(rows).to_csv(out,index=False)
    print('Attempted:',attempted,'successful:',len(rows),'failed:',sum(errors.values()))
    print('Success by adsorbate:',dict(Counter(r['adsorbate'] for r in rows)))
    print('FeNi only:',sum(r['contains_FeNi'] for r in rows))
    print('Saved:',out)
    for error,n in errors.most_common(8):print(f'  ERROR {n}: {error} | example: {examples[error]}')
    if not rows: raise SystemExit('No geometry-backed records extracted; inspect _structures schema before ML')

if __name__=='__main__':main()
