#!/usr/bin/env python3
"""Inspect whether CatBench JSON includes geometries; never invent sites."""
import json, argparse
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--input',default='raw_data/MamunHighT2019_adsorption.json');p.add_argument('--examples',type=int,default=3);a=p.parse_args()
path=Path(a.input)
with path.open() as f: obj=json.load(f)
print('Records:',len(obj))
for key in list(obj)[:a.examples]:
 rec=obj[key];print('\nRecord:',key);print('Top-level fields:',list(rec))
 for field,value in rec.items():
  print(' ',field,':',list(value)[:12] if isinstance(value,dict) else str(value)[:150])
print('\nCAUTION: references to geometry IDs are not coordinates. No local site descriptors generated.')
