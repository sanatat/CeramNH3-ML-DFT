#!/usr/bin/env python3
"""Inspect only; never auto-promote database energies as validated adsorption energies."""
from pathlib import Path
import csv
root=Path(__file__).resolve().parents[1]
for f in sorted((root/'data'/'cathub_raw').glob('*.csv')):
 with f.open(newline='') as fh: rows=list(csv.DictReader(fh))
 print(f.name, ':',len(rows),'records')
 for r in rows[:12]:print(' ',r['id'],r['chemicalComposition'],r['Equation'],r['reactionEnergy'])
print('No entries added automatically to vetted adsorption dataset.')
