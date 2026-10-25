#!/usr/bin/env python3
"""Inspect published outputs; no dependency beyond Python standard library."""
import csv, re
from pathlib import Path
root = Path(__file__).resolve().parents[1]
with (root/'outputs/local_comparison/summary.csv').open() as f:
    for row in csv.DictReader(f):
        print(f"{row['representation']:<30} {row['model']:<4} CV MAE={float(row['CV_MAE_mean_eV']):.4f} eV")
qe = (root/'outputs/qe_Cu3Fe9_N2_Fe/n2_restart4.out').read_text()
assert 'bfgs converged' in qe and 'JOB DONE' in qe, 'Converged QE record missing'
energies = re.findall(r'^\s*!\s+total energy\s*=\s*([-\d.]+)\s+Ry', qe, flags=re.M)
print('Final converged slab+N2 energy:', energies[-1], 'Ry')
print((root/'outputs/qe_Cu3Fe9_N2_Fe/Cu3Fe9_N2_Fe_final_report.txt').read_text())
