#!/usr/bin/env python3
"""Extract actual CatBench Cu3Fe9 slab; generate unrelaxed N2/NNH guesses.
No QE input, DFT energies, or kinetics are implied by these geometries.
"""
import argparse, json, io
from pathlib import Path
import numpy as np
from ase.io import read, write
from ase.constraints import FixAtoms

CLEAN = '8706a1518968cf214cbd0683c3e1751e'
NREFS = ('9cb22dea48c76ae844e1733daa1d22b0', 'e99bb53b07cbd02c7cc8215ca8d9fbab')

def atoms_of(struct):
    if isinstance(struct, dict) and 'atoms_json' in struct:
        struct = struct['atoms_json']
    if isinstance(struct,str):
        return read(io.StringIO(struct), format='json')
    return read(io.StringIO(json.dumps(struct)), format='json')

def save_atoms(atoms, dest):
    write(str(dest)+'.traj', atoms)
    write(str(dest)+'.xyz', atoms, format='extxyz')

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--input',default='raw_data/MamunHighT2019_adsorption.json')
    ap.add_argument('--output',default='outputs/nrr_structures/Cu3Fe9')
    ap.add_argument('--distance',type=float,default=2.0,help='initial anchor-to-metal height (Angstrom), unrelaxed')
    args=ap.parse_args()
    data=json.loads(Path(args.input).read_text())
    structures=data.get('_structures',{})
    if CLEAN not in structures:raise SystemExit(f'Clean reference {CLEAN} missing in _structures')
    slab=atoms_of(structures[CLEAN]); out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    save_atoms(slab,out/'clean_as_deposited')
    for j,ref in enumerate(NREFS,1):
        if ref in structures:save_atoms(atoms_of(structures[ref]),out/f'Nstar_reference_{j}')
    sym=slab.get_chemical_symbols();print('Clean atoms:',len(slab),'composition:',{s:sym.count(s) for s in sorted(set(sym))})
    print('Cell lengths (Angstrom):',np.round(slab.cell.lengths(),3).tolist())
    print('PBC:',slab.pbc.tolist(), 'z range:', np.round([min(slab.positions[:,2]),max(slab.positions[:,2])],3).tolist())
    if len(slab)!=12:
        print('WARNING: expected 12 substrate atoms; review before use')
    if not slab.pbc[0] or not slab.pbc[1]:
        print('WARNING: slab lacks xy periodic boundary conditions; review cell')
    coords=slab.positions; zmax=coords[:,2].max()
    top=[i for i,s in enumerate(sym) if s in ('Fe','Cu') and zmax-coords[i,2] <= 1.3]
    if not top:raise SystemExit('Cannot identify exposed Fe/Cu surface atoms')
    for metal in ('Fe','Cu'):
        surface=[i for i in top if sym[i]==metal]
        if not surface:
            print('No top-layer',metal,'sites');continue
        anchor=surface[0]; x,y,z=coords[anchor]
        for kind in ('N2','NNH'):
            guess=slab.copy()
            # vertical end-on initial configurations; no claim of binding stability
            # place N nearest to metal, second N ~1.13 A above it
            guess += __import__('ase').Atoms('N2',positions=[[x,y,z+args.distance],[x,y,z+args.distance+1.13]])
            if kind=='NNH':
                guess += __import__('ase').Atoms('H',positions=[[x+0.94,y,z+args.distance+1.13+0.3]])
            # Preserve slab constraints, if present; otherwise no invented fixed-layer protocol
            save_atoms(guess,out/f'{kind}_on_{metal}_unrelaxed')
        print(f'Placed starting N2/NNH above {metal} index={anchor} z={z:.3f}')
    print('\nIMPORTANT: Geometries are UNRELAXED guesses, not optimized DFT results.')
    print('Inspect vacuum direction, spin, pseudopotentials, site orientation and original DFT settings before QE.')
    print('No barriers, NRR performance, or energies generated.')
    print('Saved to:',out)

if __name__=='__main__':main()
