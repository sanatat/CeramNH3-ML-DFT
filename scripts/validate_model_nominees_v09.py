#!/usr/bin/env python3
"""Audit geometry-informed OOF predictions for pre-registered v0.8 nominees.

This is retrospective adsorption validation, NOT confirmation of NH3 catalytic activity.
Never rank candidates by measured DFT labels during nomination.
"""
from pathlib import Path
import argparse, json
import numpy as np
import pandas as pd
from sklearn.metrics import mean_absolute_error, r2_score

NOMINEES=('Al3Fe9','Cu3Fe9','Fe9Ga3','Fe3Sn9')
ADS=('H*','N*','NH*')


def checked_csv(p, required):
    if not p.exists(): raise SystemExit(f'Missing {p}; run v0.7 and v0.6 first')
    df=pd.read_csv(p)
    missing=set(required)-set(df.columns)
    if missing:raise SystemExit(f'Missing columns in {p}: {sorted(missing)}')
    return df


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--predictions',default='outputs/local_nominations/site_predictions.csv')
    ap.add_argument('--sites',default='outputs/local_sites.csv')
    ap.add_argument('--outdir',default='outputs/nominee_validation')
    args=ap.parse_args(); o=Path(args.outdir);o.mkdir(parents=True,exist_ok=True)
    p=checked_csv(Path(args.predictions),['record_id','material','adsorbate','actual_eV','predicted_eV','provenance'])
    s=checked_csv(Path(args.sites),['record_id','site_type_proxy','nearest_metal_A','site_n_metal'])
    if p.record_id.duplicated().any() or s.record_id.duplicated().any():raise SystemExit('Duplicate record_id; cannot safely evaluate')
    p=p[p.material.isin(NOMINEES)].copy()
    if p.empty:raise SystemExit('No nominee rows; run scripts/nominate_local_rf.py first')
    p=p.merge(s[['record_id','site_type_proxy','nearest_metal_A','site_n_metal']],on='record_id',how='left',validate='one_to_one')
    if p.site_type_proxy.isna().any():raise SystemExit('Unmatched geometry rows')
    p['absolute_error_eV']=(p.actual_eV-p.predicted_eV).abs()
    p['signed_error_eV']=p.predicted_eV-p.actual_eV
    p['source_dataset']='MamunHighT2019 via CatBench; reference adsorption reaction energy'
    p.to_csv(o/'nominee_site_audit.csv',index=False)
    summary=[]; mat=[]
    for (material,ads),g in p.groupby(['material','adsorbate']):
        summary.append(dict(material=material,adsorbate=ads,n_sites=len(g),
            predicted_median_eV=g.predicted_eV.median(),dft_reference_median_eV=g.actual_eV.median(),
            predicted_min_eV=g.predicted_eV.min(),predicted_max_eV=g.predicted_eV.max(),
            dft_reference_min_eV=g.actual_eV.min(),dft_reference_max_eV=g.actual_eV.max(),
            site_level_MAE_eV=g.absolute_error_eV.mean(),
            median_absolute_error_eV=g.absolute_error_eV.median(),
            site_types=';'.join(sorted(g.site_type_proxy.astype(str).unique())),
            prediction_provenance=';'.join(sorted(g.provenance.astype(str).unique()))))
    t=pd.DataFrame(summary).sort_values(['material','adsorbate']);t.to_csv(o/'nominee_adsorbate_comparison.csv',index=False)
    for m in NOMINEES:
        g=p[p.material.eq(m)];have=set(g.adsorbate)
        for ad in ADS:
            if ad not in have:print(f'WARNING: {m} missing {ad}')
        mat.append(dict(material=m,sites=len(g),adsorbates=len(have),site_level_MAE_eV=float(g.absolute_error_eV.mean()) if len(g) else np.nan,
            R2=float(r2_score(g.actual_eV,g.predicted_eV)) if len(g)>1 and g.actual_eV.nunique()>1 else np.nan,
            geometry_data='post-DFT relaxed',N2_barrier_available=False,
            NH3_experimental_FE_available=False,verified_N2_electrochemical_match=False))
    pd.DataFrame(mat).to_csv(o/'nominee_material_audit.csv',index=False)
    print('=== NOMINEE ADSORPTION DFT CHECK (NOT NH3 PERFORMANCE) ===')
    print(t[['material','adsorbate','n_sites','predicted_median_eV','dft_reference_median_eV','site_level_MAE_eV']].to_string(index=False,float_format=lambda x:f'{x:.3f}'))
    print('\n=== MATERIAL-LEVEL DIAGNOSTICS ===')
    print(pd.DataFrame(mat)[['material','sites','site_level_MAE_eV','R2']].to_string(index=False,float_format=lambda x:f'{x:.3f}'))
    print('\nWARNING: Predicted adsorption energy and ground-truth DFT both originate in the same parent dataset. This is retrospective group-held-out ML evaluation, not experimental validation.')
    print('WARNING: The N-vs-H Pareto front is an unvalidated screening proxy, not proven NH3 yield ordering.')
    print('WARNING: Energies are reaction-reference energies; do not equate them uncritically to literal N binding free energies.')
    print('Wrote:',o)

if __name__=='__main__':main()
