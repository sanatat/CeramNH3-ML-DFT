#!/usr/bin/env python3
"""Compare composition and geometry-informed regressors with group-safe folds.
CAUTION: Relaxed adsorption geometries are post-DFT descriptors and cannot be
used as inputs for pre-DFT screening without an independently accessible geometry.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold
from sklearn.metrics import mean_absolute_error
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.svm import SVR
import sys

# Reuse exactly the descriptors from v0.5
sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_descriptor_study import prepare, parse_material

SITE_NUMERIC = [
    'site_n_metal', 'site_n_Fe', 'site_n_Ni', 'nearest_metal_A',
    'mean_site_distance_A', 'site_Fe_fraction', 'site_Ni_fraction',
]

def model_factory(name):
    if name == 'RF':
        return RandomForestRegressor(n_estimators=180, min_samples_leaf=3,
                                     random_state=42, n_jobs=1)
    if name == 'SVR':
        return make_pipeline(StandardScaler(),SVR(kernel='rbf',C=10,epsilon=0.1))
    raise ValueError(name)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--composition',default='outputs/feni_adsorption.csv')
    ap.add_argument('--sites',default='outputs/local_sites.csv')
    ap.add_argument('--outdir',default='outputs/local_comparison')
    args=ap.parse_args()
    df, _, physics = prepare(pd.read_csv(args.composition))
    site = pd.read_csv(args.sites)
    if site.record_id.duplicated().any():
        raise ValueError('Duplicate record IDs in site data')
    m=df[['record_id','material','adsorbate','ref_ads_eng_eV']].merge(site,on='record_id',how='inner',validate='one_to_one',suffixes=('','_site'))
    if len(m)!=len(df):
        raise ValueError(f'Missing site records: {len(m)} / {len(df)} matched. Resolve missing rows before comparison.')
    for col in ('material','adsorbate'):
        if not (m[col].astype(str)==m[col+'_site'].astype(str)).all():
            raise ValueError('Mismatch in '+col+' between datasets')
    if not np.allclose(m.ref_ads_eng_eV,m.energy_eV,atol=1e-7,rtol=0):
        raise ValueError('Energy mismatch between source datasets')
    order=pd.Index(df.record_id).get_indexer(m.record_id)
    if (order<0).any(): raise ValueError('Invalid row alignment')
    base=physics.iloc[order].reset_index(drop=True)
    if m[SITE_NUMERIC].isna().any().any():
        raise ValueError('Missing numeric site features')
    loc=pd.concat([base,m[SITE_NUMERIC].reset_index(drop=True).astype(float),
                  pd.get_dummies(m.site_type_proxy,prefix='site_type',dtype=float).reset_index(drop=True)],axis=1)
    y=m.ref_ads_eng_eV.to_numpy()
    fe_ni=m.material.map(lambda s: {'Fe','Ni'}.issubset(parse_material(s))).to_numpy()
    valid=np.where(~fe_ni)[0]
    if m.material.iloc[valid].nunique()<5: raise ValueError('Not enough groups')
    folds=list(GroupKFold(n_splits=5).split(base.iloc[valid],y[valid],m.material.iloc[valid]))
    rows=[]; preds=[]
    for representation,X in [('composition_plus_physics',base),('plus_local_relaxed_geometry',loc)]:
      for model_name in ('RF','SVR'):
        for fold,(itr,iva) in enumerate(folds,1):
            tr,va=valid[itr],valid[iva]
            reg=model_factory(model_name)
            reg.fit(X.iloc[tr],y[tr])
            yp=reg.predict(X.iloc[va]); errors=np.abs(y[va]-yp)
            rows.append(dict(representation=representation,model=model_name,fold=fold,
                             MAE_eV=float(np.mean(errors)),n=len(va)))
            for ix,pred,err in zip(va,yp,errors):
                preds.append(dict(representation=representation,model=model_name,fold=fold,
                                  record_id=m.record_id.iloc[ix],material=m.material.iloc[ix],
                                  adsorbate=m.adsorbate.iloc[ix],actual_eV=y[ix],
                                  predicted_eV=pred,abs_error_eV=err))
    o=Path(args.outdir);o.mkdir(parents=True,exist_ok=True)
    folds_df=pd.DataFrame(rows)
    predictions=pd.DataFrame(preds)
    summ=folds_df.groupby(['representation','model'],as_index=False).agg(
        CV_MAE_mean_eV=('MAE_eV','mean'),CV_MAE_std_eV=('MAE_eV','std'))
    by_ads=predictions.groupby(['representation','model','adsorbate'],as_index=False).agg(
        MAE_eV=('abs_error_eV','mean'),n=('abs_error_eV','size'))
    folds_df.to_csv(o/'fold_results.csv',index=False)
    summ.to_csv(o/'summary.csv',index=False)
    by_ads.to_csv(o/'by_adsorbate.csv',index=False)
    predictions.to_csv(o/'cv_predictions.csv',index=False)
    print('Matched records:',len(m),'FeNi reserved:',int(fe_ni.sum()))
    print('\n=== GROUPED CV COMPARISON ===')
    print(summ.to_string(index=False,float_format=lambda x:f'{x:.4f}'))
    print('\n=== BY ADSORBATE ===')
    print(by_ads.to_string(index=False,float_format=lambda x:f'{x:.4f}'))
    print('\nCAUTION: Relaxed site geometry may be unavailable before DFT, so this is a retrospective benchmark, not proof of pre-DFT screening capability.')
    print('Saved:',o)

if __name__=='__main__': main()
