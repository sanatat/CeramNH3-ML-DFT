#!/usr/bin/env python3
"""Retrospective catalyst-composition nominations from group-held-out, geometry-informed RF.
Not prospective NH3 activity predictions; needs relaxed DFT adsorption structures.
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from sklearn.metrics import mean_absolute_error, r2_score
from sklearn.ensemble import RandomForestRegressor

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from run_descriptor_study import prepare, parse_material
from compare_local_models import SITE_NUMERIC, model_factory

outdir = Path('outputs/local_nominations')
outdir.mkdir(parents=True, exist_ok=True)

# Reuse the precise out-of-fold predictions that produced R2=0.8089.
pred_path = Path('outputs/local_comparison/cv_predictions.csv')
p = pd.read_csv(pred_path)
p = p[(p['representation']=='plus_local_relaxed_geometry') & (p['model']=='RF')].copy()
assert not p.empty
assert not p.record_id.duplicated().any(), 'Duplicate OOF row'
p['provenance'] = 'GroupKFold OOF (non-FeNi)'
print('Non-FeNi OOF records:',len(p),'R2:',round(r2_score(p.actual_eV,p.predicted_eV),4), 'MAE:',round(mean_absolute_error(p.actual_eV,p.predicted_eV),4))

# Add FeNi using model fit on NON-FeNi exclusively (only for exploratory comparison).
df, _, physics = prepare(pd.read_csv('outputs/feni_adsorption.csv'))
sites = pd.read_csv('outputs/local_sites.csv')
assert not sites.record_id.duplicated().any()
m = df[['record_id','material','adsorbate','ref_ads_eng_eV']].merge(sites,on='record_id',how='inner',validate='one_to_one',suffixes=('','_site'))
if len(m)!=len(df): raise RuntimeError('Incomplete local_sites.csv')
for c in ['material','adsorbate']:
    if not (m[c].astype(str)==m[c+'_site'].astype(str)).all():
        raise RuntimeError(f'Mismatch on {c}')
if not np.allclose(m.ref_ads_eng_eV,m.energy_eV,atol=1e-7,rtol=0):
    raise RuntimeError('Energy mismatch')
order=pd.Index(df.record_id).get_indexer(m.record_id)
base=physics.iloc[order].reset_index(drop=True)
loc=pd.concat([base, m[SITE_NUMERIC].reset_index(drop=True).astype(float),
               pd.get_dummies(m.site_type_proxy,prefix='site_type',dtype=float).reset_index(drop=True)],axis=1)
mask=m.material.map(lambda s: {'Fe','Ni'}.issubset(parse_material(s))).to_numpy()
if mask.any():
    reg=model_factory('RF')
    reg.fit(loc.loc[~mask],m.ref_ads_eng_eV.to_numpy()[~mask])
    fe=m.loc[mask,['record_id','material','adsorbate','ref_ads_eng_eV']].copy()
    fe['actual_eV']=fe.ref_ads_eng_eV
    fe['predicted_eV']=reg.predict(loc.loc[mask])
    fe['provenance']='FeNi held-out (exploratory)'
    fe['fold']=-1
    fe=fe[['record_id','material','adsorbate','actual_eV','predicted_eV','provenance','fold']]
    print('FeNi local holdout n:',len(fe),'MAE:',round(mean_absolute_error(fe.actual_eV,fe.predicted_eV),4),'R2:',round(r2_score(fe.actual_eV,fe.predicted_eV),4))
    p=pd.concat([p[['record_id','material','adsorbate','actual_eV','predicted_eV','provenance','fold']],fe],ignore_index=True)
else:
    p=p[['record_id','material','adsorbate','actual_eV','predicted_eV','provenance','fold']]
p.to_csv(outdir/'site_predictions.csv',index=False)

# Site-level predictions can differ for same material; use transparent median + range.
stats=p.groupby(['material','adsorbate'],as_index=False).agg(
    predicted_median=('predicted_eV','median'),
    predicted_min=('predicted_eV','min'),
    predicted_max=('predicted_eV','max'),
    site_count=('record_id','size'))
wide=stats.pivot(index='material',columns='adsorbate',values='predicted_median')
counts=stats.pivot(index='material',columns='adsorbate',values='site_count')
complete=wide.reindex(columns=['H*','N*','NH*']).dropna().copy()
if complete.empty: raise RuntimeError('No compositions with H*, N*, NH* predictions')
# N more negative vs H more positive is an intentionally limited, unvalidated proxy.
vals=complete[['N*','H*']].to_numpy()
dominated=np.array([np.any((vals[:,0]<=n)&(vals[:,1]>=h)&((vals[:,0]<n)|(vals[:,1]>h))) for n,h in vals])
complete['pareto_n_vs_h']=~dominated
complete['FeNi_group']=[{'Fe','Ni'}.issubset(parse_material(s)) for s in complete.index]
for target in ['H*','N*','NH*']:
    complete[f'n_{target.replace("*", "")}_sites']=counts.loc[complete.index,target].astype(int)
complete=complete.rename(columns={'H*':'median_pred_H_eV','N*':'median_pred_N_eV','NH*':'median_pred_NH_eV'})
complete=complete.reset_index()
complete['screening_note']='Retrospective, relaxed-geometry N/H Pareto proxy; not NH3 activity'
complete=complete.sort_values(['pareto_n_vs_h','median_pred_N_eV'],ascending=[False,True])
complete.to_csv(outdir/'all_composition_candidates.csv',index=False)
short=complete[complete.pareto_n_vs_h].copy()
short.to_csv(outdir/'pareto_nominees.csv',index=False)
print('\nComplete (H,N,NH) materials:',len(complete),'Pareto nominees:',len(short))
print('\n=== Geometry-informed RF nominees (N/H proxy, NOT NH3 activity) ===')
print(short[['material','median_pred_H_eV','median_pred_N_eV','median_pred_NH_eV','FeNi_group']].to_string(index=False,float_format=lambda x:f'{x:.3f}'))
print('\n=== FeNi materials ===')
print(complete[complete.FeNi_group][['material','median_pred_H_eV','median_pred_N_eV','median_pred_NH_eV','pareto_n_vs_h']].to_string(index=False,float_format=lambda x:f'{x:.3f}'))
print('\nSaved:',outdir)
print('WARNING: No N2 activation barriers, ammonia FE/rate, or site-to-site reaction network. Do not claim an optimal ammonia catalyst.')
