#!/usr/bin/env python3
"""Group-safe comparison of composition-only vs physically informed descriptors.
No site descriptor is fabricated. Uses a user-local outputs/feni_adsorption.csv.
"""
import argparse, json, re, sys
from pathlib import Path
import numpy as np
import pandas as pd
from ase.data import atomic_numbers, atomic_masses, covalent_radii
from sklearn.base import clone
from sklearn.dummy import DummyRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import GroupKFold, GroupShuffleSplit, GridSearchCV
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR
from sklearn.neural_network import MLPRegressor

PATTERN = re.compile(r'([A-Z][a-z]?)(\d+)')

def parse_material(s):
    pairs=PATTERN.findall(str(s))
    if not pairs or ''.join(f'{a}{b}' for a,b in pairs)!=s:
        raise ValueError(f'Unsupported material formula {s!r}')
    counts={}
    for el,n in pairs:
        if el not in atomic_numbers: raise ValueError(f'Unknown element: {el}')
        counts[el]=counts.get(el,0)+int(n)
    return counts

def descriptors(c):
    total=sum(c.values())
    elems=list(c)
    weights=np.array([c[e]/total for e in elems])
    result={}
    for label, values in {
        'Z':[atomic_numbers[e] for e in elems],
        'atomic_mass':[atomic_masses[atomic_numbers[e]] for e in elems],
        'covalent_radius':[covalent_radii[atomic_numbers[e]] for e in elems],
    }.items():
        v=np.array(values,dtype=float)
        if not np.isfinite(v).all(): raise ValueError(f'Missing descriptor: {elems}, {label}')
        avg=float(weights@v)
        result[f'{label}_weighted_mean']=avg
        result[f'{label}_weighted_std']=float(np.sqrt(weights@((v-avg)**2)))
        result[f'{label}_min']=float(v.min())
        result[f'{label}_max']=float(v.max())
    result['element_count']=len(elems)
    return result

def prepare(df):
    if not {'record_id','material','adsorbate','ref_ads_eng_eV'}.issubset(df.columns):
        raise ValueError('Missing required columns')
    df=df.copy().dropna(subset=['ref_ads_eng_eV','material','adsorbate'])
    df['ref_ads_eng_eV']=pd.to_numeric(df['ref_ads_eng_eV'],errors='raise')
    df=df[np.isfinite(df['ref_ads_eng_eV'])].drop_duplicates('record_id').reset_index(drop=True)
    counts=[parse_material(m) for m in df.material]
    elems=sorted({e for row in counts for e in row})
    composition=pd.DataFrame([{e:c.get(e,0)/sum(c.values()) for e in elems} for c in counts], dtype=float)
    physical=pd.DataFrame([descriptors(c) for c in counts])
    adsorption=pd.get_dummies(df.adsorbate,prefix='ads',dtype=float)
    return df, pd.concat([composition,adsorption],axis=1),pd.concat([composition,physical,adsorption],axis=1)

def fresh_model(kind):
    if kind=='SVR':
        return Pipeline([('scale',StandardScaler()),('model',SVR(kernel='rbf',C=10,epsilon=0.1))])
    if kind=='RF':
        return RandomForestRegressor(n_estimators=180,min_samples_leaf=3,random_state=42,n_jobs=1)
    if kind=='MLP':
        return Pipeline([('scale',StandardScaler()),('model',MLPRegressor(hidden_layer_sizes=(32,16),alpha=1.0,max_iter=1200,random_state=42,early_stopping=False))])
    raise KeyError(kind)

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--input',default='outputs/feni_adsorption.csv')
    p.add_argument('--output',default='outputs/descriptor_study')
    p.add_argument('--holdout-feni',action='store_true',help='Exploratory ONLY: FeNi has been repeatedly used for selection; not an independent final test')
    p.add_argument('--models',nargs='+',choices=['SVR','RF','MLP'],default=['SVR','RF'])
    a=p.parse_args()
    df,base,phys=prepare(pd.read_csv(a.input))
    out=Path(a.output); out.mkdir(parents=True,exist_ok=True)
    (out/'descriptor_columns.json').write_text(json.dumps({'composition':list(base),'physical':list(phys)},indent=2))
    mask=np.array([{'Fe','Ni'}.issubset(parse_material(m)) for m in df.material])
    fitidx=np.where(~mask)[0] if a.holdout_feni else np.arange(len(df))
    if len(fitidx)<20: raise SystemExit('Too few training records')
    groups=df.material.iloc[fitidx].to_numpy()
    if len(set(groups))<5: raise SystemExit('Need >=5 independent material groups')
    y=df.ref_ads_eng_eV.to_numpy()
    folds=list(GroupKFold(n_splits=5).split(base.iloc[fitidx],y[fitidx],groups))
    results=[]; per_fold=[]
    for space,X in [('composition',base),('composition_plus_physics',phys)]:
        for name in ['AdsorbateMean']+a.models:
            scores=[]
            for fold,(tr,va) in enumerate(folds,1):
                train_i=fitidx[tr]; val_i=fitidx[va]
                if name=='AdsorbateMean':
                    means=df.iloc[train_i].groupby('adsorbate')['ref_ads_eng_eV'].mean()
                    pred=df.iloc[val_i].adsorbate.map(means).fillna(float(y[train_i].mean())).to_numpy()
                else:
                    model=fresh_model(name)
                    model.fit(X.iloc[train_i],y[train_i]); pred=model.predict(X.iloc[val_i])
                mae=mean_absolute_error(y[val_i],pred); scores.append(mae)
                per_fold.append(dict(feature_set=space,model=name,fold=fold,MAE_eV=float(mae),n_validation=len(val_i),n_materials=df.material.iloc[val_i].nunique()))
            results.append(dict(feature_set=space,model=name,MAE_mean_eV=float(np.mean(scores)),MAE_std_eV=float(np.std(scores)),N_train=len(fitidx),N_groups=len(set(groups))))
    summary=pd.DataFrame(results).sort_values('MAE_mean_eV')
    summary.to_csv(out/'grouped_cv_summary.csv',index=False)
    pd.DataFrame(per_fold).to_csv(out/'grouped_cv_folds.csv',index=False)
    print('=== GROUPED CV (NO FeNi HOLDOUT USED FOR MODEL SELECTION) ===')
    print(summary.to_string(index=False,float_format=lambda f:f'{f:.4f}'))
    print('\nNOTE: features are GLOBAL composition descriptors, NOT local adsorption-site descriptors.')
    if a.holdout_feni:
        print(f'FeNi reserved from CV: {mask.sum()} records, {df.material[mask].nunique()} materials. NOT evaluated in this script.')
    print('Saved:',out)
if __name__=='__main__': main()
