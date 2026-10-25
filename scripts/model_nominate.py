#!/usr/bin/env python3
"""Out-of-group ML candidate nomination using composition-only inputs.
This is a screening proxy, NOT a claim about ammonia yield/activity.
"""
import argparse
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, r2_score
sys.path.insert(0,str(Path(__file__).resolve().parent))
from run_descriptor_study import prepare, parse_material

TARGETS=('H*','N*','NH*')

def new_model():
    return RandomForestRegressor(n_estimators=180,min_samples_leaf=3,random_state=42,n_jobs=1)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--input',default='outputs/feni_adsorption.csv')
    ap.add_argument('--outdir',default='outputs/model_nominations')
    ap.add_argument('--top',type=int,default=15)
    args=ap.parse_args()
    df,_,X=prepare(pd.read_csv(args.input))
    mask=df.material.map(lambda m:{'Fe','Ni'}.issubset(parse_material(m))).to_numpy()
    non=np.flatnonzero(~mask); feni=np.flatnonzero(mask)
    if df.material.iloc[non].nunique()<5:raise ValueError('At least 5 training compositions needed')
    y=df.ref_ads_eng_eV.to_numpy()
    predicted=np.full(len(df),np.nan)
    fold_id=np.zeros(len(df),dtype=int)
    folds=GroupKFold(5)
    for f,(tr,va) in enumerate(folds.split(X.iloc[non],y[non],df.material.iloc[non]),1):
        tr_idx,va_idx=non[tr],non[va]
        model=new_model();model.fit(X.iloc[tr_idx],y[tr_idx]);predicted[va_idx]=model.predict(X.iloc[va_idx]);fold_id[va_idx]=f
    model=new_model();model.fit(X.iloc[non],y[non])
    if len(feni): predicted[feni]=model.predict(X.iloc[feni]); fold_id[feni]=-1
    assert np.isfinite(predicted).all()
    out=Path(args.outdir);out.mkdir(parents=True,exist_ok=True)
    p=df[['record_id','material','adsorbate','ref_ads_eng_eV']].copy()
    p['prediction_eV']=predicted
    p['fold_id']=fold_id
    p['prediction_provenance']=np.where(mask,'FeNi held out from training (exploratory)', 'OOF grouped 5-fold')
    p['absolute_error_eV']=np.abs(y-predicted)
    p.to_csv(out/'out_of_group_predictions.csv',index=False)
    print('COMPOSITION-ONLY predictions from models that did not train on the candidate material')
    print('Records:',len(p),'unique materials:',p.material.nunique())
    print('OOF excluding FeNi MAE:',round(mean_absolute_error(y[non],predicted[non]),4), 'R2:',round(r2_score(y[non],predicted[non]),4))
    if len(feni):print('FeNi exploratory MAE:',round(mean_absolute_error(y[feni],predicted[feni]),4))

    # Composition-only RF predicts the same mean adsorption energy for all sites
    # of a given material/adsorbate. Aggregate rather than multiplying site counts.
    table=(p.groupby(['material','adsorbate'],as_index=False)
            .agg(prediction_eV=('prediction_eV','mean'),n_sites=('record_id','size')))
    matrix=table.pivot(index='material',columns='adsorbate',values='prediction_eV')
    counts=table.pivot(index='material',columns='adsorbate',values='n_sites')
    complete=matrix.reindex(columns=list(TARGETS)).dropna().copy()
    if complete.empty:
        print('No materials with all H*, N*, NH* records. Report the coverage instead.')
        print(matrix.notna().sum()); return
    # Transparent two-objective Pareto screening heuristic: lower N* binding
    # reference energy and higher H* reference energy. Not kinetic activity.
    # Extreme N* overbinding may poison a surface; this is NOT an ammonia optimum.
    vals=complete[['N*','H*']].to_numpy()
    dominated=np.zeros(len(complete),bool)
    for i,(n,h) in enumerate(vals):
        dominated[i]=bool(np.any((vals[:,0]<=n)&(vals[:,1]>=h)&((vals[:,0]<n)|(vals[:,1]>h))))
    shortlisted=complete.copy()
    shortlisted['pareto_n_vs_h']=~dominated
    shortlisted['FeNi_group']=[{'Fe','Ni'}.issubset(parse_material(v)) for v in shortlisted.index]
    shortlisted['n_H_records']=counts.loc[shortlisted.index,'H*'].astype(int)
    shortlisted['n_N_records']=counts.loc[shortlisted.index,'N*'].astype(int)
    shortlisted['n_NH_records']=counts.loc[shortlisted.index,'NH*'].astype(int)
    shortlisted=shortlisted.reset_index().rename(columns={'H*':'pred_H_eV','N*':'pred_N_eV','NH*':'pred_NH_eV'})
    shortlisted['caution']='ADSORPTION PROXY ONLY: not NH3 rate/FE; N overbinding possible'
    shortlisted.sort_values(['pareto_n_vs_h','pred_N_eV','pred_H_eV'],ascending=[False,True,False]).to_csv(out/'all_candidate_materials.csv',index=False)
    pareto=shortlisted[shortlisted.pareto_n_vs_h].sort_values(['pred_N_eV','pred_H_eV'],ascending=[True,False])
    pareto.to_csv(out/'pareto_nominees.csv',index=False)
    print('\nCoverage: candidates with H*, N*, NH*:',len(complete),'Pareto nominees:',len(pareto))
    print('\n=== MODEL-NOMINATED PARETO CANDIDATES (N*/H* ADSORPTION PROXY; NOT AMMONIA ACTIVITY) ===')
    cols=['material','pred_H_eV','pred_N_eV','pred_NH_eV','FeNi_group']
    print(pareto[cols].head(args.top).to_string(index=False,float_format=lambda v:f'{v:.3f}'))
    print('\n=== FeNi candidates, if available ===')
    print(shortlisted[shortlisted.FeNi_group][cols+['pareto_n_vs_h']].to_string(index=False,float_format=lambda v:f'{v:.3f}'))
    print('\nSaved:',out)
    print('IMPORTANT: Do not select a single optimal NH3 catalyst from this heuristic alone.')
    print('Check chemical formulas, adsorption reference definitions, independent sites, reaction barriers and experimentally comparable NH3 FE/rate.')

if __name__=='__main__':main()
