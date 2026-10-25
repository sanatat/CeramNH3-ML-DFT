#!/usr/bin/env python3
"""Conservative demonstrator: no trained model unless >=30 observations from >=5 distinct papers."""
from pathlib import Path
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import OneHotEncoder
from sklearn.pipeline import make_pipeline
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GroupKFold, cross_val_predict
from sklearn.metrics import mean_absolute_error
ROOT=Path(__file__).resolve().parents[1]
df=pd.read_csv(ROOT/'data/adsorption_template.csv')
required=['record_id','paper_id','material','facet','site','adsorbate','xc_functional','energy_eV','energy_definition','source_location']
if any(c not in df.columns for c in required): raise SystemExit('Missing required columns')
df=df.dropna(subset=required).copy()
if df.empty:
    print('NOT TRAINED: no vetted adsorption-energy examples. Populate data/adsorption_template.csv from primary sources.')
    raise SystemExit(0)
# Avoid silently mixing energies with different definitions / DFT functionals.
sets=df.groupby(['energy_definition','xc_functional'],dropna=False)
valid=[(key,g) for key,g in sets if len(g)>=30 and g.paper_id.nunique()>=5]
if not valid:
    print('NOT TRAINED: need >=30 comparable records from >=5 independent papers for one XC + energy-definition group.')
    for key,g in sets: print(key, 'rows',len(g),'papers',g.paper_id.nunique())
    raise SystemExit(0)
key,df=max(valid,key=lambda x:len(x[1]))
features=['material','facet','site','adsorbate']
X=df[features]; y=df.energy_eV.astype(float)
pre=ColumnTransformer([('cat',OneHotEncoder(handle_unknown='ignore'),features)])
model=make_pipeline(pre,RandomForestRegressor(n_estimators=200,random_state=42,min_samples_leaf=2,n_jobs=1))
n_splits=min(5,df.paper_id.nunique())
pred=cross_val_predict(model,X,y,groups=df.paper_id,cv=GroupKFold(n_splits=n_splits),n_jobs=1)
print('Held-out-paper MAE (eV):',round(mean_absolute_error(y,pred),4), 'for group',key)
print('MAE is not proof of physical transferability to proton-conducting ceramic cells.')
