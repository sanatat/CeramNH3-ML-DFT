#!/usr/bin/env python3
"""Only plot an entire, explicitly curated, internally consistent CHE profile.
Cumulative step ΔG° values MUST come from one explicitly specified study/system.
No adsorption energies are automatically converted to step energies.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

STEPS=['N2(g)+*','N2*','NNH*','NNH2*','NHNH2*','NH2NH2*','NH2*+NH3','2NH3+*']

def main():
    a=argparse.ArgumentParser()
    a.add_argument('--input',default='data/nrr_step_evidence_template.csv')
    a.add_argument('--outdir',default='outputs/nrr_v10')
    args=a.parse_args(); out=Path(args.outdir);out.mkdir(parents=True,exist_ok=True)
    df=pd.read_csv(args.input,dtype=str).fillna('')
    required={'material','source_doi','source_figure_or_table','pathway_id','step_index','state','delta_G_eV','reference_conditions','evidence_type'}
    if missing:=required-set(df.columns):raise SystemExit(f'Missing columns: {sorted(missing)}')
    drawn=0; status=[]
    for (material,source,pathway),g in df.groupby(['material','source_doi','pathway_id'],dropna=False):
        if not source or source.upper()=='TBD':status.append((material,pathway,'missing source'));continue
        g=g.copy()
        # No silent mixing of conditions or simulated vs experimental provenance.
        if g.reference_conditions.nunique()!=1 or not g.reference_conditions.iloc[0]:
            status.append((material,pathway,'mixed/missing conditions'));continue
        if any(x=='' for x in g.delta_G_eV) or len(g)<7:
            status.append((material,pathway,'incomplete step energies'));continue
        try:
            g['idx']=pd.to_numeric(g.step_index,errors='raise').astype(int)
            g['dg']=pd.to_numeric(g.delta_G_eV,errors='raise').astype(float)
        except (ValueError,TypeError):status.append((material,pathway,'nonnumeric index/energy'));continue
        g=g.sort_values('idx')
        if list(g.idx)!=list(range(1,len(g)+1)) or not np.isfinite(g.dg).all():
            status.append((material,pathway,'invalid steps'));continue
        if any(x not in ('DFT_CHE','EXPERIMENT_DERIVED') for x in g.evidence_type):
            status.append((material,pathway,'unsupported evidence type'));continue
        if g.evidence_type.nunique()!=1:status.append((material,pathway,'mixed evidence types'));continue
        cumulative=np.r_[0.,np.cumsum(g.dg.values)]
        states=['start']+g.state.tolist()
        fig,ax=plt.subplots(figsize=(10,4.7))
        ax.plot(range(len(cumulative)),cumulative,marker='o')
        ax.set_xticks(range(len(cumulative)),states,rotation=40,ha='right')
        ax.set_ylabel('Cumulative ΔG (eV)')
        ax.set_title(f'{material} • {pathway}\n{source} | {g.reference_conditions.iloc[0]}')
        fig.tight_layout()
        filename=out/(f'{material}_{pathway}'.replace('/','_')+'.png')
        fig.savefig(filename,dpi=170);plt.close(fig)
        drawn+=1;status.append((material,pathway,'plotted'))
    pd.DataFrame(status,columns=['material','pathway_id','status']).to_csv(out/'profile_plot_status.csv',index=False)
    print('Verified complete profiles plotted:',drawn)
    if not drawn:print('No comparable complete CHE pathway entered; no diagram fabricated.')
    print('Status:',out/'profile_plot_status.csv')

if __name__=='__main__':main()
