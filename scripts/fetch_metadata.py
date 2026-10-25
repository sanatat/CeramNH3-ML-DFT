#!/usr/bin/env python3
"""Optional online DOI metadata retrieval via Crossref; does not download paper text or invent data."""
import json,urllib.request
from pathlib import Path
import pandas as pd
root=Path(__file__).resolve().parents[1]
dois=pd.read_csv(root/'data/papers.csv').doi.to_list()
result=[]
for doi in dois:
    url='https://api.crossref.org/works/'+doi
    try:
        req=urllib.request.Request(url,headers={'User-Agent':'CeramNH3LiteratureML/0.3 (academic research)'})
        with urllib.request.urlopen(req,timeout=15) as r: msg=json.load(r)['message']
        result.append({'doi':doi,'title':msg.get('title',[''])[0],'published':str(msg.get('published',{}).get('date-parts','')),'url':msg.get('URL','')})
    except Exception as exc:
        result.append({'doi':doi,'error':str(exc)})
out=root/'outputs/crossref_metadata.json';out.parent.mkdir(exist_ok=True)
out.write_text(json.dumps(result,indent=2,ensure_ascii=False))
print('Saved',out)
