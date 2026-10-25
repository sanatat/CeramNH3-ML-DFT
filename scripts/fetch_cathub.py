#!/usr/bin/env python3
"""Fetch limited, attributed Catalysis-Hub reaction records (no fabricated energies).
Requires network on user's computer. GraphQL schema and availability may change.
Usage: python scripts/fetch_cathub.py --first 50 --surface FeNi
"""
import argparse,csv,json,sys,urllib.request,urllib.error
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
API='https://api.catalysis-hub.org/graphql'
p=argparse.ArgumentParser();p.add_argument('--first',type=int,default=30);p.add_argument('--surface',default='FeNi');p.add_argument('--url',default=API);a=p.parse_args()
if not (1<=a.first<=200):p.error('--first must be 1..200')
if not all(ch.isalnum() or ch in '-_' for ch in a.surface):p.error('invalid surface filter')
# Documentary filter syntax: chemicalComposition; some servers may require a different match convention.
query="""{ reactions(first:%d, chemicalComposition:"~%s") { totalCount edges { node { id Equation chemicalComposition reactionEnergy activationEnergy } } } }"""%(a.first,a.surface)
request=urllib.request.Request(a.url,data=json.dumps({'query':query}).encode(),headers={'Content-Type':'application/json','User-Agent':'CeramNH3-LiteratureML/0.4'},method='POST')
try:
 with urllib.request.urlopen(request,timeout=35) as r: obj=json.load(r)
except (urllib.error.URLError,TimeoutError,ValueError) as exc:
 print('Fetch failed (database/API/network):',exc,file=sys.stderr);sys.exit(2)
if obj.get('errors'):
 print('GraphQL errors:',json.dumps(obj['errors'],indent=2),file=sys.stderr);sys.exit(3)
edges=((obj.get('data') or {}).get('reactions') or {}).get('edges') or []
out=ROOT/'data'/'cathub_raw';out.mkdir(exist_ok=True)
path=out/f'reactions_{a.surface}.csv'
fields=['id','Equation','chemicalComposition','reactionEnergy','activationEnergy','database','source_query','verification_status']
with path.open('w',newline='') as fh:
 w=csv.DictWriter(fh,fieldnames=fields);w.writeheader()
 for e in edges:
  n=(e or {}).get('node') or {}
  w.writerow({**{k:n.get(k) for k in fields[:5]},'database':'Catalysis-Hub','source_query':query,'verification_status':'unvetted'})
print('Fetched',len(edges),'unvetted reaction records to',path)
print('IMPORTANT: reactionEnergy is NOT necessarily adsorption energy; do not copy to adsorption_template.csv without verifying equation, reference, site, XC, publication.')
