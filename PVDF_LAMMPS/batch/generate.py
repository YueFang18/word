#!/usr/bin/env python3
"""Generate independent per-model/case/seed commands; never runs LAMMPS."""
import argparse,csv,hashlib,json,re,shlex
from collections import Counter
from pathlib import Path


def descriptors(seq):
    if not seq: raise ValueError('Empty sequence')
    n=len(seq);counts=Counter(seq);runs=[]
    for x in seq:
        if runs and runs[-1][0]==x:runs[-1][1]+=1
        else:runs.append([x,1])
    out={'DP':n,'monomer_fractions':{k:v/n for k,v in sorted(counts.items())},
         'heterodyad_fraction':sum(a!=b for a,b in zip(seq,seq[1:]))/(n-1) if n>1 else 0,
         'mean_run_length':n/len(runs),'max_run_length':max(v for _,v in runs),
         'mean_run_by_monomer':{k:sum(v for s,v in runs if s==k)/sum(s==k for s,v in runs) for k in counts}}
    for size,name in [(2,'dyads'),(3,'triads')]:
        c=Counter(tuple(seq[i:i+size]) for i in range(n-size+1))
        out[name]={','.join(k):v/max(n-size+1,1) for k,v in sorted(c.items())}
    return out


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('manifest');p.add_argument('--output',default='batch/generated')
    p.add_argument('--lmp',default='lmp');p.add_argument('--require-models',action='store_true');a=p.parse_args()
    root=Path(__file__).resolve().parents[1];dest=(root/a.output).resolve();dest.mkdir(parents=True,exist_ok=True)
    manifest=json.loads(Path(a.manifest).read_text());records=[];commands=[];used=set()
    for m in manifest['models']:
        if not re.fullmatch(r'[A-Za-z0-9_-]+',m['id']):raise ValueError('Use safe model IDs')
        desc=descriptors(m['sequence_tokens'])
        for c in m['cases']:
            if not re.fullmatch(r'[A-Za-z0-9_-]+',c['name']):raise ValueError('Use safe case names')
            for seed in manifest['seeds']:
                key=f"{m['id']}__{c['name']}__s{seed}"
                if key in used:raise ValueError('Duplicate model/case/seed')
                used.add(key)
                out=f'results/{key}'
                args=[a.lmp,'-in','lammps.in','-var','task',c['task'],'-var','seed',str(seed),'-var','out',out]
                if 'mode' in c:args+=['-var','mode',c['mode']]
                checks={}
                for k in ['data','styles','coeffs','groups','hooks']:
                    v=m[k]
                    if any(ch.isspace() for ch in v) or any(ch in v for ch in '$\"\''):
                        raise ValueError('Model paths must have no whitespace, quotes or $ for LAMMPS substitution')
                    file=(root/v).resolve();checks[k]={'path':str(file),'exists':file.is_file()}
                    if file.is_file():checks[k]['sha256']=hashlib.sha256(file.read_bytes()).hexdigest()
                    if a.require_models and not file.is_file():raise FileNotFoundError(file)
                    args+=['-var',k,v]
                reserved={'data','styles','coeffs','groups','hooks','task','mode','out','seed'}
                if reserved.intersection(c.get('vars',{})):raise ValueError('Reserved variable in case vars')
                for k,v in c.get('vars',{}).items():args+=['-var',k,str(v)]
                # No shell interpolation of user text; shell quoting preserves argv boundaries.
                command=shlex.join(args)
                commands.append(f'if [ -e {shlex.quote(out)} ]; then echo {shlex.quote("Refusing overwrite: "+out)} >&2; exit 1; fi\n{command}')
                records.append({'run_id':key,'model':m,'case':c,'seed':seed,'descriptors':desc,'files':checks,'command':command})
    script='#!/usr/bin/env bash\nset -euo pipefail\ncd '+shlex.quote(str(root))+'\n\n'+'\n\n'.join(commands)+'\n'
    (dest/'run_all.sh').write_text(script);(dest/'run_all.sh').chmod(0o755)
    (dest/'run_manifest.json').write_text(json.dumps(records,ensure_ascii=False,indent=2)+'\n')
    with (dest/'sequence_features.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['model_id','DP','heterodyad_fraction','mean_run_length','max_run_length','monomer_fractions_json','dyads_json','triads_json'])
        for m in manifest['models']:
            d=descriptors(m['sequence_tokens']);w.writerow([m['id'],d['DP'],d['heterodyad_fraction'],d['mean_run_length'],d['max_run_length'],json.dumps(d['monomer_fractions']),json.dumps(d['dyads']),json.dumps(d['triads'])])
    print(f'Generated {len(records)} independent jobs in {dest}; not executed.')

if __name__=='__main__':main()
