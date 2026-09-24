#!/usr/bin/env python3
"""Per-chain mass-weighted Rg, end-to-end vector and backbone bond P2 along z.
Mapping: {"chains":[{"molecule_id":1,"backbone_ids":[1,4,7,...]}]}.
Requires complete polymer molecules and correct image flags. Entanglement is NOT inferred.
"""
import argparse,csv,json
import numpy as np
from postprocess import dumps


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('dump');p.add_argument('mapping');p.add_argument('--output',default='chain_structure.csv')
    p.add_argument('--stride',type=int,default=1);p.add_argument('--ppa-export',help='Optional JSONL backbone coordinate/topology exchange file, NOT native Z1 input')
    a=p.parse_args();chains=json.load(open(a.mapping))['chains'];ppa=open(a.ppa_export,'w') if a.ppa_export else None
    with open(a.output,'w',newline='') as f:
        w=csv.writer(f);w.writerow(['step','molecule_id','Rg2_A2','Re2_A2','Re_x_A','Re_y_A','Re_z_A','P2_bond_z'])
        for i,(step,cell,d) in enumerate(dumps(a.dump)):
            if i%a.stride:continue
            idx={int(k):j for j,k in enumerate(d['id'])};r=np.c_[d['xu'],d['yu'],d['zu']];export=[]
            for ch in chains:
                mid=ch['molecule_id'];mask=d['mol']==mid;b=np.array([idx[k] for k in ch['backbone_ids']])
                if not mask.any() or len(b)<2 or not np.all(mask[b]):raise ValueError('Invalid molecule/backbone mapping')
                rc=r[mask];mass=d['mass'][mask];com=np.average(rc,axis=0,weights=mass)
                rg2=np.average(np.sum((rc-com)**2,axis=1),weights=mass);bb=r[b];re=bb[-1]-bb[0]
                bonds=np.diff(bb,axis=0);lens=np.linalg.norm(bonds,axis=1)
                if (lens<=0).any():raise ValueError('Zero backbone bond length')
                p2=np.mean((3*(bonds[:,2]/lens)**2-1)/2)
                w.writerow([step,mid,rg2,re@re,*re,p2])
                export.append({'molecule_id':mid,'atom_ids':ch['backbone_ids'],'xyz_unwrapped_A':bb.tolist()})
            if ppa:ppa.write(json.dumps({'step':step,'cell_rows_A':cell.tolist(),'chains':export})+'\n')
    if ppa:ppa.close()

if __name__=='__main__':main()
