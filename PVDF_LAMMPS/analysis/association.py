#!/usr/bin/env python3
"""Li--partner contacts and intermittent/continuous contact survival from whole molecules.
Use partner donor atom TYPES and a cutoff chosen from the corresponding RDF minimum.
No-direct-anion-contact is NOT automatically a solvent-separated ion pair.
Orthogonal bulk boxes only; sample sufficiently fast to resolve contact breaking.
"""
import argparse,csv
import numpy as np
from postprocess import dumps,ids


def contact_survival(h,lags):
    out=[]
    for k in lags:
        initial=h[:-k]
        denom=initial.sum()
        if denom==0:out.append([k,np.nan,np.nan]);continue
        intermittent=(initial & h[k:]).sum()/denom
        continuous=initial.copy()
        for j in range(1,k+1):continuous &= h[j:len(h)-k+j]
        out.append([k,intermittent,continuous.sum()/denom])
    return np.asarray(out)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('dump');p.add_argument('--li-ids',required=True)
    p.add_argument('--donor-types',type=int,nargs='+',required=True);p.add_argument('--cutoff',type=float,required=True)
    p.add_argument('--dt-fs',type=float,required=True);p.add_argument('--maxlag-frames',type=int,default=100)
    p.add_argument('--stride',type=int,default=1);p.add_argument('--output',default='association.csv');p.add_argument('--survival',default='contact_survival.dat')
    a=p.parse_args();liids=ids(a.li_ids);history=[];steps=[];reference=None
    with open(a.output,'w',newline='') as f:
        w=csv.writer(f);w.writerow(['step','donor_CN','partner_molecule_CN','fraction_no_contact','fraction_one_partner','fraction_multi_partner'])
        for i,(step,cell,d) in enumerate(dumps(a.dump)):
            if i%a.stride:continue
            if not np.allclose(cell,np.diag(np.diag(cell))):raise ValueError('Orthogonal cells only')
            if reference is None:reference=d['id'].copy()
            elif not np.array_equal(reference,d['id']):raise ValueError('Atom IDs change')
            li=np.array([int(k) in liids for k in d['id']]);don=np.isin(d['type'].astype(int),a.donor_types)
            if not li.any() or not don.any():raise ValueError('Empty Li/donor selection')
            mol=d['mol'][don].astype(int);unique=np.unique(mol)
            if np.any(unique<=0):raise ValueError('Partners need positive molecule IDs')
            r=np.c_[d['xu'],d['yu'],d['zu']];dr=r[li,None,:]-r[None,don,:]
            dr-=np.rint(dr/np.diag(cell))*np.diag(cell)
            contact=np.linalg.norm(dr,axis=-1)<a.cutoff
            bymol=np.stack([np.any(contact[:,mol==m],axis=1) for m in unique],axis=1)
            nc=bymol.sum(axis=1)
            w.writerow([step,contact.sum(axis=1).mean(),nc.mean(),(nc==0).mean(),(nc==1).mean(),(nc>1).mean()])
            history.append(bymol);steps.append(step)
    if len(steps)<4 or not np.all(np.diff(steps)==np.diff(steps)[0]):raise ValueError('Need >=4 uniform frames')
    lagmax=min(a.maxlag_frames,len(steps)//3)
    lags=np.unique(np.linspace(1,lagmax,min(100,lagmax)).astype(int))
    s=contact_survival(np.array(history),lags);s[:,0]*=np.diff(steps)[0]*a.dt_fs
    np.savetxt(a.survival,s,header='lag_fs C_intermit C_continuous; donor definition/cutoff/frame interval dependent; NaN=no initial contacts')

if __name__=='__main__':main()
