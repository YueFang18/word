#!/usr/bin/env python3
"""Export a selected MD frame's COMPLETE molecules and paired ORCA vertical-IP inputs.
No molecule truncation/capping, no automatic solvent shell, no inferred formal charge.
Type mapping example {"1":"C","2":"H","3":"F"}. Gas-phase diagnostic by default.
"""
import argparse,json
from pathlib import Path
import numpy as np
from postprocess import dumps


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('dump');p.add_argument('--type-map',required=True)
    p.add_argument('--frame',type=int,default=0);p.add_argument('--molecules',type=int,nargs='+',required=True)
    p.add_argument('--center-id',type=int,required=True);p.add_argument('--charge',type=int,required=True)
    p.add_argument('--mult-initial',type=int,required=True);p.add_argument('--mult-oxidized',type=int,required=True)
    p.add_argument('--method',default='wB97X-D3 def2-TZVP TightSCF');p.add_argument('--nprocs',type=int,default=4)
    p.add_argument('--output',required=True);a=p.parse_args();out=Path(a.output)
    if out.exists():raise FileExistsError('Use a new output directory to preserve QM inputs')
    found=False
    for i,(step,cell,d) in enumerate(dumps(a.dump)):
        if i!=a.frame:continue
        found=True;break
    if not found:raise ValueError('Frame not found')
    typemap=json.loads(Path(a.type_map).read_text());atomid=d['id'].astype(int);mol=d['mol'].astype(int)
    if any(m<=0 for m in a.molecules):raise ValueError('Positive molecular IDs required; surfaces need a dedicated periodic QM model')
    if any(m not in mol for m in a.molecules):raise ValueError('Requested molecule missing')
    center=np.where(atomid==a.center_id)[0]
    if len(center)!=1 or mol[center[0]] not in a.molecules:raise ValueError('Center atom must be in selected molecules')
    xyz=np.c_[d['xu'],d['yu'],d['zu']];origin=xyz[center[0]].copy();mask=np.isin(mol,a.molecules)
    # Translate complete molecular images to a common local neighborhood.
    # Check very extended polymers manually; no cutting of a macromolecule occurs here.
    for m in a.molecules:
        g=mol==m;com=np.average(xyz[g],axis=0,weights=d['mass'][g]);frac=(com-origin)@np.linalg.inv(cell)
        xyz[g]-=np.rint(frac)@cell
    xyz-=origin
    lines=[];selected=[]
    for j in np.where(mask)[0]:
        el=typemap[str(int(d['type'][j]))]
        lines.append(f'{el:3s} {xyz[j,0]: .10f} {xyz[j,1]: .10f} {xyz[j,2]: .10f}')
        selected.append({'id':int(atomid[j]),'molecule_id':int(mol[j]),'element':el})
    out.mkdir(parents=True)
    (out/'snapshot.xyz').write_text(f'{len(lines)}\nMD step {step}; full molecules; inspect cluster boundaries\n'+'\n'.join(lines)+'\n')
    for name,chg,mult in [('initial_vertical',a.charge,a.mult_initial),('oxidized_vertical',a.charge+1,a.mult_oxidized)]:
        (out/(name+'.inp')).write_text(f'! {a.method}\n%pal nprocs {a.nprocs} end\n%maxcore 2000\n* xyzfile {chg} {mult} snapshot.xyz\n')
    (out/'selection.json').write_text(json.dumps({'step':step,'atoms':selected,'initial_charge':a.charge,'oxidized_charge':a.charge+1,'initial_multiplicity':a.mult_initial,'oxidized_multiplicity':a.mult_oxidized,'method':a.method,'warning':'Default is a gas-phase vertical IP diagnostic. Both calculations MUST use identical geometry/environment. Explicit environment, embedding/solvation, spin-state validation and reference calibration are needed for condensed-phase oxidation potentials. No electron removal occurs in LAMMPS.'},indent=2)+'\n')
    print(out)

if __name__=='__main__':main()
