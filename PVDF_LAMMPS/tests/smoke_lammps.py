#!/usr/bin/env python3
"""Run short NONPHYSICAL toy systems through fixed-charge script branches.
These are software smoke tests, not parameter validation or PVDF predictions.
ReaxFF/NEB are excluded because no validated reactive force field is supplied.
"""
import argparse,json,shutil,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]


def fixture(folder):
    folder.mkdir(parents=True,exist_ok=True);atoms=[];bonds=[];angles=[];torsions=[];grip=[];anchor=[]
    def atom(mol,typ,q,x,y,z):
        i=len(atoms)+1;atoms.append(f'{i} {mol} {typ} {q} {x} {y} {z} 0 0 0');return i
    for m in range(4):
        chain=[atom(m+1,1,0,10+8*(m%2)+(.3 if j%2 else 0),12+8*(m//2)+(.2 if j%3 else 0),16+2*j) for j in range(8)]
        grip.append(chain[-1])
        bonds.extend([[1,*chain[j:j+2]] for j in range(7)])
        angles.extend([[1,*chain[j:j+3]] for j in range(6)])
        torsions.extend([[1,*chain[j:j+4]] for j in range(5)])
    for z in [6,10]:
        for x in [8,16,24,32]:
            for y in [8,16,24,32]:
                i=atom(0,2,0,x,y,z)
                if z==6:anchor.append(i)
    for j in range(8):atom(100+j,3,0,6+4*j,34,35)
    atom(201,4,.1,8,8,40);atom(202,4,.1,28,8,40)
    atom(203,5,-.1,10,10,42);atom(204,5,-.1,30,10,42)
    txt=f'TOY ONLY - not PVDF\n\n{len(atoms)} atoms\n{len(bonds)} bonds\n{len(angles)} angles\n{len(torsions)} dihedrals\n0 impropers\n\n5 atom types\n1 bond types\n1 angle types\n1 dihedral types\n\n0 40 xlo xhi\n0 40 ylo yhi\n0 60 zlo zhi\n\nMasses\n\n1 12\n2 63\n3 18\n4 7\n5 19\n\nAtoms # full\n\n'+'\n'.join(atoms)+'\n'
    for name,rows in [('Bonds',bonds),('Angles',angles),('Dihedrals',torsions)]:
        txt+='\n'+name+'\n\n'+'\n'.join(str(i+1)+' '+' '.join(map(str,row)) for i,row in enumerate(rows))+'\n'
    (folder/'toy.data').write_text(txt)
    (folder/'styles.in').write_text('pair_style lj/cut/coul/long 12.0\nbond_style harmonic\nangle_style harmonic\ndihedral_style opls\nspecial_bonds lj/coul 0 0 0.5\npair_modify mix geometric\nkspace_style pppm 1e-5\n')
    (folder/'coeffs.in').write_text('pair_coeff * * 0.01 2.0\nbond_coeff 1 1.0 2.02\nangle_coeff 1 1.0 165\ndihedral_coeff 1 0.01 0.01 0.01 0.01\n')
    (folder/'groups.in').write_text('group polymer type 1\ngroup backbone type 1\ngroup substrate type 2\ngroup particles type 2\ngroup solvent type 3\ngroup Li type 4\ngroup anion type 5\ngroup donor_poly type 1\ngroup donor_solvent type 3\ngroup donor_anion type 5\ngroup probe molecule 1\ngroup anchor id '+' '.join(map(str,anchor))+'\ngroup grip id '+' '.join(map(str,grip))+'\nvariable Li_type index 4\nvariable poly_donor_type index 1\nvariable solv_donor_type index 3\nvariable anion_donor_type index 5\n')
    (folder/'hooks.in').write_text('# none\n')
    return folder


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--lmp',required=True);p.add_argument('--work-dir',required=True);p.add_argument('--generate-only',action='store_true');a=p.parse_args()
    work=Path(a.work_dir).resolve();folder=fixture(work/'model');results=[]
    jobs=[('equilibrate',None,{}),('adhesion','equilibrium',{}),('adhesion','pull_normal',{'pull_distance':.002}),('adhesion','pull_shear',{'pull_distance':.002}),('adhesion','umbrella',{'zcenter':17}),('mechanics','structure',{}),('mechanics','tensile',{'maxstrain':.000001}),('mechanics','shear',{'maxstrain':.000001}),('rheology','structure',{}),('rheology','adsorption',{}),('rheology','gk',{}),('rheology','nemd',{}),('transport',None,{}),('oxidation',None,{})]
    jobs += [('mechanics','elastic',{'component':i,'epsilon':eps}) for i in range(1,7) for eps in [.005,-.005]]
    executable=shutil.which(a.lmp) or a.lmp
    for j,(task,mode,extra) in enumerate(jobs):
        out=work/f'{j:02d}_{task}_{mode or "default"}';out.mkdir(exist_ok=True)
        vars={'task':task,'data':str(folder/'toy.data'),'styles':str(folder/'styles.in'),'coeffs':str(folder/'coeffs.in'),'groups':str(folder/'groups.in'),'hooks':str(folder/'hooks.in'),'out':str(out),'neq':10,'nprod':100,'nrelax':10,'sample':2,'dump_every':20,'ion_every':10,'corr_every':2,'corr_points':10,'dtfs':.1,'T':10}
        if mode:vars['mode']=mode
        vars.update(extra)
        cmd=[executable,'-in','lammps.in']
        for k,v in vars.items():cmd+=['-var',k,str(v)]
        if a.generate_only:
            results.append({'command':cmd,'status':'not_run'});continue
        r=subprocess.run(cmd,cwd=ROOT,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=180)
        (out/'console.txt').write_text(r.stdout);results.append({'task':task,'mode':mode,'extra':extra,'returncode':r.returncode})
        print(task,mode,extra,'PASS' if r.returncode==0 else 'FAIL')
    (work/'smoke_results.json').write_text(json.dumps(results,indent=2)+'\n')
    if not a.generate_only and any(r['returncode'] for r in results):raise SystemExit(1)

if __name__=='__main__':main()
