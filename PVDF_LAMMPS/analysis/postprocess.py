#!/usr/bin/env python3
"""PVDF protocol postprocessing. Python >=3.10, numpy. No scipy required.
Columns are zero-based including LAMMPS timestep column. See --help and README.
"""
import argparse
import json
from pathlib import Path
import numpy as np

KB = 1.380649e-23
QE = 1.602176634e-19
R = 0.00198720425864083  # kcal/mol/K


def table(path):
    a = np.loadtxt(path, comments='#', ndmin=2)
    if not np.all(np.isfinite(a)):
        raise ValueError(f'Nonfinite values in {path}')
    return a


def report(obj, path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False)+'\n')
    print(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False))


def block_stats(x, blocks=5):
    if blocks < 2 or len(x) < 2*blocks:
        raise ValueError('Need >=2 blocks and >=2 samples per block')
    chunks = np.array_split(np.asarray(x), blocks)
    means = np.array([np.mean(c, axis=0) for c in chunks])
    return {'mean': np.mean(x, axis=0).tolist(),
            'block_sem': (means.std(axis=0, ddof=1)/np.sqrt(blocks)).tolist(),
            'block_means': means.tolist(), 'n': len(x),
            'caution': 'SEM valid only if blocks exceed the relevant correlation time; increase block length.'}


def fit(x, y, lo, hi):
    mask = (x >= lo) & (x <= hi)
    if mask.sum() < 4 or hi <= lo:
        raise ValueError('Select a fit window with at least 4 points')
    xx, yy = x[mask], y[mask]
    if np.ptp(xx) <= 0:
        raise ValueError('Fit x has no range')
    m,b = np.polyfit(xx, yy, 1)
    residual = np.sum((yy-(m*xx+b))**2)
    denom = np.sum((yy-yy.mean())**2)
    return {'slope': float(m), 'intercept': float(b),
            'R2': float(1-residual/denom) if denom > 0 else None,
            'n': int(mask.sum()), 'fit_range': [lo,hi]}


def running_integral(y, dt):
    y = np.asarray(y)
    return np.concatenate((np.zeros((1,)+y.shape[1:]),
                           np.cumsum(0.5*(y[1:]+y[:-1])*dt, axis=0)), axis=0)


def acf_fft(x, maxlag):
    x = np.asarray(x, float)
    x = x - x.mean(axis=0)
    n = len(x)
    if maxlag >= n or maxlag < 1:
        raise ValueError('Invalid correlation lag')
    size = 1 << (2*n-1).bit_length()
    f = np.fft.rfft(x, n=size, axis=0)
    c = np.fft.irfft(f*np.conj(f), n=size, axis=0)[:maxlag+1]
    return c/np.arange(n,n-maxlag-1,-1).reshape((-1,)+(1,)*(x.ndim-1))


def lag_msd(coords, lags):
    # coords shape frame, particle, xyz; time origins all available per lag.
    return np.array([np.mean(np.sum((coords[k:]-coords[:-k])**2, axis=-1)) for k in lags])


def dumps(path):
    """Read standard orthogonal/restricted triclinic LAMMPS custom dump."""
    with open(path) as f:
        while True:
            line = f.readline()
            if not line:
                return
            if line.strip() != 'ITEM: TIMESTEP':
                raise ValueError('Expected ITEM: TIMESTEP')
            step = int(f.readline())
            if f.readline().strip() != 'ITEM: NUMBER OF ATOMS':
                raise ValueError('Expected atom dump, not local dump')
            n = int(f.readline())
            bh = f.readline().split()
            bounds = np.array([[float(z) for z in f.readline().split()] for _ in range(3)])
            if 'xy' in bh:
                xb,yb,zb = bounds
                xy,xz,yz = xb[2],yb[2],zb[2]
                xlo = xb[0]-min(0,xy,xz,xy+xz); xhi = xb[1]-max(0,xy,xz,xy+xz)
                ylo = yb[0]-min(0,yz); yhi = yb[1]-max(0,yz)
                cell = np.array([[xhi-xlo,0,0],[xy,yhi-ylo,0],[xz,yz,zb[1]-zb[0]]])
            else:
                cell = np.diag(bounds[:,1]-bounds[:,0])
            head = f.readline().split()
            if head[:2] != ['ITEM:','ATOMS']:
                raise ValueError('Expected custom ATOMS header')
            cols = head[2:]
            a = np.array([[float(z) for z in f.readline().split()] for _ in range(n)])
            a = a[np.argsort(a[:,cols.index('id')])]
            yield step,cell,{k:a[:,i] for i,k in enumerate(cols)}


def ids(path):
    return set(next(dumps(path))[2]['id'].astype(int))


def load_transport(path, li_ids, anion_ids, stride=1, discard=0):
    xyzs=[]; moments=[]; times=[]; cells=[]; ref=None; q=None; mass=None
    li_set, an_set = ids(li_ids), ids(anion_ids)
    for i,(step,cell,d) in enumerate(dumps(path)):
        if i < discard or (i-discard)%stride:
            continue
        atomids=d['id'].astype(int)
        if ref is None:
            ref=atomids
            q=d['q']; mass=d['mass']; mol=d['mol'].astype(int)
            if abs(q.sum())>1e-4:
                raise ValueError('Full trajectory must have zero net charge')
            li=np.array([k in li_set for k in atomids]); an=np.array([k in an_set for k in atomids])
            if not li.any() or not an.any() or (li&an).any():
                raise ValueError('Invalid Li/anion ID selections')
            if not li_set.issubset(set(atomids)) or not an_set.issubset(set(atomids)):
                raise ValueError('Selection IDs absent from trajectory')
            amols=np.unique(mol[an])
            if (amols<=0).any():
                raise ValueError('Each complete anion needs unique positive molecule ID')
            an_groups=[]
            for m in amols:
                mask=mol==m
                if not np.all(an[mask]):
                    raise ValueError(f'Anion molecule {m} is incomplete in anion group')
                an_groups.append(mask)
            ionq=np.r_[q[li], [q[g].sum() for g in an_groups]]
            nli=int(li.sum())
            if np.any(ionq[:nli] <= 0) or np.any(ionq[nli:] >= 0):
                raise ValueError('Li charges must be positive; complete anion charges negative')
            remaining=~(li|an)
            for m in np.unique(mol[remaining]):
                if abs(q[remaining & (mol==m)].sum())>1e-4:
                    raise ValueError('Unselected charged molecules: extend mobile-species analysis before comparing NE and collective conductivity')
        elif not np.array_equal(ref,atomids) or not np.allclose(d['q'],q,atol=1e-10,rtol=0):
            raise ValueError('Atom IDs/charges change; fixed-charge constant-N analysis only')
        xyz=np.column_stack([d[k] for k in ['xu','yu','zu']])
        # Common SYSTEM COM frame for self diffusion; neutral total charge moment unchanged.
        xyz-=np.average(xyz,axis=0,weights=mass)
        ions=np.vstack([xyz[li], [np.average(xyz[g],axis=0,weights=mass[g]) for g in an_groups]])
        xyzs.append(ions); moments.append((q[:,None]*xyz).sum(axis=0)); times.append(step); cells.append(cell)
    if len(times)<20:
        raise ValueError('Need >=20 retained frames')
    cells=np.array(cells)
    if not np.allclose(cells,cells[0],rtol=1e-7,atol=1e-7):
        raise ValueError('Use constant-cell production for transport analysis')
    return np.array(times),np.array(xyzs),np.array(moments),float(np.linalg.det(cells[0])),nli,ionq


def cmd_transport(a):
    steps,xyz,M,V,nli,q=load_transport(a.input,a.li_ids,a.anion_ids,a.stride,a.discard)
    diff=np.diff(steps)
    if not np.all(diff==diff[0]):
        raise ValueError('Nonuniform trajectory sampling')
    dt=diff[0]*a.dt_fs
    maxlag=min(len(steps)//3,int(a.maxlag_fs/dt))
    lags=np.unique(np.linspace(1,maxlag,min(maxlag,a.nlag)).astype(int))
    if len(lags)<4:
        raise ValueError('Insufficient lag range')
    t=lags*dt
    msdli=lag_msd(xyz[:,:nli],lags); msdan=lag_msd(xyz[:,nli:],lags)
    msdq=lag_msd(M[:,None,:],lags)
    # Individual effective charges are retained, including charge-scaled force fields.
    selfq=np.array([np.mean(np.sum((xyz[k:]-xyz[:-k])**2,axis=-1),axis=0)@q**2 for k in lags])
    fli=fit(t,msdli,a.fit_min,a.fit_max); fan=fit(t,msdan,a.fit_min,a.fit_max)
    fc=fit(t,msdq,a.fit_min,a.fit_max); fn=fit(t,selfq,a.fit_min,a.fit_max)
    factor=QE**2*1e25/(6*V*KB*a.temperature)
    sigma=fc['slope']*factor; ne=fn['slope']*factor
    np.savetxt(a.curve,np.column_stack([t,msdli,msdan,msdq,selfq]),
               header='lag_fs Li_MSD_A2 anion_COM_MSD_A2 collective_e2A2 NE_self_e2A2')
    report({'D_Li_m2_s':fli['slope']*1e-5/6,'D_anion_m2_s':fan['slope']*1e-5/6,
            'sigma_collective_S_m':sigma,'sigma_NE_effective_charges_S_m':ne,
            'sigma_over_NE':sigma/ne if ne>0 else None,'volume_A3':V,
            'Li_fit':fli,'anion_fit':fan,'collective_fit':fc,
            'warning':'Check diffusive regime (log slope near 1), positive stable slopes, independent replicas and fit-window sensitivity. Bound polarization causes short-time curvature. Negative conductivity slope is unconverged, not physical. Molecular image flags must be consistent. Charge-scaled results are model conductivities; NE uses the same effective charges.'},a.output)


def cmd_gk(a):
    if a.blocks < 2:
        raise ValueError('Need at least 2 statistical blocks')
    d=table(a.input)[a.discard:]
    t=d[:,1]; x=d[:,2:5]; V=d[:,5]
    dt=np.diff(t)
    if not np.allclose(dt,dt[0]) or not np.allclose(V,V.mean(),rtol=1e-7):
        raise ValueError('GK needs uniform times and constant volume')
    lag=int(a.maxlag_fs/dt[0])
    chunks=np.array_split(x,a.blocks)
    if lag<1 or min(map(len,chunks))<5*(lag+1):
        raise ValueError('Each block must be at least 5 correlation windows long; reduce maxlag or lengthen data')
    pref=V.mean()*1e-30/(KB*a.temperature)*101325.0**2*1e-15
    curves=np.array([running_integral(acf_fft(c,lag),dt[0])*pref for c in chunks])
    mean=curves.mean(axis=0)
    sem=curves.mean(axis=2).std(axis=0,ddof=1)/np.sqrt(a.blocks)
    tau=np.arange(lag+1)*dt[0]
    np.savetxt(a.curve,np.column_stack([tau,mean,mean.mean(axis=1),sem]),
               header='lag_fs eta_xy eta_xz eta_yz eta_mean_Pa_s block_sem_Pa_s')
    report({'end_integral_Pa_s':float(mean[-1].mean()),'end_block_sem_Pa_s':float(sem[-1]),
            'mean_pressure_components_atm':x.mean(axis=0).tolist(),'correlation_cutoff_fs':float(tau[-1]),
            'warning':'End integral is NOT automatically viscosity. Inspect plateau and cutoff/block/thermostat sensitivity. Isotropic averaging is appropriate only for isotropic equilibrium fluids; solid elastic plateaus do not define a finite zero-shear viscosity.'},a.output)


def cmd_traction(a):
    d=table(a.input)[a.discard:]; z=d[:,2]; stress=d[:,3 if a.mode=='normal' else 4]
    if len(z)<3 or np.any(np.diff(z)<=0):
        raise ValueError('Need monotonic loading displacement')
    # MPa * angstrom = 1e-4 J/m2
    value=float(np.trapezoid(stress,z)*1e-4)
    report({'peak_traction_MPa':float(stress.max()),'work_over_sampled_interval_J_m2':value,
            'interval_A':[float(z[0]),float(z[-1])],
            'warning':'Integral covers sampled interval only; confirm initial zero load and complete separation. Includes finite-rate/bulk dissipation; inspect failure position.'},a.output)


def cmd_modulus(a):
    d=table(a.input)[a.discard:]
    result=fit(d[:,a.strain_col],d[:,a.stress_col],a.fit_min,a.fit_max)
    result['modulus_MPa']=result.pop('slope')
    result['warning']='Finite-rate apparent modulus unless quasistatic/time-scale convergence demonstrated. Linear range must be selected physically.'
    report(result,a.output)


def cmd_elastic(a):
    jobs=json.loads(Path(a.manifest).read_text())
    C=np.empty((6,6))
    for j in range(1,7):
        plus=jobs[str(j)]['plus']; minus=jobs[str(j)]['minus']; eps=jobs[str(j)]['epsilon']
        if eps<=0:
            raise ValueError('epsilon must be positive magnitude')
        p=table(plus)[a.discard:,2:8]; m=table(minus)[a.discard:,2:8]
        if len(p)<4 or len(m)<4:
            raise ValueError('Too few stress samples')
        C[:,j-1]=(p.mean(axis=0)-m.mean(axis=0))/(2*eps)
    sym=(C+C.T)/2
    eig=np.linalg.eigvalsh(sym)
    result={'C_raw_MPa':C.tolist(),'C_symmetric_MPa':sym.tolist(),
            'relative_asymmetry':float(np.linalg.norm(C-C.T)/max(np.linalg.norm(C),1e-12)),
            'symmetric_eigenvalues_MPa':eig.tolist(),
            'warning':'Same reference cell and morphology required. Repeat +/-strain at multiple magnitudes and replicas. Large antisymmetry or negative eigenvalues indicate instability/statistical/model problems; do not hide them by symmetrization.'}
    if np.all(eig>0):
        S=np.linalg.inv(sym)
        result['Young_xyz_MPa']=(1/np.diag(S)[:3]).tolist()
        result['shear_yz_xz_xy_MPa']=(1/np.diag(S)[3:]).tolist()
        result['nu_xy']=-float(S[1,0]/S[0,0])
        result['nu_xz']=-float(S[2,0]/S[0,0])
    report(result,a.output)


def logsumexp(x,axis=None):
    mx=np.max(x,axis=axis,keepdims=True)
    y=mx+np.log(np.exp(x-mx).sum(axis=axis,keepdims=True))
    return np.squeeze(y,axis=axis)


def wham(hist,centers,k,z,T,tol=1e-9,maxiter=100000):
    hist=np.asarray(hist,float); N=hist.sum(axis=1); H=hist.sum(axis=0)
    if (N<=0).any():
        raise ValueError('Empty umbrella window')
    # Histogram graph connectivity is necessary, not sufficient, for adequate overlap.
    adjacent=(hist>0).astype(int)@(hist>0).astype(int).T > 0
    seen={0}
    for _ in N:
        seen|={j for i in list(seen) for j in np.where(adjacent[i])[0]}
    if len(seen)!=len(N):
        raise ValueError('Disconnected umbrella histograms; add/equilibrate windows')
    use=H>0; zz=z[use]; HH=H[use]
    u=0.5*np.asarray(k)[:,None]*(zz[None,:]-np.asarray(centers)[:,None])**2/(R*T)
    f=np.zeros(len(N))
    for iteration in range(maxiter):
        logp=np.log(HH)-logsumexp(np.log(N)[:,None]+f[:,None]-u,axis=0)
        logp-=logsumexp(logp)
        new=-logsumexp(logp[None,:]-u,axis=1); new-=new[0]
        if np.max(abs(new-f))<tol:
            pmf=-R*T*logp; pmf-=pmf.min()
            return zz,pmf,iteration+1
        f=new
    raise ValueError('WHAM did not converge')


def cmd_wham(a):
    windows=json.loads(Path(a.manifest).read_text())
    samples=[]; centers=[]; ks=[]
    for w in windows:
        z=table(w['file'])[w.get('discard',0)::w.get('stride',1),2]
        if not len(z):
            raise ValueError('No retained samples')
        samples.append(z); centers.append(w['center_A']); ks.append(w['K_kcal_mol_A2'])
    lo=min(map(np.min,samples)); hi=max(map(np.max,samples))
    edges=np.arange(lo,hi+a.bin_width,a.bin_width)
    H=np.array([np.histogram(s,edges)[0] for s in samples]); z=(edges[1:]+edges[:-1])/2
    z,G,n=wham(H,centers,ks,z,a.temperature)
    np.savetxt(a.curve,np.c_[z,G],header='COM_z_minus_anchor_A PMF_kcal_mol_arbitrary_zero')
    report({'iterations':n,'windows':len(windows),
            'warning':'1D PMF, not an area-normalized adhesion free energy or standard adsorption free energy. Match lateral/conformational restraints across windows, use decorrelated samples, check hidden slow modes, perform independent repeats/bootstrap. Histogram graph connectivity alone is insufficient.'},a.output)


def cmd_dihedral(a):
    rows=[]
    with open(a.input) as f:
        while True:
            line=f.readline()
            if not line: break
            if line.strip()!='ITEM: TIMESTEP':
                raise ValueError('Unexpected local dump format')
            step=int(f.readline()); f.readline(); n=int(f.readline())
            header=f.readline()
            if header.startswith('ITEM: BOX BOUNDS'):
                for _ in range(3):f.readline()
                header=f.readline()
            d=np.array([[float(x) for x in f.readline().split()] for _ in range(n)])
            if not n: raise ValueError('No backbone dihedrals')
            # angle is final column of common/chain_stats.in dump.
            phi=(d[:,-1]+180)%360-180
            trans=np.abs(np.abs(phi)-180)<=a.tolerance
            gauche=np.minimum(abs(phi-60),abs(phi+60))<=a.tolerance
            rows.append([step,trans.mean(),gauche.mean(),(~(trans|gauche)).mean()])
    np.savetxt(a.curve,rows,header='step trans_fraction gauche_fraction other_fraction')
    report({'mean_trans_fraction':float(np.mean(np.array(rows)[:,1])),
            'warning':'Verify convention against a known all-trans geometry. Torsion fractions are NOT crystal phase fractions. Copolymer torsions must be stratified by chemical motif; this report pools selected backbone torsions.'},a.output)


def main():
    p=argparse.ArgumentParser(description=__doc__); sp=p.add_subparsers(dest='cmd',required=True)
    for name in ['mean','traction','modulus','elastic','gk','transport','wham','dihedrals']:
        s=sp.add_parser(name);s.add_argument('--output',default=f'{name}_summary.json')
        if name not in ['elastic','wham']:s.add_argument('input')
        else:s.add_argument('manifest')
        if name in ['mean','traction','modulus','elastic','gk','transport']:s.add_argument('--discard',type=int,default=0)
        if name in ['gk','transport','wham','dihedrals']:s.add_argument('--curve',default=f'{name}_curve.dat')
        if name in ['gk','transport','wham']:s.add_argument('--temperature',type=float,required=True)
        if name in ['modulus','transport']:
            s.add_argument('--fit-min',type=float,required=True);s.add_argument('--fit-max',type=float,required=True)
        if name in ['mean','gk']:s.add_argument('--blocks',type=int,default=5)
        if name=='mean':s.add_argument('--column',type=int,required=True)
        if name=='traction':s.add_argument('--mode',choices=['normal','shear'],default='normal')
        if name=='modulus':
            s.add_argument('--strain-col',type=int,default=2);s.add_argument('--stress-col',type=int,default=5)
        if name=='gk':s.add_argument('--maxlag-fs',type=float,required=True)
        if name=='transport':
            s.add_argument('--li-ids',required=True);s.add_argument('--anion-ids',required=True)
            s.add_argument('--dt-fs',type=float,required=True);s.add_argument('--stride',type=int,default=1)
            s.add_argument('--maxlag-fs',type=float,required=True);s.add_argument('--nlag',type=int,default=150)
        if name=='wham':s.add_argument('--bin-width',type=float,default=0.2)
        if name=='dihedrals':s.add_argument('--tolerance',type=float,default=30)
    a=p.parse_args()
    if a.cmd=='mean':report(block_stats(table(a.input)[a.discard:,a.column],a.blocks),a.output)
    else:globals()['cmd_'+a.cmd](a)

if __name__=='__main__':main()
