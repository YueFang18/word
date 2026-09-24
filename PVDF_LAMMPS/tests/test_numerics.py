import sys,unittest,tempfile,json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'analysis'));sys.path.insert(0,str(ROOT/'batch'))
from postprocess import fit,acf_fft,lag_msd,running_integral,wham,R,KB,QE,dumps,load_transport
from generate import descriptors
from association import contact_survival

class NumericalTests(unittest.TestCase):
    def test_regression_modulus(self):
        e=np.linspace(0,.02,21);r=fit(e,1750*e+2,.002,.01)
        self.assertAlmostEqual(r['slope'],1750);self.assertAlmostEqual(r['intercept'],2)

    def test_units_work(self):
        d=np.linspace(0,10,101);stress=np.full_like(d,100)
        self.assertAlmostEqual(float(np.trapezoid(stress,d)*1e-4),.1)
        # kcal/mol/A2 and kcal/mol/A3 conversion.
        self.assertAlmostEqual(4184/6.02214076e23/1e-20,.6947695457,places=8)
        self.assertAlmostEqual(4184/6.02214076e23/1e-30/1e6,6947.695457,places=5)

    def test_fft_against_direct(self):
        x=np.random.default_rng(12).normal(size=(123,3))+7
        c=acf_fft(x,30);y=x-x.mean(0)
        direct=np.array([(y[:len(y)-k]*y[k:]).mean(0) for k in range(31)])
        np.testing.assert_allclose(c,direct,atol=1e-13)

    def test_exponential_integral(self):
        t=np.arange(0,20,.001);y=np.exp(-t/2)
        self.assertAlmostEqual(running_integral(y,.001)[-1],2,places=3)

    def test_multi_origin_msd(self):
        # Opposite deterministic translations: exact squared displacement and cross correlations.
        t=np.arange(20.);xyz=np.zeros((20,2,3));xyz[:,0,0]=t;xyz[:,1,0]=-t
        lags=np.arange(1,7);np.testing.assert_allclose(lag_msd(xyz,lags),lags**2)
        M=xyz[:,0]-xyz[:,1]
        np.testing.assert_allclose(lag_msd(M[:,None,:],lags),4*lags**2)

    def test_diffusion_and_conductivity_units(self):
        # D=1e-10 m2/s -> slope 6e-5 A2/fs.
        self.assertAlmostEqual(6e-5*1e-5/6,1e-10)
        # 2 singly-charged species each 1 ion in 1000 A3 and D=1e-10.
        ne=2*QE**2*1e-10/(1000e-30*KB*300)
        slope=2*6e-5
        collective=slope*QE**2*1e25/(6*1000*KB*300)
        self.assertAlmostEqual(ne,collective)

    def test_wham_known_harmonic(self):
        z=np.linspace(-3,3,121);centers=np.linspace(-2,2,9);k=np.full(9,2.0);T=300
        physical=.5*.6*z*z
        weights=np.exp(-(physical[None,:]+.5*k[:,None]*(z[None,:]-centers[:,None])**2)/(R*T))
        # Deterministic exact histograms eliminate random estimation noise.
        hist=1e6*weights/weights.sum(axis=1)[:,None]
        zz,g,_=wham(hist,centers,k,z,T)
        target=.5*.6*zz*zz;target-=target.min()
        np.testing.assert_allclose(g,target,atol=1e-7)

    def test_wham_reject_disconnected(self):
        with self.assertRaisesRegex(ValueError,'Disconnected'):
            wham(np.array([[2,2,0,0],[0,0,2,2]]),[0,3],[1,1],np.arange(4.),300)

    def test_sequence_matched_composition(self):
        random=descriptors(list('ABAABABA'));block=descriptors(list('AAAAABBB'))
        self.assertEqual(random['monomer_fractions'],block['monomer_fractions'])
        self.assertGreater(random['heterodyad_fraction'],block['heterodyad_fraction'])
        self.assertGreater(block['max_run_length'],random['max_run_length'])

    def test_contact_survival(self):
        h=np.array([1,1,0,1,1],dtype=bool)[:,None,None]
        s=contact_survival(h,[1,2])
        self.assertAlmostEqual(s[0,1],2/3);self.assertAlmostEqual(s[0,2],2/3)
        self.assertAlmostEqual(s[1,1],1/2);self.assertAlmostEqual(s[1,2],0)

    def test_dump_parser_sort_and_box(self):
        text='ITEM: TIMESTEP\n20\nITEM: NUMBER OF ATOMS\n2\nITEM: BOX BOUNDS pp pp pp\n0 10\n0 20\n0 30\nITEM: ATOMS id mol type q mass xu yu zu\n2 2 1 -1 5 1 2 3\n1 1 1 1 5 4 5 6\n'
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'dump';p.write_text(text);step,cell,d=next(dumps(p))
            self.assertEqual(step,20);self.assertAlmostEqual(np.linalg.det(cell),6000)
            self.assertEqual(list(d['id']),[1,2]);self.assertEqual(d['xu'][0],4)

    def test_drift_not_removed_per_species(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp=Path(tmp);alltxt=[];li='';an=''
            for t in range(24):
                header=f'ITEM: TIMESTEP\n{t*10}\nITEM: NUMBER OF ATOMS\n'
                tail='ITEM: BOX BOUNDS pp pp pp\n0 100\n0 100\n0 100\nITEM: ATOMS id mol type q mass xu yu zu\n'
                # common drift 10t plus opposite species motion t and -t.
                rows=[f'1 1 1 1 1 {11*t} 0 0\n',f'2 2 2 -1 1 {9*t} 0 0\n']
                alltxt.append(header+'2\n'+tail+''.join(rows))
                if not t:li=header+'1\n'+tail+rows[0];an=header+'1\n'+tail+rows[1]
            for name,text in [('all',''.join(alltxt)),('li',li),('an',an)]: (tmp/name).write_text(text)
            steps,xyz,M,V,nli,q=load_transport(tmp/'all',tmp/'li',tmp/'an')
            np.testing.assert_allclose(xyz[:,0,0],np.arange(24))
            np.testing.assert_allclose(xyz[:,1,0],-np.arange(24))
            np.testing.assert_allclose(M[:,0],2*np.arange(24))

if __name__=='__main__':unittest.main()
