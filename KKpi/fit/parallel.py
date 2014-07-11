from multiprocessing import Pool, cpu_count

from math import sqrt
import base64

from tools import *
from variables import *
from cuts import *
from data import selection7
from model import *
import uuid

logger = getLogger(__name__)

from Ostap.Selectors import SelectorWithVarsCached

sel_Bu = SelectorWithVarsCached(
    variables=selector_variables,
    selection=cuts_Bu,
    files=selection7.files
)

if not sel_Bu._loaded_from_cache:
    selection7.chain.process(sel_Bu)

ds_Bu = sel_Bu.dataset()
ds_Bu.Print('v')


# all_KKK = db['KKK']['MC'].values() + db['KKK']['RD'].values() 
all_Kpipi = db['Kpipi']['MC'].values() + db['Kpipi']['RD'].values()

#var_names = ["BBu", "S2Bu", "S3Bu", "SBu", "mean_Bu1", "phi1_BBu", "sigma_Bu1", "tau_BBu", "covqual"]
var_names = ["BBu", "S2Bu", "SBu", "mean_Bu1", "sigma_Bu1", "tau_BBu", "covqual", "minNll", "apolonios"]




ROOT.gROOT.SetBatch (True)

canvas_size = 1024, 768



def count_significance(model_Bu, ds_Bu, nbin_Bu):
    with RooSilent():
        try:
            model_Bu.s.setMax(1.2 * len(ds_Bu))

            ru, fu = model_Bu.fitTo(ds_Bu, draw=False, nbins=nbin_Bu)

            model_Bu.signal.sigma.release()
            model_Bu.signal.mean.release()


            ru, fu = model_Bu.fitTo(ds_Bu, draw=False, nbins=nbin_Bu)
            Sfixed = str(ru("SBu")[0])

            vals = [
                model_Bu.b,
                model_Bu.s,
                model_Bu.signal2,
                model_Bu.signal3,
                model_Bu.background.tau,
                model_Bu.signal.mean,
                model_Bu.signal.sigma,
            ] + model_Bu.background.phis

            for x in vals:
                x.fix(x.getVal())

            ru, fu = model_Bu.fitTo(ds_Bu, draw=False, nbins=nbin_Bu)
            Lfixed = ru.minNll()

            model_Bu.s.fix(0)

            ru, fu = model_Bu.fitTo(ds_Bu, draw=False, nbins=nbin_Bu)
            Ls0 = ru.minNll()
            return {'significance': sqrt(2 * (Ls0 - Lfixed)), 'Nsig': Sfixed}
        except ValueError:
            return -1



def perform_fit(i, model, ds, nbin):
    model.s.setMax(1.2 * len(ds))
    ru, fu = model.fitTo(ds, draw=True, nbins=nbin)

    model.signal.mean.release()
    ru, fu = model.fitTo(ds, draw=True, nbins=nbin)

    model.signal.sigma.release()
    ru, fu = model.fitTo(ds, draw=True, nbins=nbin)

    model.legend.AddEntry("", "N_{sig} \, = \," + str(ru("SBu")[0]), "")
    model.legend.Draw()

    fu.SetTitle("")
    
    fu.SetXTitle('Inv.\,mass(J/\psi KK\pi), GeV/c^{2}}')
    fu.SetYTitle("Events / %.1f MeV/c^{2}" % binning_b)

    fu.Draw()

    model.legend.Draw()
    
    return ru, fu





def worker(i, kpipi_hist, apolonios, power):
    global ntuple

    print "wolverine ", kpipi_hist

    if not apolonios:
        s1_Bu = Models.CB2_pdf(
            'Bu1',
            m_Bu.getMin(),
            m_Bu.getMax(),
            fixMass=5.2792e+0,
            fixSigma=0.006,
            fixAlphaL=2.1018e+00,
            fixAlphaR=1.9818e+00,
            fixNL=6.1464e-01,
            fixNR=2.1291e+00,
            mass=m_Bu
        )
    else:
        s1_Bu = Models.Apolonios_pdf(
            'Bu1',
            m_Bu.getMin(),
            m_Bu.getMax(),
            fixMass=5.2792e+0,
            fixSigma=0.008499e+0,
            fixAlpha=3.6,
            fixN=0.08,
            fixB=1.369,
            mass=m_Bu
        )
        s1_Bu.alpha.relese()
        s1_Bu.sigma.release()
        s1_Bu.n.release()
        s1_Bu.b.release()


    bkg_Bu = Models.Bkg_pdf('BBu', mass=m_Bu, power=power)
    kpipi = Models.H1D_pdf(name=kpipi_hist.GetName(), mass=m_Bu, histo=smear_kpipi(kpipi_hist))
    
    model_B = Charm3_pdf(
        signal=s1_Bu,
        signal2=kpipi,
        background=bkg_Bu,
        suffix="Bu"
    )

    bkg_Bu.tau.setMax(0.0)

    # model_Bu.background.tau.setMax(-2.0)
    # model_Bu.background.tau.setVal(-1.0)
    c = ROOT.TCanvas(str(uuid.uuid4()), '', *canvas_size)

    ru,fu = perform_fit(i, model_B, ds_Bu, nbin_Bu)

    additional = [ru.covQual(), ru.minNll(), int(apolonios), power]
    output = [float(ru(v)[0]) for v in var_names[:-len(additional)]] + additional

    fu.Draw()
    model_B.legend.Draw()

    c.SaveAs('pics/{}{}.png'.format(i, base64.b64encode(":".join([str(x) for x in output]))))


    return output




def star(x):
    return worker(*x)

def main():
    global all_Kpipi, ntuple


    jobs = []

    for i, kpipi_hist in enumerate(all_Kpipi):
        for power in (0, 1, 2, 3):
            jobs.append((i, kpipi_hist, False, power))
            jobs.append((i, kpipi_hist, True, power))


    # Parallel:
    workers  = cpu_count()
    p = Pool(processes=workers)
    results = p.map(star, jobs)
    
    # Iterative:
    # results = []
    # for j in jobs:
    #     results.append(star(j))

    ntuple = ROOT.TNtuple("ntuple","ntuple", ":".join(var_names))


    for r in results:
        ntuple.Fill(*r)

    ntuple.SaveAs("fits.root")




main()
