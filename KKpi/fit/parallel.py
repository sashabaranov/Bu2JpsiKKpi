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
var_names = ["phi1_BBu", "BBu", "S2Bu", "SBu", "mean_Bu1", "sigma_Bu1", "tau_BBu", "covqual", "minNll", "apolonios", "power", "i", "error"]




ROOT.gROOT.SetBatch (True)

canvas_size = 1024, 768




class Numerator(object):
    def __init__(self):
        self.n = 0

    def __iter__(self):
        return self

    # Python 3 compatibility
    def __next__(self):
        return self.next()

    def next(self):
        self.n += 1
        return self.n

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


def perform_fit_apolonios(i, model, ds, nbin):
    model.s.setMax(1.2 * len(ds))
    ru, fu = model.fitTo(ds, draw=True, nbins=nbin)

    model.signal.mean.release()
    ru, fu = model.fitTo(ds, draw=True, nbins=nbin)

    model.signal.sigma.release()
    ru, fu = model.fitTo(ds, draw=True, nbins=nbin)
 
    # model.signal.alpha.release()
    # ru, fu = model.fitTo(ds_Bu, draw=True, nbins=nbin_Bu)

    # model.signal.n.release()
    # ru, fu = model.fitTo(ds_Bu, draw=True, nbins=nbin_Bu)


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
            fixMass=5.28,
            fixSigma=5.75e-3,
            fixAlphaL=2.02,
            fixAlphaR=1.97,
            fixNL=0.83,
            fixNR=2.56,
            mass=m_Bu
        )
    else:
        s1_Bu = Models.Apolonios_pdf(
            'Bu1',
            m_Bu.getMin(),
            m_Bu.getMax(),
            fixMass=5.2792e+0,
            fixSigma=0.008499e+0,
            fixAlpha=9.52,
            fixN=19.9,
            fixB=1.0,
            mass=m_Bu
        )


    bkg_Bu = Models.Bkg_pdf('BBu', mass=m_Bu, power=power)
    kpipi = Models.H1D_pdf(name=kpipi_hist.GetName(), mass=m_Bu, histo=smear_kpipi(kpipi_hist))
    
    model_B = Charm3_pdf(
        signal=s1_Bu,
        signal2=kpipi,
        background=bkg_Bu,
        suffix="Bu"
    )

    bkg_Bu.tau.setMax(0.0)
    bkg_Bu.tau.setMin(-10.0)

    model_B.b.setMax(1800)
    # model_Bu.background.tau.setMax(-2.0)
    # model_Bu.background.tau.setVal(-1.0)
    c = ROOT.TCanvas(str(uuid.uuid4()), '', *canvas_size)

    if apolonios:
        ru, fu = perform_fit_apolonios(i, model_B, ds_Bu, nbin_Bu)
    else:
        ru, fu = perform_fit(i, model_B, ds_Bu, nbin_Bu)

    additional = [ru.covQual(), ru.minNll(), int(apolonios), power, bkg_Bu.phis[-1], i, ru("SBu")[0].error()]
    output = [float(ru(v)[0]) for v in var_names[:-len(additional)]] + additional

    fu.Draw()
    model_B.legend.Draw()
    prefix = "apol" if apolonios else "cb2"
    c.SaveAs('pics/{}_{}{}.png'.format(prefix, i, ":".join([str(x) for x in output])))


    return output




def star(x):
    return worker(*x)

def main():
    global all_Kpipi, ntuple

    numerate = Numerator()
    jobs = []
    for i, kpipi_hist in enumerate(all_Kpipi):
        for power in (1,2):
            jobs.append((numerate.next(), kpipi_hist, False, power))
            # jobs.append((numerate.next(), kpipi_hist, True, power))


    # Parallel:
    # workers  = cpu_count()
    # p = Pool(processes=workers)
    # results = p.map(star, jobs)
    
    # Iterative:
    results = []
    for j in jobs:
        results.append(star(j))

    ntuple = ROOT.TNtuple("ntuple","ntuple", ":".join(var_names))


    for r in results:
        ntuple.Fill(*r)

    ntuple.SaveAs("fits.root")




main()
