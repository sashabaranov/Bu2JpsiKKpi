from multiprocessing import Pool, cpu_count

from math import sqrt
import base64

from tools import *
from variables import *
from cuts import *
from data import tSelection7 as tBu
from model import *


logger = getLogger(__name__)

sel_Bu = SelectorWithVars(
    variables=selector_variables,
    selection=cuts_Bu
)

tBu.process(sel_Bu)

ds_Bu = sel_Bu.dataset()
ds_Bu.Print('v')


all_KKK = db['KKK']['RD'].values() + db['KKK']['MC'].values()
all_Kpipi = db['Kpipi']['MC'].values() # + db['Kpipi']['RD'].values()

#var_names = ["BBu", "S2Bu", "S3Bu", "SBu", "mean_Bu1", "phi1_BBu", "sigma_Bu1", "tau_BBu", "covqual"]
var_names = ["BBu", "S2Bu", "S3Bu", "SBu", "mean_Bu1", "sigma_Bu1", "tau_BBu", "covqual", "minNll", "i", "correlation"]

def perform_fit(model, ds, nbin):
    model.s.setMax(1.2 * len(ds))
    ru, fu = model.fitTo(ds, draw=True, nbins=nbin)

    model.signal.mean.release()
    ru, fu = model.fitTo(ds, draw=True, nbins=nbin)

    model.signal.sigma.release()
    ru, fu = model.fitTo(ds, draw=True, nbins=nbin)

    model.legend.AddEntry("", "N_{sig} \, = \," + str(ru("SBu")[0]), "")
    model.legend.Draw()

    fu.Draw()
    model.legend.Draw()

    return ru, fu


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


def random_canvas(size=(1024, 768)):
    """Helper method for creating canvas with random name"""

    import uuid
    name = str(uuid.uuid4())

    # Check if icanvas already exists
    canvas = ROOT.gROOT.FindObject(name)
    assert len(size) == 2
    if canvas:
        raise Exception('uuid found')
    else:
        return ROOT.TCanvas(name, name, size[0], size[1])



def worker(i, kkk_hist, kpipi_hist):
    global ntuple

    c1 = random_canvas()
    s1_Bu = Models.CB2_pdf(
        'Bu1',
        m_Bu.getMin(),
        m_Bu.getMax(),
        fixMass=5.2792e+0,
        fixSigma=0.008499e+0,
        fixAlphaL=2.0266,
        fixAlphaR=2.00875,
        fixNL=0.839,
        fixNR=2.34195,
        mass=m_Bu
    )

    kkk = Models.H1D_pdf(name=kkk_hist.GetName(), mass=m_Bu, histo=smear_kkk(kkk_hist))
    kpipi = Models.H1D_pdf(name=kpipi_hist.GetName(), mass=m_Bu, histo=smear_kpipi(kpipi_hist))

    model_Bu = Charm3_pdf(
        signal=s1_Bu,
        signal2=kkk,
        signal3=kpipi,
        background=Models.Bkg_pdf('BBu', mass=m_Bu, power=0), suffix='Bu'
    )

    model_Bu.background.tau.setMax(-2.0)
    model_Bu.background.tau.setVal(-1.0)

    model_Bu.s.setMin(100)
    model_Bu.s2.setMax(220)
    model_Bu.s2.setVal(0.1)


    ru, fu = perform_fit(model_Bu, ds_Bu, nbin_Bu)

    additional = [ru.covQual(), ru.minNll(), i, ru.correlationMatrix()[1][6]]

    output = [float(ru(v)[0]) for v in var_names[:-len(additional)]] + additional

    c1.SaveAs('pics/{}.png'.format(base64.b64encode(":".join([str(x) for x in output]))))

    return output

def star(x):
    return worker(*x)

def main():
    global all_KKK, all_Kpipi, ntuple


    jobs = []

    i = 0
    for kkk_hist in all_KKK:
        for kpipi_hist in all_Kpipi:
            jobs.append((i, kkk_hist, kpipi_hist))
            i += 1

    workers  = cpu_count()
    p = Pool(processes=workers)
    results = p.map(star, jobs)


    ntuple = ROOT.TNtuple("ntuple","ntuple", ":".join(var_names))

    for r in results:
        ntuple.Fill(*r)

    ntuple.SaveAs("fits.root")




if __name__ == "__main__":
    main()
