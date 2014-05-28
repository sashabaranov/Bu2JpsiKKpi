from math import sqrt

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
all_Kpipi = db['Kpipi']['RD'].values() + db['Kpipi']['MC'].values()


def perform_fit(model, ds, nbin):
    model.s.setMax(1.2 * len(ds))
    ru, fu = model.fitTo(ds, draw=True, nbins=nbin)

    model.signal.mean.release()
    ru, fu = model.fitTo(ds, draw=True, nbins=nbin)

    model.signal.sigma.release()
    ru, fu = model.fitTo(ds, draw=True, nbins=nbin)

    model.legend.AddEntry("", "Nsig = " + str(ru("SBu")[0]), "")
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

print len(all_KKK) + len(all_Kpipi)

var_names = ["BBu", "S2Bu", "S3Bu", "SBu", "mean_Bu1", "phi1_BBu", "phi2_BBu", "phi3_BBu", "sigma_Bu1", "tau_BBu"]
ntuple = ROOT.TNtuple("ntuple","ntuple", ":".join(var_names))


for kkk_hist in all_KKK:
    for kpipi_hist in all_Kpipi:
        s1_Bu = Models.CB2_pdf(
            'Bu1',
            m_Bu.getMin(),
            m_Bu.getMax(),
            fixMass=5.2792e+0,
            fixSigma=0.008499e+0,
            fixAlphaL=2.1018e+00,
            fixAlphaR=1.9818e+00,
            fixNL=6.1464e-01,
            fixNR=2.1291e+00,
            mass=m_Bu
        )

        kkk = Models.H1D_pdf(name="KKK", mass=m_Bu, histo=smear_kkk(kkk_hist))
        kpipi = Models.H1D_pdf(name="Kpipi", mass=m_Bu, histo=smear_kpipi(kpipi_hist))

        model_Bu = Charm3_pdf(
            signal=s1_Bu,
            signal2=kkk,
            signal3=kpipi,
            background=Models.Bkg_pdf('BBu', mass=m_Bu, power=3), suffix='Bu'
        )

        ru, fu = perform_fit(model_Bu, ds_Bu, nbin_Bu)

        ntuple.Fill(*[float(ru(v)[0]) for v in var_names])

        canvas.SaveAs('pics/{}_{}.png'.format(kkk_hist.GetName(), kpipi_hist.GetName()))

ntuple.SaveAs("fits.root")
