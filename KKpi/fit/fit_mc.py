from tools import *


from variables import m_Bu, nbin_Bu
from cuts import cuts_Bu, prntCuts, mctrue
from histograms import h1, h2


from model import model_Bu_mc as model_Bu
from data import mc_Pythia6, mc_Pythia8, mc_total


cuts_Bu += mctrue
tBu = mc_total.data

logger.info('DATA chain name is %s ' % (tBu.GetName()))


for i in prntCuts(cuts_Bu, "  CUTS B+  "):
    logger.info(i)


# logger.info('Fill control B+  histogram (takes some time)')
# with timing():
#     tBu.Project(h1.GetName(), 'DTFm_b', cuts_Bu)

# model_Bu.b.fix(0)

# with rooSilent():
#     logger.info('Fit Bc+ & B+ histogram (check the model)')
#     r, f = model_Bu.fitHisto(h1)


from PyPAW.Selectors import SelectorWithVars
from variables import selector_variables

sel_Bu = SelectorWithVars(
    variables=selector_variables,
    selection=cuts_Bu
)

logger.info(
    'Build RooFit dataset for B+ , it could take as long as 3-5 minutes')

tBu.process(sel_Bu)

ds_Bu = sel_Bu.dataset()
ds_Bu.Print('v')



logger.info('Make unbinned fit for B+')
#
model_Bu.s.setMax(1.2 * len(ds_Bu))
ru, fu = model_Bu.fitTo(ds_Bu, draw=True, nbins=nbin_Bu)

model_Bu.signal.sigma.release()
model_Bu.signal.mean.release()
ru, fu = model_Bu.fitTo(ds_Bu, draw=True, nbins=nbin_Bu)

model_Bu.signal.aR.release()
model_Bu.signal.aL.release()
ru, fu = model_Bu.fitTo(ds_Bu, draw=True, nbins=nbin_Bu)

model_Bu.signal.nR.release()
model_Bu.signal.nL.release()

ru, fu = model_Bu.fitTo(ds_Bu, draw=True, nbins=nbin_Bu)



# logger.info('running sPlot')
# model_Bu.sPlot(ds_Bu)





# ===============================================
# Drawing
# ===============================================


# h_1 = ROOT.TH1F("h_1", '',  30, 5.16, 5.45)
# h_1.Sumw2()
# h_2 = ROOT.TH1F("h_2", '',  30, 5.16, 5.45)
# h_2.Sumw2()


# h_1.SetXTitle('#Inv.\,mass(J/\psi\,K^+\,K^-\,K^+) \,\, with \,\, misid, GeV/c^2')
# h_2.SetXTitle('#Inv.\,mass(J/\psi\,K^+\,K^-\,K^+) \,\, with \,\, misid, GeV/c^2')

# ds_Bu.project(h_1, "mass_k1aspi", "SBu_sw && ann_kaon_PI_2 > 0.1")# && ann_kaon0 > 0.3")
# ds_Bu.project(h_2, "mass_k3aspi", "SBu_sw && ann_kaon_PI_0 > 0.1")# && ann_kaon2 > 0.3")

# h_1.red()
# h_2.blue()

# for j in xrange(0, h_1.GetNbinsX()):
#     if h_1.GetBinContent(j) < 0:
#         h_1.SetBinContent(j, 0)

# for j in xrange(0, h_2.GetNbinsX()):
#     if h_2.GetBinContent(j) < 0:
#         h_2.SetBinContent(j, 0)


# h_1.Scale(1.0/h_1.Integral())
# h_2.Scale(1.0/h_2.Integral())



# h_1.Draw("")
# h_2.Draw("same")
