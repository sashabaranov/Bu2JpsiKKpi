from tools import *

from Variables import m_Bu, nbin_Bu
from cuts import cuts_Bu, prntCuts, mctrue
from cuts import h1

cuts_Bu += mctrue


logger.info('Import models/PDFs')
from model import model_Bu


logger.info('Import data (for the first time it could take some time)')
from data import mc_Pythia6, mc_Pythia8, mc_total



tBu = mc_Pythia6.data


logger.info('DATA chain name is %s ' % (tBu.GetName()))


for i in prntCuts(cuts_Bu, "  CUTS B+  "):
    logger.info(i)


logger.info('Fill control B+  histogram (takes some time)')
with timing():
    tBu.Project(h1.GetName(), 'DTFm_b', cuts_Bu)

model_Bu.b.fix(0)
model_Bu.background.tau.fix(0)

with rooSilent():
    logger.info('Fit Bc+ & B+ histogram (check the model)')
    r, f = model_Bu.fitHisto(h1)


from PyPAW.Selectors import SelectorWithVars
from Variables import selector_variables

sel_Bu = SelectorWithVars(
    variables=selector_variables,
    selection=cuts_Bu
)

logger.info('Build RooFit dataset for B+ , it could take as long as 3-5 minutes')

tBu.process(sel_Bu)

ds_Bu = sel_Bu.dataset()
ds_Bu.Print('v')



logger.info('Make unbinned fit for B+')

model_Bu.s.setMax(1.2 * len(ds_Bu))
ru, fu = model_Bu.fitTo(ds_Bu, draw=True, nbins=nbin_Bu)

model_Bu.signal.sigma.release()
ru, fu = model_Bu.fitTo(ds_Bu, draw=True, nbins=nbin_Bu)

model_Bu.signal.mean.release()
ru, fu = model_Bu.fitTo(ds_Bu, draw=True, nbins=nbin_Bu)

model_Bu.signal.aR.release()
ru, fu = model_Bu.fitTo(ds_Bu, draw=True, nbins=nbin_Bu)

model_Bu.signal.aL.release()
ru, fu = model_Bu.fitTo(ds_Bu, draw=True, nbins=nbin_Bu)

model_Bu.signal.nR.release()
ru, fu = model_Bu.fitTo(ds_Bu, draw=True, nbins=nbin_Bu)

model_Bu.signal.nL.release()
ru, fu = model_Bu.fitTo(ds_Bu, draw=True, nbins=nbin_Bu)

logger.info('running sPlot')
model_Bu.sPlot(ds_Bu)


# print 'FIT#2 results for B+  ', ru(model_Bu.s_name)[0]
# print 'FIT#2 precision:', ru("SBu")[0].prec()



# hists = [
#     ("p6_pi1", "mass_pi1ask", "SBu_sw"),
#     ("p6_pi1_cuts", "mass_pi1ask", "SBu_sw && ann_pion_K > 0.1"),
# ]



# logger.info('Writing histos')
# db = shelve.open('$KKpidir/fit/histos.shelve')

# d = db['Kpipi']
# d['MC'] = {param[0]: make_hist(*param) for param in hists}
# db['Kpipi'] = d

# db.sync()
# db.close()



# h_1 = ROOT.TH1F("h_1", '',  30, 5.16, 5.45)
# h_1.Sumw2()

# h_1.SetXTitle('#Inv.\,mass(J/\psi\,K^+\,\pi^-\,\pi^+) \,\, with \,\, misid, GeV/c^2')

# ds_Bu.project(h_1, "mass_pi1ask", "SBu_sw") # && ann_kaon_PI_2 > 0.1")# && ann_kaon0 > 0.3")

# h_1.red()

# for j in xrange(0, h_1.GetNbinsX()):
#     if h_1.GetBinContent(j) < 0:
#         h_1.SetBinContent(j, 0)

# h_1.Scale(1.0/h_1.Integral())



# h_1.Draw("")
