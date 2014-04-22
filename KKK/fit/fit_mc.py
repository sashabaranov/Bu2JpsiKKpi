import ROOT
from AnalysisPython.PyRoot import *
from AnalysisPython.PyRoUts import *
from AnalysisPython.Utils import timing
from AnalysisPython.Utils import rooSilent
from AnalysisPython.Logger import getLogger
logger = getLogger(__name__)


from cuts import m_Bu, nbin_Bu
from cuts import cuts_Bu, prntCuts
from cuts import h1, h2


logger.info('Import models/PDFs')
from model import model_Bu


logger.info('Import data (for the first time it could take some time)')
from data import mc_Pythia6, mc_Pythia8, mc_total


cuts_s6 = " && MIPCHI2DV_k1 > 9 && MIPCHI2DV_k2 > 9 && MIPCHI2DV_k3 > 9"
cuts_Bu += cuts_s6


tBu = mc_Pythia6.data


logger.info('DATA chain name is %s ' % (tBu.GetName()))


for i in prntCuts(cuts_Bu, "  CUTS B+  "):
    logger.info(i)


logger.info('Fill control B+  histogram (takes some time)')
with timing():
    tBu.Project(h1.GetName(), 'DTFm_b', cuts_Bu)

model_Bu.b.fix(0)

with rooSilent():
    logger.info('Fit Bc+ & B+ histogram (check the model)')
    r, f = model_Bu.fitHisto(h1)


from PyPAW.Selectors import SelectorWithVars
from Variables import selector_variables

sel_Bu = SelectorWithVars(
    variables=selector_variables,
    selection=cuts_Bu
)

logger.info(
    'Build RooFit dataset for B+ , it could take as long as 3-5 minutes')

tBu.process(sel_Bu)

ds_Bu = sel_Bu.dataset()
ds_Bu.Print('v')

# frame = m_Bu.frame(nbin_Bu)

# ds_Bu.plotOn(frame)

# frame.SetMinimum(1e-5)
# frame.Draw()

logger.info('Make unbinned fit for B+')

model_Bu.s.setMax(1.2 * len(ds_Bu))
ru, fu = model_Bu.fitTo(ds_Bu, draw=True, nbins=nbin_Bu)

model_Bu.signal.sigma.release()
model_Bu.signal.mean.release()
ru, fu = model_Bu.fitTo(ds_Bu, nbins=nbin_Bu)

model_Bu.signal.aR.release()
model_Bu.signal.aL.release()
ru, fu = model_Bu.fitTo(ds_Bu, nbins=nbin_Bu)

model_Bu.signal.nR.release()
model_Bu.signal.nL.release()

ru, fu = model_Bu.fitTo(ds_Bu, draw=True, nbins=nbin_Bu)

# fu.SetXTitle('#Inv.\,mass(J/\psi\,K^+\,K^-\,K^+), GeV/c^2')
# fu.SetYTitle('')


# print 'FIT#2 results for B+  ', ru(model_Bu.s_name)[0]
# print 'FIT#2 precision:', ru("SBu")[0].prec()


# logger.info('running sPlot')

# model_Bu.sPlot(ds_Bu)

# h = ROOT.TH1F("h", '',  30, 5.16, 5.45)
# h.Sumw2()

# ds_Bu.project(h, "m_b_misid", "SBu_sw")
# # h.Draw()

# d = shelve.open('$KKpidir/fit/histos.shelve')

# d['KKK_30bin_sel6'] = h


# d.close()


# plotting against MC


  
logger.info('end of the module')
