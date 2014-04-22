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
from data import tSelection7 as tBu


cuts_s6 = " && MIPCHI2DV_k1 > 9 && MIPCHI2DV_k2 > 9 && MIPCHI2DV_k3 > 9"
cuts_Bu += cuts_s6



logger.info('DATA chain name is %s ' % (tBu.GetName()))

data_cuts_Bu = cuts_Bu

for i in prntCuts(data_cuts_Bu, "  CUTS B+  "):
    logger.info(i)


logger.info('Fill control B+  histogram (takes some time)')
with timing():
    tBu.Project(h1.GetName(), 'DTFm_b', data_cuts_Bu)

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




logger.info('Make unbinned fit for B+')
#
model_Bu.s.setMax(1.2 * len(ds_Bu))
ru, fu = model_Bu.fitTo(ds_Bu, draw=True, nbins=nbin_Bu)

model_Bu.signal.sigma.release()
model_Bu.signal.mean.release()
ru, fu = model_Bu.fitTo(ds_Bu, draw=True, nbins=nbin_Bu)

# model_Bu.signal.aR.release()
# model_Bu.signal.aL.release()
# ru, fu = model_Bu.fitTo(ds_Bu, draw=True, nbins=nbin_Bu)

# model_Bu.signal.nR.release()
# model_Bu.signal.nL.release()

# ru, fu = model_Bu.fitTo(ds_Bu, draw=True, nbins=nbin_Bu)



# print 'FIT#2 results for B+  ', ru(model_Bu.s_name)[0]
# print 'FIT#2 precision:', ru("SBu")[0].prec()


logger.info('running sPlot')
model_Bu.sPlot(ds_Bu)

h = ROOT.TH1F("h", '',  30, 5.16, 5.45)
h.Sumw2()

ds_Bu.project(h, "m_b_misid1", "SBu_sw")
h.Draw()

# d = shelve.open('$KKpidir/fit/histos.shelve')

# d['KKK_30bin_sel6'] = h


# d.close()


# plotting against MC


  
logger.info('end of the module')
