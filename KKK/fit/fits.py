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
from data_Bc import tSelection5 as tBu
from data_Bc import tBu_mc, tBu_mc_p6, tBu_mc_p8


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


from Selectors import SBT

sel_Bu = SBT(m_Bu, cuts_Bu)


logger.info(
    'Build RooFit dataset for B+ , it could take as long as 3-5 minutes')

tBu.process(sel_Bu)

ds_Bu = sel_Bu.dataset()
ds_Bu.Print('v')



with rooSilent():
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

print 'FIT#2(data 5) results for B+  ', ru(model_Bu.s_name)[0]
print 'FIT#2 precision:', ru("SBu")[0].prec()



logger.info('running sPlot')
model_Bu.sPlot(ds_Bu)

h = ROOT.TH1F("h", '',  30, 5.16, 5.45)
h.Sumw2()

ds_Bu.project(h, "m_b_misid", "SBu_sw")
h.Draw()

d = shelve.open('$KKpidir/fit/histos.shelve')

d['KKK_30bin_sel5'] = h


d.close()

# logger.info('running sPlot')
# model_Bu.sPlot(ds_Bu)

# hh_rd = ROOT.TH1F("test", '', 10, 1.9685 - 6 * 0.00032, 1.9685 + 6 * 0.00032)
# hh_rd.Sumw2()

# hh_rd.SetXTitle("M(KKK), GeV/c^2")

# ds_Bu.project(hh_rd, "DTFm_KKK", "SBu_sw")

# hh_rd.Draw()

# h = ROOT.TH1F("h", '',  30, 5.1, 5.45)
# h.Sumw2()

# ds_Bu.project(h, "m_b_misid", "SBu_sw")
# h.Draw()

# d = shelve.open('new_result.shelve')

# d['KKK_hist'] = h

# d.close()


# plotting against MC


  
logger.info('end of the module')
