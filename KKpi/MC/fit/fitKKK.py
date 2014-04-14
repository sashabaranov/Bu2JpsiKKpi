import ROOT
from AnalysisPython.PyRoUts import *
from AnalysisPython.Utils import timing
from AnalysisPython.Utils import rooSilent
from AnalysisPython.Logger import getLogger
logger = getLogger(__name__)


from cutsKKK import m_Bu, nbin_Bu


logger.info('Import models/PDFs')
from modelKKK import model_Bu


logger.info('Import data (for the first time it could take some time)')
from cutsKKK import cuts_Bu, prntCuts


from data_KKK import tBu2 as tBu


logger.info('DATA chain name is %s ' % (tBu.GetName()))

#
# make some demo-plots
#
h1 = ROOT.TH1F(hID(), '', nbin_Bu, m_Bu.getMin(), m_Bu.getMax())
h1.Sumw2()

h2 = ROOT.TH1F(hID(), '', nbin_Bu, m_Bu.getMin(), m_Bu.getMax())
h2.Sumw2()

data_cuts_Bu = cuts_Bu

for i in prntCuts(data_cuts_Bu, "  CUTS B+  "):
    logger.info(i)

logger.info('Fill control B+  histogram (takes some time)')
with timing():
    tBu.Project(h1.GetName(), 'm_b_misid')
    tBu.Project(h2.GetName(), 'm_b_misid', data_cuts_Bu)

h1.red()
h2.blue()

h1.scale(1.0)
h2.scale(1.0)

h2.Draw()
h1.Draw("same")


# with rooSilent () :
#     logger.info ('Fit Bc+ & B+ histogram (check the model)')
#     r,f = model_Bu.fitHisto ( h1 )
#     model_Bu.signal.mean .release()
#     model_Bu.signal.sigma.release()
#     r,f = model_Bu.fitHisto ( h1 )
#     r,f = model_Bu.fitHisto ( h1 , draw = True )

# from SelectorsKKK import SBT

# sel_Bu = SBT(m_Bu, cuts_Bu)


# logger.info('Build RooFit dataset for B+ , it could take as long as 3-5 minutes')

# with timing():
#     ## logger.warning('almost skip B+')
#     ## with rooSilent () : tBu.process ( sel_Bu , 100000 )
#     with rooSilent():
#         tBu.process(sel_Bu)

# ds_Bu = sel_Bu.dataset()
# ds_Bu.Print('v')

# with rooSilent():
#     logger.info('Make unbinned fit for B+')
#     #
#     model_Bu.s.setMax(1.2 * len(ds_Bu))
#     ru, fu = model_Bu.fitTo(ds_Bu, draw=True, nbins=nbin_Bu)
#     print ru

# print ru
# # print 'FIT results for B+  ', ru ( model_Bu.s_name )[0]
# # print 'FIT precision:', ru("S")[0].prec()

# # ru,fu = model_Bu.fitTo    ( ds_Bu , draw = True , nbins = nbin_Bu )

# # print 'FIT#2 results for B+  ', ru ( model_Bu.s_name )[0]
# # print 'FIT#2 precision:', ru("S")[0].prec()

# logger.info('end of the module')
# # raw_input()
