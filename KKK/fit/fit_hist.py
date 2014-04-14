import ROOT
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
from data_Bc import tBu2 as tBu
from data_Bc import tBu_mc


logger.info('DATA chain name is %s ' % (tBu.GetName()))

tBu.project(h2, "pt_kaon[2]")
h2.Draw()

data_cuts_Bu = cuts_Bu

for i in prntCuts(data_cuts_Bu, "  CUTS B+  "):
    logger.info(i)


logger.info('Fill control B+  histogram (takes some time)')
with timing():
    tBu.Project(h1.GetName(), 'DTFm_b', data_cuts_Bu)

with rooSilent():
    logger.info('Fit Bc+ & B+ histogram (check the model)')
    r, f = model_Bu.fitHisto(h1)
    model_Bu.signal.mean .release()
    model_Bu.signal.sigma.release()
    r, f = model_Bu.fitHisto(h1)
    r, f = model_Bu.fitHisto(h1, draw=True)
