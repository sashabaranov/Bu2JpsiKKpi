import ROOT
from AnalysisPython.PyRoot import *
from AnalysisPython.PyRoUts import *

from data_Kpipi import tBu2 as tBu
from cutsKpipi import cuts_Bu, prntCuts

from AnalysisPython.Logger import getLogger
logger = getLogger(__name__)


logger.info('DATA chain name is %s ' % (tBu.GetName()))


m_Bu = ROOT.RooRealVar('m_Bu', 'mass(J/psi K K K~misid~pi)', 5.16, 5.45)
nbin_Bu = 30

h1 = ROOT.TH1F(hID(), '', nbin_Bu, m_Bu.getMin(), m_Bu.getMax())
h1.Sumw2()

tBu.Project(h1.GetName(), 'm_b_misid', cuts_Bu)

d = shelve.open('result_lowbins.shelve')

d['Kpipi_hist'] = h1

d.close()
