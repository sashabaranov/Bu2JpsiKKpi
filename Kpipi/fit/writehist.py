import ROOT
from AnalysisPython.PyRoUts import *
from AnalysisPython.Utils import timing
from AnalysisPython.Utils import rooSilent
from AnalysisPython.Logger import getLogger
logger = getLogger(__name__)


from cuts import m_Bu, nbin_Bu
from cuts import cuts_Bu, prntCuts
from cuts import h1, h_mc


logger.info('Import models/PDFs')
from model import model_Bu


logger.info('Import data (for the first time it could take some time)')
from data_Bc import tBu2 as tBu
from data_Bc import tBu_mc


logger.info('DATA chain name is %s ' % (tBu.GetName()))


data_cuts_Bu = cuts_Bu


f = ROOT.TFile('pi_minpt_mc_final.root','recreate')
h = ROOT.TH1F("minpt", '', 40, 0.2, 2.2)
h.SetBins(40, 0.2, 2.2)
# t = ROOT.TTree('t1','tree with histos')

# t.Branch('pi_minpt', 40, 'pi_minpt/F')
# tBu_mc.Project(h.GetName(), 'min(pt_pion[0], pt_pion[1])', data_cuts_Bu)
tBu_mc.project(h, 'min(pt_pion[0],pt_pion[1])', data_cuts_Bu)
h.SetBins(40, 0.2, 2.2)
h.Draw()

# f.Write()
# f.Close()