import ROOT
from AnalysisPython.PyRoUts import *
import AnalysisPython.ZipShelve as ZipShelve
# ========================================================================
from AnalysisPython.Logger import getLogger
logger = getLogger(__name__)

m_Bu = ROOT.RooRealVar('m_Bu', 'mass(J/psi K pi pi~misid~K)', 5.25, 5.55)
nbin_Bu = 40

h1 = ROOT.TH1F(hID(), '', nbin_Bu, m_Bu.getMin(), m_Bu.getMax())
h1.Sumw2()

cuts_ = " DTFchi2ndof > 0"
cuts_ += "&& DTFchi2ndof < 5"
cuts_ += "&& DTFctau > 0.2"
cuts_ += "&& pt_kaon[0] > 0.6 && pt_pion[1] > 0.6"
cuts_ += "&& pt_pion[0] > 0.3"
cuts_ += "&& m_jpsi    > 3.020 && m_jpsi    < 3.135"
cuts_ += "&& minann_K  > 0.2"

cuts_Bu = cuts_


def prntCuts(cuts, prefix=""):
    for cut in cuts.split("&&"):
        yield prefix + cut
