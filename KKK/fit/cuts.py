import ROOT
from AnalysisPython.PyRoUts import *
import AnalysisPython.ZipShelve as ZipShelve
# ========================================================================
from AnalysisPython.Logger import getLogger
logger = getLogger(__name__)


# Bu meson
m_Bu = ROOT.RooRealVar('m_Bu', 'mass(J/psiKKK)', 5.16, 5.45)
nbin_Bu = 75

# Histograms
h1 = ROOT.TH1F(hID(), '', nbin_Bu, m_Bu.getMin(), m_Bu.getMax())
h1.Sumw2()


h2 = ROOT.TH1F(hID(), '', nbin_Bu, 0.0, 3.0)
h2.Sumw2()


# Cuts
cuts_ = " DTFchi2ndof > 0"
cuts_ += "&& DTFchi2ndof < 5"
cuts_ += "&& DTFctau > 0.2"
cuts_ += "&& pt_kaon[0] > 0.6 && pt_kaon[1] > 0.6 && pt_kaon[2] > 0.6 "
cuts_ += "&& m_jpsi    > 3.020 && m_jpsi    < 3.135"
cuts_ += "&& minann_K  > 0.2"
cuts_ += "&& ((psi_l0tos & 2) == 2)"
# cuts_ += "&& ((psi_l1tos & 2) == 2)"
# cuts_ += "&& ((psi_l2tos & 2) == 2)"

# cuts_ += "&& ((psi3_l0tos & 2) == 2)"


# Weak cuts
# cuts_ = " DTFchi2ndof > 0"
# cuts_ += "&& DTFchi2ndof < 5"
# cuts_ += "&& DTFctau > 0.2"
# cuts_ += "&& pt_kaon[0] > 0.6 && pt_kaon[1] > 0.6" # && pt_kaon[2] > 0.6 "
# cuts_ += "&& m_jpsi    > 3.020 && m_jpsi    < 3.135"
# cuts_ += "&& minann_K  > 0.2"
# cuts_ += " && ((psi3_l0tos & 2) == 2)"


cuts_Bu = cuts_


def prntCuts(cuts, prefix=""):
    for cut in cuts.split("&&"):
        yield prefix + cut
