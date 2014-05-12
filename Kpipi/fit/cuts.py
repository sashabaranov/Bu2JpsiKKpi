import ROOT
from AnalysisPython.PyRoUts import *
import AnalysisPython.ZipShelve as ZipShelve
# ========================================================================
from AnalysisPython.Logger import getLogger
logger = getLogger(__name__)


# Bu meson
m_Bu = ROOT.RooRealVar('m_Bu', 'mass(J/psiKpipi)', 5.18, 5.45)
nbin_Bu = 75



# Histograms
h1 = ROOT.TH1F(hID(), '', nbin_Bu, m_Bu.getMin(), m_Bu.getMax())
h1.Sumw2()

h_mc = ROOT.TH1F(hID(), '', nbin_Bu, m_Bu.getMin(), m_Bu.getMax())
h_mc.Sumw2()


mctrue = " && mcTrueB && mcTruePsi && mcTruePi1 && mcTruePi2 && mcTrueK && mcTrueMu1 && mcTrueMu2"

# Cuts
cuts_ = " DTFchi2ndof > 0"
cuts_ += "&& DTFchi2ndof < 5"
cuts_ += "&& DTFctau > 0.2"
cuts_ += "&& pt_kaon > 0.6"
cuts_ += "&& pt_pion[0] > 0.3 && pt_pion[1] > 0.3"
cuts_ += "&& m_jpsi    > 3.020 && m_jpsi    < 3.135"
cuts_ += "&& minann_K  > 0.2"
cuts_ += "&& ((psi_l0tos & 2) == 2)"
cuts_ += "&& minann_pi > 0.15 && MIPCHI2DV_k > 9 && MIPCHI2DV_pi1 > 9 && MIPCHI2DV_pi2 > 9" # selection6+
# cuts_ += "&& ((psi_l1tos & 2) == 2)"
# cuts_ += "&& ((psi_l2tos & 2) == 2)"
# # cuts_ += "&& minann_pi  > 0.4"
# cuts_ += "&& ((psi3_l0tos & 2) == 2)"


# Victors fancy cuts

# cuts_ = " DTFchi2ndof > 0"
# cuts_ += "&& DTFchi2ndof < 5"
# cuts_ += "&& DTFctau > 0.15"
# cuts_ += "&& minann_pi  > 0.4"
# cuts_ += "&& ((psi3_l0tos & 2) == 2)"
# cuts_ += "&& mcTrueB == 1"

# cuts_ += "&& pt_pion[0] > 0.2 && pt_pion[0] < 2.2"
# cuts_ += "&& pt_pion[1] > 0.2 && pt_pion[1] < 2.2"


cuts_Bu = cuts_


def prntCuts(cuts, prefix=""):
    for cut in cuts.split("&&"):
        yield prefix + cut
