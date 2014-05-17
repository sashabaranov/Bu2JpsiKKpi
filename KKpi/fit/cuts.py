import ROOT
from AnalysisPython.PyRoUts import *
import AnalysisPython.ZipShelve as ZipShelve
# ========================================================================
from AnalysisPython.Logger import getLogger
logger = getLogger(__name__)


# Cuts
cuts_ = " DTFchi2ndof > 0"
cuts_ += "&& DTFchi2ndof < 5"
cuts_ += "&& DTFctau > 0.2"
cuts_ += "&& pt_kaon[0] > 0.6 && pt_kaon[1] > 0.6"
cuts_ += "&& pt_pion > 0.3"
cuts_ += "&& m_jpsi    > 3.020 && m_jpsi    < 3.135"
cuts_ += "&& minann_K  > 0.2"
cuts_ += "&& ((psi_l0tos & 2) == 2)"
cuts_ += "&& ann_pion > 0.15 && minann_K > 0.15 && MIPCHI2DV_k1 > 9 && MIPCHI2DV_k2 > 9 && MIPCHI2DV_pi > 9"

mctrue = "&& mcTrueB && mcTruePsi && mcTrueK1 && mcTrueK2 && mcTruePi && mcTrueMu1 && mcTrueMu2"

Nsigma = 5


reflection1 = " && (mass_PI_as_K < {} || mass_PI_as_K > {}) ".format(5.28 - Nsigma*0.005, 5.28 + Nsigma*0.005)
reflection2 = " && (mass_K_as_PI < {} || mass_K_as_PI > {}) ".format(5.28 - Nsigma*0.007, 5.28 + Nsigma*0.007)


# cuts_ += "&& ((psi_l1tos & 2) == 2)"
# cuts_ += "&& ((psi_l2tos & 2) == 2)"

# cuts_ += "&& ann_kaon[0] > 0.2 && ann_kaon[1] < 0.4"
# cuts_ += "&& ann_pion > 0.7"
# cuts_ += "&& ann_kaon[1] * (1 - ann_pion[0]) > 0.025"
# cuts_ += "&& ((psi3_l0tos & 2) == 2)"
# cuts_  += "&& ((psi3_l1tos & 2) == 2)"
# cuts_  += "&& ((psi3_l2tos & 2) == 2)"
cuts_Bu = cuts_


def prntCuts(cuts, prefix=""):
    for cut in cuts.split("&&"):
        yield prefix + cut
