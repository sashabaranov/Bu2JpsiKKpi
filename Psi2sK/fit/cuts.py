import ROOT
from AnalysisPython.PyRoUts import *
# ========================================================================
from AnalysisPython.Logger import getLogger
logger = getLogger(__name__)


# Trigger
cuts_trg = "((psi_l0tos & 2) == 2)"
cuts_trg += "&& ((psi_l1tos & 2) == 2)"
cuts_trg += "&& ((psi_l2tos & 2) == 2)"

# Cuts
cuts_ = " DTFchi2ndof > 0"
cuts_ += "&& DTFchi2ndof < 5"
cuts_ += "&& DTFctau > 0.25"
cuts_ += "&& vchi2_b < 20"
cuts_ += "&& pt_pion[0] > 0.3 && pt_pion[1] > 0.3"
cuts_ += "&& pt_kaon > 0.6"
cuts_ += "&& m_jpsi    > 3.020 && m_jpsi    < 3.135"
cuts_ += "&& minann_K  > 0.3 && minann_pi > 0.3"
cuts_ += "&& MIPCHI2DV_k > 12. && MIPCHI2DV_pi1 > 12. && MIPCHI2DV_pi2 > 12."

mctrue = " && mcTrueB && mcTruePsi && mcTruePi1 && mcTruePi2 && mcTrueK && mcTrueMu1 && mcTrueMu2"

cuts_Bu = cuts_ + "&& " + cuts_trg
cuts_mc = cuts_Bu + mctrue


def prntCuts(cuts, prefix=""):
    for cut in cuts.split("&&"):
        yield prefix + cut
