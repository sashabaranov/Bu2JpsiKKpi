import ROOT
from AnalysisPython.PyRoUts import *
import AnalysisPython.ZipShelve as ZipShelve
# ========================================================================
from AnalysisPython.Logger import getLogger
logger = getLogger(__name__)


# Bu meson
m_Bu = ROOT.RooRealVar('m_Bu', 'mass(J/psiKKpi)', 5.16, 5.45)
nbin_Bu = 75


# K*(892)
m_Kstar = ROOT.RooRealVar('m_Kstar', 'mass(K+ pi-)', 0.634, 1.6)
nbin_Kstar = 30

mm_Kstar = ROOT.RooRealVar('mm_Kstar', 'mass(Kstar0)', 0.85, 0.9)
mm_Kstar.fix(0.892)

width_Kstar = ROOT.RooRealVar('width_Kstar', 'width(Kstar0)', 0.0, 1.0)

low_kpi = ROOT.RooRealVar('low', 'low_KstarBg', 0.6, 0.7)
low_kpi.fix(0.634)
high_kpi = ROOT.RooRealVar('high', 'high_KstarBg', 1.3, 2.0)



# Phi(1020)
m_Phi = ROOT.RooRealVar('m_phi', 'mass(KK)', 1.0, 1.08)
nbin_Phi = 30

mm_Phi = ROOT.RooRealVar('mm_Phi', 'mass(Phi)', 1.0, 1.1)
mm_Phi.fix(1.02)

width_Phi = ROOT.RooRealVar('width_Phi', 'width(Phi)', 0.0, 1.0)

low_kk = ROOT.RooRealVar('low', 'low_PhiBg', 0.9, 1.0)
low_kk.fix(0.986)
high_kk = ROOT.RooRealVar('high', 'high_PhiBg', 2.0, 2.1)
#high_kk.fix(2.04)



# Histograms
h1 = ROOT.TH1F(hID(), '', nbin_Bu, m_Bu.getMin(), m_Bu.getMax())
h1.Sumw2()

h2 = ROOT.TH1F(hID(), '', nbin_Kstar, m_Kstar.getMin(), m_Kstar.getMax())
h2.Sumw2()

h3 = ROOT.TH1F(hID(), '', nbin_Phi, m_Phi.getMin(), m_Phi.getMax())
h3.Sumw2()

h4 = ROOT.TH1F(hID(), '', nbin_Bu, 1.0, 3.5)
h4.Sumw2()

# Cuts
cuts_ = " DTFchi2ndof > 0"
cuts_ += "&& DTFchi2ndof < 5"
cuts_ += "&& DTFctau > 0.2"
cuts_ += "&& pt_kaon[0] > 0.6 && pt_kaon[1] > 0.6"
cuts_ += "&& pt_pion > 0.3"
cuts_ += "&& m_jpsi    > 3.020 && m_jpsi    < 3.135"
cuts_ += "&& minann_K  > 0.2"
cuts_ += "&& ((psi_l0tos & 2) == 2)"
cuts_ += "&& ((psi_l1tos & 2) == 2)"
cuts_ += "&& ((psi_l2tos & 2) == 2)"

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
