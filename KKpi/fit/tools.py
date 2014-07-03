import ROOT
from AnalysisPython.PyRoot import *
from PyPAW.PyRoUts import *
from PyPAW.Utils import timing
from PyPAW.Utils import rooSilent
from AnalysisPython.Logger import getLogger
from PyPAW.Selectors import SelectorWithVars

logger = getLogger(__name__)


interval = (5.2, 5.4)


def smear_kkk(h):
    h2 = ROOT.TH1F('h2', '', 50, *interval)
    h2 += h

    h3 = ROOT.TH1F('h3', '', 300, *interval)
    h3 += h2
    h3.Smooth(20)

    return h3
    # h2.Smooth(20000)
    # return h2

def smear_kpipi(h):
    h2 = ROOT.TH1F('h2', '', 50, *interval)
    h2 += h

    h3 = ROOT.TH1F('h3', '', 300, *interval)
    h3 += h2

    h3.Smooth(20)
    return h3
    # h2.Smooth(20000)
    # return h2



def make_hist_mc(tBu, name, variable, cuts):
    h = ROOT.TH1F(name, '', 30, 5.2, 5.4)
    h.Sumw2()

    tBu.Project(h.GetName(), variable, cuts)

    # for j in xrange(0, h.GetNbinsX()):
    #     if h.GetBinContent(j) < 0:
    #         h.SetBinContent(j, 0)

    h.Scale(1.0/h.Integral())

    return h

def make_legend(title):
    "p6_k1_cuts -> misid K1 + cuts(Pythia6)"

    vals = title.split('_')
    if len(vals) == 2:
        return "misid " + vals[1] + "(Pythia{})".format(vals[0][1])
    if len(vals) == 3:
        return "misid " + vals[1] + " + cut (Pythia{})".format(vals[0][1])
