import ROOT
from AnalysisPython.PyRoot import *
from PyPAW.PyRoUts import *
from PyPAW.Utils import timing
from PyPAW.Utils import rooSilent
from AnalysisPython.Logger import getLogger

from PyPAW.Selectors import SelectorWithVars

logger = getLogger(__name__)

def make_hist(ds, name, variable, cuts):
    h = ROOT.TH1F(name, '', 30, 5.16, 5.45)
    h.Sumw2()

    ds.project(h, variable, cuts)

    # for j in xrange(0, h.GetNbinsX()):
    #     if h.GetBinContent(j) < 0:
    #         h.SetBinContent(j, 0)

    h.Scale(1.0/h.Integral())

    return h

def make_hist_mc(tBu, name, variable, cuts):
    h = ROOT.TH1F(name, '', 30, 5.16, 5.45)
    h.Sumw2()

    tBu.Project(h.GetName(), variable, cuts)

    # for j in xrange(0, h.GetNbinsX()):
    #     if h.GetBinContent(j) < 0:
    #         h.SetBinContent(j, 0)

    h.Scale(1.0/h.Integral())

    return h

