import ROOT
from AnalysisPython.PyRoot import *
from PyPAW.PyRoUts import *
from PyPAW.Utils import timing
from PyPAW.Utils import rooSilent
from PyPAW.LHCbStyle import *
from AnalysisPython.Logger import getLogger

from PyPAW.Selectors import SelectorWithVars

logger = getLogger(__name__)

def make_hist(ds, name, variable, cuts):
    h = ROOT.TH1F(name, '', 100, 5.0, 5.5)
    h.Sumw2()

    ds.project(h, variable, cuts)

    # for j in xrange(0, h.GetNbinsX()):
    #     if h.GetBinContent(j) < 0:
    #         h.SetBinContent(j, 0)

    # h.Scale(1.0/h.Integral())

    return h

def make_hist_mc(tBu, name, variable, cuts):
    h = ROOT.TH1F(name, '', 100, 5.0, 5.5)
    h.Sumw2()

    tBu.Project(h.GetName(), variable, cuts)

    # for j in xrange(0, h.GetNbinsX()):
    #     if h.GetBinContent(j) < 0:
    #         h.SetBinContent(j, 0)

    # h.Scale(1.0/h.Integral())

    return h

def make_legend(title, mc):
    "p6_k1_cuts -> misid K1 + cuts(Pythia6)"

    vals = title.split('_')
    if mc:
        if len(vals) == 2:
            return "misid " + vals[1] + "(Pythia{})".format(vals[0][1])
        if len(vals) == 3:
            return "misid " + vals[1] + " + cut (Pythia{})".format(vals[0][1])
    else:
        if len(vals) == 1:
            return "misid " + vals[0] + "(RD)"
        if len(vals) == 2:
            return "misid " + vals[0] + " + cut (RD)"