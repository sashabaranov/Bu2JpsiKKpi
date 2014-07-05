import ROOT
from Ostap.PyRoUts import *
from Ostap.Utils import timing
from Ostap.Utils import rooSilent
from AnalysisPython.Logger import getLogger

from PyPAW.Selectors import SelectorWithVars

logger = getLogger(__name__)

interval = (5.1, 6.5)
binning = 150

def make_hist(ds, name, variable, cuts):
    h = ROOT.TH1F(name, '', binning, *interval)
    h.Sumw2()

    ds.project(h, variable, cuts)

    # for j in xrange(0, h.GetNbinsX()):
    #     if h.GetBinContent(j) < 0:
    #         h.SetBinContent(j, 0)

    # h.Scale(1.0/h.Integral())

    return h

def make_hist_mc(tBu, name, variable, cuts):
    h = ROOT.TH1F(name, '', binning, *interval)
    h.Sumw2()

    tBu.Project(h.GetName(), variable, cuts)

    # for j in xrange(0, h.GetNbinsX()):
    #     if h.GetBinContent(j) < 0:
    #         h.SetBinContent(j, 0)

    # h.Scale(1.0/h.Integral())

    return h

def make_legend(title, mc=True):
    "p6_k1_cuts -> misid K1 + cuts(Pythia6)"

    vals = title.split('_')
    if mc:
        if len(vals) == 2:
            return "Misid - Pythia{}".format(vals[0][1])
        if len(vals) == 3:
            return "Misid + cut (Pythia{})".format(vals[0][1])
    else:
        if len(vals) == 1:
            return "Misid (RD)"
        if len(vals) == 2:
            return "Misid + cut (RD)"
