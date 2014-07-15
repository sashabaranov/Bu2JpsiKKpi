#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
import ROOT
from ROOT import SetOwnership

from PyPAW.PyRoUts import *
import PyPAW.FitModels as Models

from AnalysisPython.Logger import getLogger

logger = getLogger(__name__)

# =============================================================================
_nemax = 80000  ## number of events per CPU-core
_ncmax =     6  ## maximal number of CPUs: there are some problems with >= 7
                ## @see https://sft.its.cern.ch/jira/browse/ROOT-4897
#
_ncpus = []

def ncpu (  events ) :
    #
    #
    ### return  ROOT.RooFit.Save()
    #
    n  = events // _nemax
    if n       <= 1 : return ROOT.RooFit.Save() ## fake!!!
    #
    import multiprocessing
    n_cores = multiprocessing.cpu_count()
    if n_cores <= 1 : return ROOT.RooFit.Save () ## fake!!!
    #
    num = min ( n , n_cores , _ncmax )
    if not _ncpus :
        _ncpus.append ( num )
    #
    return ROOT.RooFit.NumCPU ( num )


def legend_kkk(title):
    "p6_k1_cuts -> misid K1 + cuts(Pythia6)"

    vals = title.split('_')
    if 'p6' in title or 'p8' in title:
        if len(vals) == 2:
            return "KKK \, misid \, of \, \, " + vals[1].upper() + " \, (Pythia{})".format(vals[0][1])
        if len(vals) == 3:
            return "KKK \, misid \, of \,  \, " + vals[1].upper() + " \, + \, cut \, (Pythia{})".format(vals[0][1])
    else:
        if len(vals) == 1:
            return "KKK \, misid \, of \, \, " + vals[0].upper() + "\,(RD)"
        if len(vals) == 2:
            return "KKK \, misid \, of \, \, " + vals[0].upper() + " \, + \, cut \, (RD)"


def legend_kpipi(title):
    "p6_k1_cuts -> misid K1 + cuts(Pythia6)"

    vals = title.split('_')
    if 'p6' in title or 'p8' in title:
        if len(vals) == 2:
            return "K\\pi\\pi \, misid \, (Pythia{})".format(vals[0][1])
        if len(vals) == 3:
            return "K\\pi\\pi \, misid \, + \, cut \, (Pythia{})".format(vals[0][1])
    else:
        if len(vals) == 1:
            return "K\\pi\\pi \, misid \, (RD)"
        if len(vals) == 2:
            return "K\\pi\\pi \, misid \, + \, cut \, (RD)"


class Charm3_pdf (object):

    """
    """

    def __init__(self,
                 signal,
                 background=None,
                 signal2=None,
                 signal3=None,
                 suffix=''):

        self.signal = signal
        self.mass = self.signal.mass

        if background:
            self.background = background
        else:
            self.background = Bkg_pdf('Background' + suffix, self.mass)

        self.s = ROOT.RooRealVar(
            "S" + suffix, "Signal", 1000,  0,  1.e+8)
        self.b = ROOT.RooRealVar(
            "B" + suffix, "Background",   10,  0,  1.e+8)

        self.s_name = self.s.GetName()
        self.b_name = self.b.GetName()

        if not signal2:
            self.alist1 = ROOT.RooArgList(self.signal     .pdf, self.background .pdf)
            self.alist2 = ROOT.RooArgList(self.s, self.b)
        else:
            self.signal2 = signal2
            self.s2 = ROOT.RooRealVar(
                "S2" + suffix, "Signal(2)", 0.1,  0,  220)

            if not signal3:
                self.alist1 = ROOT.RooArgList(self.signal.pdf, self.signal2.pdf, self.background.pdf)
                self.alist2 = ROOT.RooArgList(self.s, self.s2, self.b)
            else:
                self.signal3 = signal3
                self.s3 = ROOT.RooRealVar("S3" + suffix, "Signal(3)", 0.1, 0, 1.e+9)

                self.alist1 = ROOT.RooArgList(
                    self.signal.pdf,
                    self.signal2.pdf,
                    self.signal3.pdf,
                    self.background.pdf
                )
                self.alist2 = ROOT.RooArgList(
                    self.s,
                    self.s2,
                    self.s3,
                    self.b
                )



        self.pdf = ROOT.RooAddPdf("model" + suffix,
                                  "model(%s)" % suffix,
                                  self.alist1,
                                  self.alist2)

        self._splots = []

    # fit
    def fitTo(self, dataset, draw=False, nbins=100, refit=True, *args):
        """
        Perform the fit
        """
        size = len(dataset)
        if 1 > size:
            logger.error('Empty data set!')

        nmax = max(1.2 * (size + 10 * math.sqrt(size)), 100)

        for s in self.alist2:
            s.setMin(0.0)
            s.setMax(nmax)

        self.b.setMax(1800)

        result = self.pdf.fitTo(dataset,
                                ROOT.RooFit.Save(),
                                ##ncpu ( len ( dataset ) ) ,
                                *args)

        if 0 != result.status():
            result = self.pdf.fitTo(dataset,
                                    ROOT.RooFit.Save(),
                                    ncpu ( len ( dataset ) ) ,
                                    *args)

        if draw:
            self.legend = ROOT.TLegend(0.55, 0.65, 0.9, 0.9)
            self.legend.SetFillColor(ROOT.kWhite)
            self.legend.SetTextSize(0.023)

            frame = self.mass.frame(nbins)
            dataset  .plotOn(frame, ROOT.RooFit.Name("data"))
            self.legend.AddEntry("data", "Data", "P")

            self.pdf .plotOn(frame,
                             ROOT.RooFit.Components(
                                 self.background.pdf.GetName()),
                             ROOT.RooFit.Name("background"),
                             ROOT.RooFit.LineStyle(ROOT.kDashed),
                             ROOT.RooFit.LineColor(ROOT.kBlue))
            self.legend.AddEntry("background", "Exponential  background", "L")


            self.pdf .plotOn(frame,
                             ROOT.RooFit.Components(
                                 self.signal.pdf.GetName()),
                             ROOT.RooFit.Name("signal"),
                             ROOT.RooFit.LineStyle(ROOT.kDashed),
                             ROOT.RooFit.LineColor(ROOT.kGreen))
            self.legend.AddEntry("signal", "Signal", "L")

            if hasattr(self, 'signal2'):
                self.pdf .plotOn(frame,
                                 ROOT.RooFit.Components(
                                     self.signal2.pdf.GetName()),
                                 ROOT.RooFit.Name("signal2"),
                                 ROOT.RooFit.LineStyle(ROOT.kDashed),
                                 ROOT.RooFit.LineColor(ROOT.kCyan))
                self.legend.AddEntry("signal2", legend_kpipi(self.signal2.pdf.GetName()[5:]), "L")

            if hasattr(self, 'signal3'):
                self.pdf .plotOn(frame,
                             ROOT.RooFit.Components(
                                 self.signal3.pdf.GetName()),
                             ROOT.RooFit.Name("signal3"),
                             ROOT.RooFit.LineStyle(ROOT.kDashed),
                             ROOT.RooFit.LineColor(ROOT.kMagenta))
                self.legend.AddEntry("signal3", legend_kkk(self.signal3.pdf.GetName()[5:]), "L")


            self.pdf.plotOn(frame, ROOT.RooFit.Name("total"), ROOT.RooFit.LineColor(ROOT.kRed))

            self.legend.AddEntry("total", "Total", "L")


            frame.Draw()
            self.legend.Draw()
            # SetOwnership(self.legend, 0)

            return result, frame

        return result, None

    #
    def sPlot(self,
              dataset,
              *args):

        splot = ROOT.RooStats.SPlot(rootID("sPlot_"),
                                    "sPlot",
                                    dataset,
                                    self.pdf,
                                    self.alist2)

        self._splots += [splot]

        return splot

    def fitHisto(self, histo, draw=False, refit=True, *args):
        """
        Fit the histogram
        """
        self.lst = ROOT.RooArgList(self.mass)
        self.imp = ROOT.RooFit.Import(histo)

        self.hset = ROOT.RooDataHist(
            rootID('hds_'),
            "Data set for histogram '%s'" % histo.GetTitle(),
            self.lst,
            self.imp
        )

        result = self.pdf.fitTo(self.hset,
                                ROOT.RooFit.Save(),
                                *args)

        frame = None

        if draw:
            frame = self.mass.frame()
            self.hset.plotOn(frame)

            self.pdf .plotOn(frame,
                             ROOT.RooFit.Components(
                                 self.background.pdf.GetName()),
                             ROOT.RooFit.Name("background"),
                             ROOT.RooFit.LineStyle(ROOT.kDashed),
                             ROOT.RooFit.LineColor(ROOT.kBlue))

            self.pdf .plotOn(frame,
                             ROOT.RooFit.Components(
                                 self.signal.pdf.GetName()),
                             ROOT.RooFit.Name("signal"),
                             ROOT.RooFit.LineStyle(ROOT.kDashed),
                             ROOT.RooFit.LineColor(ROOT.kGreen))

            if hasattr(self, 'signal2'):
                self.pdf .plotOn(frame,
                                 ROOT.RooFit.Components(
                                     self.signal2.GetName()),
                                 ROOT.RooFit.LineStyle(ROOT.kDashed),
                                 ROOT.RooFit.LineColor(ROOT.kPink))
                if hasattr(self, 'signal3'):
                    self.pdf .plotOn(frame,
                                 ROOT.RooFit.Components(
                                     self.signal3.GetName()),
                                 ROOT.RooFit.LineStyle(ROOT.kDashed),
                                 ROOT.RooFit.LineColor(ROOT.kMagenta))



            self.pdf .plotOn(frame, ROOT.RooFit.LineColor(ROOT.kRed))

            frame.SetXTitle('')
            frame.SetYTitle('')
            frame.SetZTitle('')

            frame.Draw()

        return result, self
