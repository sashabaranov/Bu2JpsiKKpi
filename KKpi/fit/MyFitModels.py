#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
import ROOT
from ROOT import SetOwnership

from PyPAW.PyRoUts import *
import PyPAW.FitModels as Models

from AnalysisPython.Logger import getLogger

logger = getLogger(__name__)

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
                "S2" + suffix, "Signal(2)", 1000,  0,  1.e+9)

            if not signal3:
                self.alist1 = ROOT.RooArgList(self.signal.pdf, self.signal2.pdf, self.background.pdf)
                self.alist2 = ROOT.RooArgList(self.s, self.s2, self.b)
            else:
                self.signal3 = signal3
                self.s3 = ROOT.RooRealVar("S3" + suffix, "Signal(3)", 1000, 0, 1.e+9)

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

        result = self.pdf.fitTo(dataset,
                                ROOT.RooFit.Save(),
                                ## ncpu ( len ( dataset ) ) ,
                                *args)

        if 0 != result.status():
            result = self.pdf.fitTo(dataset,
                                    ROOT.RooFit.Save(),
                                    ## ncpu ( len ( dataset ) ) ,
                                    *args)

        if draw:
            self.legend = ROOT.TLegend(0.55, 0.65, 0.86, 0.9)
            self.legend.SetFillColor(ROOT.kWhite)

            frame = self.mass.frame(nbins)
            dataset  .plotOn(frame, ROOT.RooFit.Name("data"))
            self.legend.AddEntry("data", "Data", "P")

            self.pdf .plotOn(frame,
                             ROOT.RooFit.Components(
                                 self.background.pdf.GetName()),
                             ROOT.RooFit.Name("background"),
                             ROOT.RooFit.LineStyle(ROOT.kDashed),
                             ROOT.RooFit.LineColor(ROOT.kBlue))
            self.legend.AddEntry("background", "Exponential background", "L")


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
                self.legend.AddEntry("signal2", self.signal2.pdf.GetName(), "L")

            if hasattr(self, 'signal3'):
                self.pdf .plotOn(frame,
                             ROOT.RooFit.Components(
                                 self.signal3.pdf.GetName()),
                             ROOT.RooFit.Name("signal3"),
                             ROOT.RooFit.LineStyle(ROOT.kDashed),
                             ROOT.RooFit.LineColor(ROOT.kMagenta))
                self.legend.AddEntry("signal3", self.signal3.pdf.GetName(), "L")


            self.pdf.plotOn(frame, ROOT.RooFit.Name("total"), ROOT.RooFit.LineColor(ROOT.kRed))

            self.legend.AddEntry("total", "Total", "L")

            evt_binning = int(self.mass.getBinWidth(nbins) * 1000)

            frame.SetXTitle('#Inv.\,mass(J/\psi\,K K \pi), GeV/c^{2}')
            frame.SetYTitle('Events / (%d \, MeV/c^{2})' % evt_binning)
            frame.SetZTitle('')

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
