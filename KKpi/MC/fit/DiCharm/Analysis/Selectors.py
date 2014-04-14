#!/usr/bin/env ipython
# ==========================================================================================#
# $Id:$
# ========================================================================
# @file  DiCharm/Analysis/Selectors.py
#
#  Various selectors for 2xCharm analysis
#
#  This file is a part of
#  <a href="http://cern.ch/lhcb-comp/Analysis/Bender/index.html">Bender project</a>
#  <b>``Python-based Interactive Environment for Smart and Friendly
#   Physics Analysis''</b>
#
#  The package has been designed with the kind help from
#  Pere MATO and Andrey TSAREGORODTSEV.
#  And it is based on the
#  <a href="http://cern.ch/lhcb-comp/Analysis/LoKi/index.html">LoKi project:</a>
#  ``C++ ToolKit for Smart and Friendly Physics Analysis''
#
#  By usage of this code one clearly states the disagreement
#  with the smear campaign of Dr.O.Callot et al.:
#  ``No Vanya's lines are allowed in LHCb/Gaudi software.''
#
#  @date   2011-07-22
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#
#                    $Revision$
#  Last modification $Date$
#                 by $Author$
# =============================================================================
"""

Various selectors for 2xCharm analysis

This file is a part of BENDER project:
    ``Python-based Interactive Environment for Smart and Friendly Physics Analysis''

The project has been designed with the kind help from Pere MATO and Andrey TSAREGORODTSEV. 

And it is based on the LoKi project:
    ``C++ ToolKit for Smart and Friendly Physics Analysis''

By usage of this code one clearly states the disagreement 
with the smear campain of Dr.O.Callot et al.: 
    ``No Vanya's lines are allowed in LHCb/Gaudi software.''

"""
# =============================================================================
__author__ = 'Vanya BELYAEV  Ivan.Belyaev@cern.ch'
__date__ = "2011-07-22"
__version__ = '$Revision$'
# =============================================================================
import ROOT
from AnalysisPython.PyRoUts import hID, SE
from AnalysisPython.PySelector import Selector
from AnalysisPython.progress_bar import ProgressBar
from DiCharm.Analysis.Slices import Slices
# =============================================================================
from math import atan2, pi
import sys
# =============================================================================


def phi(s, suffix):
    return getattr(s, 'phi_' + suffix)


def dphi(s, suffix1, suffix2):

    phi1 = phi(s, suffix1)
    phi2 = phi(s, suffix2)

    d = float(phi1 - phi2)
    d /= pi

    while d > 1:
        d -= 2
    while d < -1:
        d += 2

    return abs(d)

# ========================================================================
# Create&Fill the basic dataset for RooFit
#  @date   2011-07-22
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class SpsiD(Selector):

    """
    Create and fill the basic dataset for RooFit
    """

    def __init__(self,
                 cuts,
                 dvar,
                 weight,
                 jvar="psi",
                 h1=None,
                 h2=None,
                 h3=None,
                 h4=None):

        Selector.__init__(self, None, self)  # initialize the base

        self._jvar = jvar

        self._cuts = cuts
        self._weight = weight

        self.m_2c = ROOT.RooRealVar(
            "m_2c", "mass(J/psiD)",  3.0,  100)
        self.m_psi = ROOT.RooRealVar(
            "m_psi", "mass(J/psi)",  3.0,  3.2)
        self.pt_psi = ROOT.RooRealVar("pt_psi", "pt(J/psi)",  0, 12.0)
        self.y_psi = ROOT.RooRealVar("y_psi", "y(J/psi)",  2.0,  4.5)
        self.lv01_psi = ROOT.RooRealVar(
            "lv01_psi", "lv01(J/psi)", -1.01,  1.01)
        self.chi2dtf = ROOT.RooRealVar(
            "chi2dtf", "chi2(dtf)/ndf",  0,  1.e+100)

        self._dvar = dvar[0]
        self.m_D = ROOT.RooRealVar(
            "m_" + dvar[0], "mass(D)", dvar[1], dvar[2])
        self.pt_D = ROOT.RooRealVar("pt_" + dvar[0], "pt(D)", 0, 12.0)
        self.y_D = ROOT.RooRealVar("y_" + dvar[0], "y(D)", 2.0,  4.5)

        self.weight = ROOT.RooRealVar("weight", "weight", 0.0,  1.e+20)
        self.dphi = ROOT.RooRealVar(
            "dphi", "|delta(phi)|/pi", 0.0,  1.0)

        self.varset = ROOT.RooArgSet(
            #
            self.m_2c,
            self.m_psi,
            self.pt_psi,
            self.y_psi,
            #
            self.m_D,
            self.pt_D,
            self.y_D,
            # efficiency weight
            self.weight
        )

        self.varset.add(self.lv01_psi)
        self.varset.add(self.dphi)
        self.varset.add(self.chi2dtf)

        self.nPV = ROOT.RooRealVar("nPV", 'n(PV)', 0,    20)
        self.nSpd = ROOT.RooRealVar("nSpd", 'n(Spd)', 0, 20000)
        self.nBest = ROOT.RooRealVar("nBest", 'n(Best)', 0,  5000)
        self.nLong = ROOT.RooRealVar("nLong", 'n(Long)', 0,  5000)
        self.nOT = ROOT.RooRealVar("nOT", 'n(OT)', 0, 50000)

        self.varset.add(self.nPV)
        self.varset.add(self.nSpd)
        self.varset.add(self.nBest)
        self.varset.add(self.nLong)
        self.varset.add(self.nOT)

        self.data = ROOT.RooDataSet(
            #
            "SpsiD",
            "Jpsi + Charm",
            #
            self.varset
        )
        #

        #
        # prepare the structure for slice fits:
        #   fill histos only in case both histo templates are supplied
        self.slices_raw = None
        self.slices_cor = None
        if h1 and h2:
            self.slices_raw = Slices(h1, h2)
        if h3 and h4:
            self.slices_cor = Slices(h3, h4)
        #
        self._events = 0
        self._counter = SE()
        #
        self._progress = None
        #
    #

    def dataset(self):
        return self.data

    # the only one important method
    def Process(self, entry):
        """
        Fills data set 
        """
        #
        # == getting the next entry from the tree
        #
        if self.GetEntry(entry) <= 0:
            return 0  # RETURN

        # if 0 == self._events % 100 :
        if not self._progress:
            self._progress = ProgressBar(0,
                                         self.fChain.GetEntries(),
                                         77,
                                         mode='fixed')

        self._progress.increment_amount()
        print self._progress, '\r',

        self._events += 1

        #
        # == for more convenience
        #
        bamboo = self.fChain

        #
        # J/psi  acceptance
        if not 12 > bamboo.pt_psi:
            return 0
        if not 2.0 < bamboo. y_psi <= 4.5:
            return 0
        #
        # "good" J/psi
        if not 1 == bamboo.good_psi:
            return 0  # Good J/psi

        #
        # apply cuts
        if not self . _cuts(bamboo):
            return 0

        # calculate & store the efficiency weight
        w = self._weight(bamboo)
        self._counter += 0 != w
        if 0 == w:
            return 0  # skip invalid weights

        self.weight    . setVal(w)

        self.m_2c      . setVal(bamboo.m_2c)

        m1 = getattr(bamboo, 'm_' + self._jvar)
        m2 = getattr(bamboo, 'm_' + self._dvar)

        self.m_psi     . setVal(m1)
        self.pt_psi    . setVal(getattr(bamboo, 'pt_' + self._jvar))
        self.y_psi     . setVal(getattr(bamboo, 'y_' + self._jvar))
        self.lv01_psi  . setVal(getattr(bamboo, 'lv01_' + self._jvar))

        self.chi2dtf   . setVal(bamboo.dtf_2c)

        self.m_D       . setVal(m2)
        self.pt_D      . setVal(getattr(bamboo, 'pt_' + self._dvar))
        self.y_D       . setVal(getattr(bamboo, 'y_' + self._dvar))

        self.dphi      . setVal(dphi(bamboo, self._jvar, self._dvar))

        # GEC
        self.nPV       . setVal(bamboo.nPV_rec)
        self.nSpd      . setVal(bamboo.nSpd_gec)
        self.nBest     . setVal(bamboo.nBest_rec)
        self.nLong     . setVal(bamboo.nLong_rec)
        self.nOT       . setVal(bamboo.nOT_gec)

        self.data .add(self.varset)

        if self.slices_raw:
            self.slices_raw.fill(m1, m2)
        if self.slices_cor:
            self.slices_cor.fill(m1, m2, w)

        return 1


# ========================================================================
# Create&Fill  the basic dataset for RooFit
#  @date   2011-07-22
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
class DD(Selector):

    """
    Create and fill the basic dataset for RooFit
    """

    def __init__(self,
                 cuts,
                 d1var,
                 d2var,
                 weight,
                 h1=None,
                 h2=None,
                 h3=None,
                 h4=None):

        Selector.__init__(self, None, self)  # initialize the base

        self._cuts = cuts
        self._weight = weight

        self.m_2c = ROOT.RooRealVar("m_2c", "mass(DD)", 3.0,  100)
        self.chi2dtf = ROOT.RooRealVar(
            "chi2dtf", "chi2(dtf)/ndf",  0,  1.e+100)

        self._d1 = d1var[0]
        self.m_D1 = ROOT.RooRealVar(
            "m_" + d1var[0], "mass(%s)" % d1var[0], d1var[1], d1var[2])
        self.pt_D1 = ROOT.RooRealVar(
            "pt_" + d1var[0], "pt(%s)" % d1var[0], 2.0, 12.0)
        self.y_D1 = ROOT.RooRealVar(
            "y_" + d1var[0], "y(%s)" % d1var[0], 2.0,  4.5)

        self._d2 = d2var[0]
        self.m_D2 = ROOT.RooRealVar(
            "m_" + d2var[0], "mass(%s)" % d2var[0], d2var[1], d2var[2])
        self.pt_D2 = ROOT.RooRealVar(
            "pt_" + d2var[0], "pt(%s)" % d2var[0], 2.0, 12.0)
        self.y_D2 = ROOT.RooRealVar(
            "y_" + d2var[0], "y(%s)" % d2var[0], 2.0,  4.5)

        self.weight = ROOT.RooRealVar(
            "weight", "efficiency-weight", 0.0,  1.e+20)
        self.dphi = ROOT.RooRealVar("dphi", "|delta(phi)|/pi", 0.0,  1.0)

        self.varset = ROOT.RooArgSet(
            #
            self.m_2c,
            self.m_D1,
            self.pt_D1,
            self.y_D1,
            #
            self.m_D2,
            self.pt_D2,
            self.y_D2,
            #
            # efficiency weight
            self.weight,
        )

        self.varset.add(self.dphi)
        self.varset.add(self.chi2dtf)

        self.nPV = ROOT.RooRealVar("nPV", 'n(PV)', 0,    20)
        self.nSpd = ROOT.RooRealVar("nSpd", 'n(Spd)', 0, 20000)
        self.nBest = ROOT.RooRealVar("nBest", 'n(Best)', 0,  5000)
        self.nLong = ROOT.RooRealVar("nLong", 'n(Long)', 0,  5000)
        self.nOT = ROOT.RooRealVar("nOT", 'n(OT)', 0, 50000)

        self.varset.add(self.nPV)
        self.varset.add(self.nSpd)
        self.varset.add(self.nBest)
        self.varset.add(self.nLong)
        self.varset.add(self.nOT)

        self.data1 = ROOT.RooDataSet(
            "DDbar",
            "D & anti-Charm",
            self.varset)

        self.data2 = ROOT.RooDataSet(
            "DD",
            "D & Charm",
            self.varset)

        #
        # prepare the structure for slice fits:
        #   fill histos only in case both histo templates are supplied
        self.slices_raw = None
        self.slices_cor = None
        if h1 and h2:
            self.slices_raw = (Slices(h1, h2), Slices(h1, h2))
        if h3 and h4:
            self.slices_cor = (Slices(h3, h4), Slices(h3, h4))

        self._events = 0
        self._counter = SE()
        self._progress = None

    #
    def dataset(self):
        return self.data1, self.data2

    # the only one important method
    def Process(self, entry):
        """
        Fills data set 
        """
        #
        # == getting the next entry from the tree
        #
        if self.GetEntry(entry) <= 0:
            return 0  # RETURN

        self._events += 1
        # if 0 == self._events % 1000 :
        #    print self._events

        if not self._progress:
            self._progress = ProgressBar(0,
                                         self.fChain.GetEntries(),
                                         77,
                                         mode='fixed')

        self._progress.increment_amount()
        print self._progress, '\r',

        #
        # == for more convenience
        #
        bamboo = self.fChain

        #
        # apply cuts
        if not self . _cuts(bamboo):
            return 0

        # calculate & store the efficiency weight
        w = self._weight(bamboo)
        self._counter += 0 != w
        if 0 == w:
            return 0  # skip invalid weights
        self.weight    . setVal(w)

        self.m_2c      . setVal(bamboo.m_2c)

        m1 = getattr(bamboo, 'm_' + self._d1)

        self.m_D1      . setVal(m1)
        self.pt_D1     . setVal(getattr(bamboo, 'pt_' + self._d1))
        self.y_D1      . setVal(getattr(bamboo, 'y_' + self._d1))

        m2 = getattr(bamboo, 'm_' + self._d2)

        self.m_D2      . setVal(m2)
        self.pt_D2     . setVal(getattr(bamboo, 'pt_' + self._d2))
        self.y_D2      . setVal(getattr(bamboo, 'y_' + self._d2))

        self.chi2dtf   . setVal(bamboo.dtf_2c)

        self.dphi      . setVal(dphi(bamboo, self._d1, self._d2))

        pid1 = getattr(bamboo, 'pid_' + self._d1)
        pid2 = getattr(bamboo, 'pid_' + self._d2)

        ccbar = pid1 * pid2 < 0

        data = self.data1 if ccbar else self.data2

        # GEC
        self.nPV       . setVal(bamboo.nPV_rec)
        self.nSpd      . setVal(bamboo.nSpd_gec)
        self.nBest     . setVal(bamboo.nBest_rec)
        self.nLong     . setVal(bamboo.nLong_rec)
        self.nOT       . setVal(bamboo.nOT_gec)

        data .add(self.varset)

        if self.slices_raw:
            if ccbar:
                self.slices_raw[0].fill(m1, m2)
            else:
                self.slices_raw[1].fill(m1, m2)

        if self.slices_cor:
            if ccbar:
                self.slices_cor[0].fill(m1, m2, w)
            else:
                self.slices_cor[1].fill(m1, m2, w)

        return 1


# ========================================================================
# Create&Fill  the basic dataset for RooFit
#  @date   2011-07-22
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
class SpsiD0(SpsiD):

    """
    Create and fill the basic dataset for RooFit
    """

    def __init__(self,
                 cuts,
                 weight=lambda s: 1,
                 slices_raw=False,
                 slices_cor=False):

        h1 = None
        h2 = None
        h3 = None
        h4 = None

        if slices_raw:
            h1 = ROOT.TH1D(hID(), '', 20, 3.00, 3.20)
            h1.Sumw2()
            h2 = ROOT.TH1D(hID(), '', 20, 1.82, 1.92)
            h2.Sumw2()
        if slices_cor:
            h3 = ROOT.TH1D(hID(), '', 10, 3.00, 3.20)
            h3.Sumw2()
            h4 = ROOT.TH1D(hID(), '', 10, 1.82, 1.92)
            h4.Sumw2()

        SpsiD.__init__(self,
                       lambda s: 3.0 <= s.m_psi <= 3.2 and 1.82 <= s.m_D0 <= 1.92 and cuts(
                           s),
                       ('D0', 1.82, 1.92),
                       weight,
                       'psi',
                       h1, h2, h3, h4)

# ========================================================================
# Create&Fill  the basic dataset for RooFit
#  @date   2011-07-22
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class SpsiDp(SpsiD):

    """
    Create and fill the basic dataset for RooFit
    """

    def __init__(self,
                 cuts,
                 weight=lambda s: 1,
                 slices_raw=False,
                 slices_cor=False):
        h1 = None
        h2 = None
        h3 = None
        h4 = None

        if slices_raw:
            h1 = ROOT.TH1D(hID(), '', 20, 3.00, 3.20)
            h1.Sumw2()
            h2 = ROOT.TH1D(hID(), '', 18, 1.82, 1.91)
            h2.Sumw2()

        if slices_cor:
            h3 = ROOT.TH1D(hID(), '', 10, 3.00, 3.20)
            h3.Sumw2()
            h4 = ROOT.TH1D(hID(), '',  9, 1.82, 1.91)
            h4.Sumw2()

        SpsiD.__init__(self,
                       lambda s: 3.0 <= s.m_psi <= 3.2 and 1.82 <= s.m_Dp <= 1.91 and cuts(
                           s),
                       ('Dp', 1.82, 1.91),
                       weight,
                       'psi',
                       h1, h2, h3, h4)

# ========================================================================
# Create&Fill  the basic dataset for RooFit
#  @date   2011-07-22
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class SpsiDs(SpsiD):

    """
    Create and fill the basic dataset for RooFit
    """

    def __init__(self,
                 cuts,
                 weight=lambda s: 1,
                 slices_raw=False,
                 slices_cor=False):

        h1 = None
        h2 = None
        h3 = None
        h4 = None

        if slices_raw:
            h1 = ROOT.TH1D(hID(), '', 10, 3.00, 3.20)
            h1.Sumw2()
            h2 = ROOT.TH1D(hID(), '', 14, 1.90, 2.04)
            h2.Sumw2()

        if slices_cor:
            h3 = ROOT.TH1D(hID(), '', 10, 3.00, 3.20)
            h3.Sumw2()
            h4 = ROOT.TH1D(hID(), '',  7, 1.90, 2.04)
            h4.Sumw2()

        SpsiD.__init__(self,
                       lambda s: 3.0 <= s.m_psi <= 3.2 and 1.9 <= s.m_Ds <= 2.04 and cuts(
                           s),
                       ('Ds', 1.9, 2.04),
                       weight,
                       'psi',
                       h1, h2, h3, h4)

# ========================================================================
# Create&Fill  the basic dataset for RooFit
#  @date   2011-07-22
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class SpsiLc(SpsiD):

    """
    Create and fill the basic dataset for RooFit
    """

    def __init__(self,
                 cuts,
                 weight=lambda s: 1,
                 slices_raw=False,
                 slices_cor=False):

        h1 = None
        h2 = None
        h3 = None
        h4 = None

        if slices_raw:
            h1 = ROOT.TH1D(hID(), '', 10, 3.00, 3.20)
            h1.Sumw2()
            h2 = ROOT.TH1D(hID(), '',  9, 2.24, 2.33)
            h2.Sumw2()

        if slices_cor:
            h3 = ROOT.TH1D(hID(), '', 10, 3.00, 3.20)
            h3.Sumw2()
            h4 = ROOT.TH1D(hID(), '',  9, 2.24, 2.33)
            h4.Sumw2()

        SpsiD.__init__(self,
                       lambda s: 3.0 <= s.m_psi <= 3.2 and 2.24 <= s.m_Lc <= 2.33 and cuts(
                           s),
                       ('Lc', 2.24, 2.33),
                       weight,
                       'psi',
                       h1, h2, h3, h4)

# ========================================================================
# Create&Fill  the basic dataset for RooFit
#  @date   2011-07-22
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class SD0D0(DD):

    """
    Create and fill the basic dataset for RooFit
    """

    def __init__(self,
                 cuts,
                 weight=lambda s: 1,
                 slices_raw=False,
                 slices_cor=False):

        h1 = None
        h2 = None
        h3 = None
        h4 = None

        if slices_raw:
            h1 = ROOT.TH1D(hID(), '', 20, 1.82, 1.92)
            h1.Sumw2()
            h2 = ROOT.TH1D(hID(), '', 20, 1.82, 1.92)
            h2.Sumw2()

        if slices_cor:
            h3 = ROOT.TH1D(hID(), '', 10, 1.82, 1.92)
            h3.Sumw2()
            h4 = ROOT.TH1D(hID(), '', 10, 1.82, 1.92)
            h4.Sumw2()

        DD.__init__(
            self,
            lambda s: 1.8 <= s.m_D01 <= 1.92 and 1.8 <= s.m_D02 <= 1.92 and cuts(
                s),
            ('D01', 1.8,  1.92),
            ('D02', 1.8,  1.92),
            weight,
            h1, h2, h3, h4)


# ========================================================================
# Create&Fill  the basic dataset for RooFit
#  @date   2011-07-22
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
class SD0Dp(DD):

    """
    Create and fill the basic dataset for RooFit
    """

    def __init__(self,
                 cuts,
                 weight=lambda s: 1,
                 slices_raw=False,
                 slices_cor=False):

        h1 = None
        h2 = None
        h3 = None
        h4 = None

        if slices_raw:
            h1 = ROOT.TH1D(hID(), '', 20, 1.82, 1.92)
            h1.Sumw2()
            h2 = ROOT.TH1D(hID(), '', 18, 1.82, 1.91)
            h2.Sumw2()
        if slices_cor:
            h3 = ROOT.TH1D(hID(), '', 10, 1.82, 1.92)
            h3.Sumw2()
            h4 = ROOT.TH1D(hID(), '',  9, 1.82, 1.91)
            h4.Sumw2()

        DD.__init__(
            self,
            lambda s: 1.8 <= s.m_D0 <= 1.92 and 1.82 <= s.m_Dp <= 1.91 and cuts(
                s),
            ('D0', 1.8,  1.92),
            ('Dp', 1.82,  1.91),
            weight,
            h1, h2, h3, h4)


# ========================================================================
# Create&Fill  the basic dataset for RooFit
#  @date   2011-07-22
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
class SD0Ds(DD):

    """
    Create and fill the basic dataset for RooFit
    """

    def __init__(self,
                 cuts,
                 weight=lambda s: 1,
                 slices_raw=False,
                 slices_cor=False):

        h1 = None
        h2 = None
        h3 = None
        h4 = None

        if slices_raw:
            h1 = ROOT.TH1D(hID(), '', 10, 1.82, 1.92)
            h1.Sumw2()
            h2 = ROOT.TH1D(hID(), '', 14, 1.90, 2.04)
            h2.Sumw2()

        if slices_cor:
            h3 = ROOT.TH1D(hID(), '', 10, 1.82, 1.92)
            h3.Sumw2()
            h4 = ROOT.TH1D(hID(), '',  7, 1.90, 2.04)
            h4.Sumw2()

        DD.__init__(
            self,
            lambda s: 1.8 <= s.m_D0 <= 1.92 and 1.9 <= s.m_Ds <= 2.04 and cuts(
                s),
            ('D0', 1.8,  1.92),
            ('Ds', 1.9,  2.04),
            weight,
            h1, h2, h3, h4)


# ========================================================================
# Create&Fill  the basic dataset for RooFit
#  @date   2011-07-22
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
class SD0Lc(DD):

    """
    Create and fill the basic dataset for RooFit
    """

    def __init__(self,
                 cuts,
                 weight=lambda s: 1,
                 slices_raw=False,
                 slices_cor=False):

        h1 = None
        h2 = None
        h3 = None
        h4 = None

        if slices_raw:
            h1 = ROOT.TH1D(hID(), '', 10, 1.82, 1.92)
            h1.Sumw2()
            h2 = ROOT.TH1D(hID(), '',  9, 2.24, 2.33)
            h2.Sumw2()

        if slices_cor:
            h3 = ROOT.TH1D(hID(), '',  8, 1.82, 1.92)
            h3.Sumw2()
            h4 = ROOT.TH1D(hID(), '',  9, 2.24, 2.33)
            h4.Sumw2()

        DD.__init__(
            self,
            lambda s: 1.8 <= s.m_D0 <= 1.92 and 2.24 <= s.m_Lc <= 2.33 and cuts(
                s),
            ('D0', 1.8,  1.92),
            ('Lc', 2.24,  2.33),
            weight,
            h1, h2, h3, h4)


# ========================================================================
# Create&Fill  the basic dataset for RooFit
#  @date   2011-07-22
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
class SDpDp(DD):

    """
    Create and fill the basic dataset for RooFit
    """

    def __init__(self,
                 cuts,
                 weight=lambda s: 1,
                 slices_raw=False,
                 slices_cor=False):

        h1 = None
        h2 = None
        h3 = None
        h4 = None

        if slices_raw:
            h1 = ROOT.TH1D(hID(), '', 18, 1.82, 1.91)
            h1.Sumw2()
            h2 = ROOT.TH1D(hID(), '', 18, 1.82, 1.91)
            h2.Sumw2()

        if slices_cor:
            h3 = ROOT.TH1D(hID(), '',  6, 1.82, 1.91)
            h3.Sumw2()
            h4 = ROOT.TH1D(hID(), '',  6, 1.82, 1.91)
            h4.Sumw2()

        DD.__init__(
            self,
            lambda s: 1.82 <= s.m_Dp1 <= 1.91 and 1.82 <= s.m_Dp2 <= 1.91 and cuts(
                s),
            ('Dp1', 1.82,  1.91),
            ('Dp2', 1.82,  1.91),
            weight,
            h1, h2, h3, h4)

# ========================================================================
# Create&Fill  the basic dataset for RooFit
#  @date   2011-07-22
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class SDpDs(DD):

    """
    Create and fill the basic dataset for RooFit
    """

    def __init__(self,
                 cuts,
                 weight=lambda s: 1,
                 slices_raw=False,
                 slices_cor=False):

        h1 = None
        h2 = None
        h3 = None
        h4 = None

        if slices_raw:
            h1 = ROOT.TH1D(hID(), '', 18, 1.82, 1.91)
            h1.Sumw2()
            h2 = ROOT.TH1D(hID(), '', 14, 1.90, 2.04)
            h2.Sumw2()

        if slices_cor:
            h3 = ROOT.TH1D(hID(), '',  6, 1.82, 1.91)
            h3.Sumw2()
            h4 = ROOT.TH1D(hID(), '',  7, 1.90, 2.04)
            h4.Sumw2()

        DD.__init__(
            self,
            lambda s: 1.82 <= s.m_Dp <= 1.91 and 1.9 <= s.m_Ds <= 2.04 and cuts(
                s),
            ('Dp', 1.82,  1.91),
            ('Ds', 1.9,  2.04),
            weight,
            h1, h2, h3, h4)


# ========================================================================
# Create&Fill  the basic dataset for RooFit
#  @date   2011-07-22
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
class SDpLc(DD):

    """
    Create and fill the basic dataset for RooFit
    """

    def __init__(self,
                 cuts,
                 weight=lambda s: 1,
                 slices_raw=False,
                 slices_cor=False):

        h1 = None
        h2 = None
        h3 = None
        h4 = None

        if slices_raw:
            h1 = ROOT.TH1D(hID(), '',  9, 1.82, 1.91)
            h1.Sumw2()
            h2 = ROOT.TH1D(hID(), '',  9, 2.24, 2.33)
            h2.Sumw2()

        if slices_cor:
            h3 = ROOT.TH1D(hID(), '',  6, 1.82, 1.91)
            h3.Sumw2()
            h4 = ROOT.TH1D(hID(), '',  6, 2.24, 2.33)
            h4.Sumw2()

        DD.__init__(
            self,
            lambda s: 1.82 <= s.m_Dp <= 1.91 and 2.24 <= s.m_Lc <= 2.33 and cuts(
                s),
            ('Dp', 1.82,  1.91),
            ('Lc', 2.24,  2.33),
            weight,
            h1, h2, h3, h4)


# ========================================================================
# Create&Fill  the basic dataset for RooFit
#  @date   2011-07-22
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
class SDsDs(DD):

    """
    Create and fill the basic dataset for RooFit
    """

    def __init__(self,
                 cuts,
                 weight=lambda s: 1,
                 slices_raw=False,
                 slices_cor=False):

        h1 = None
        h2 = None
        h3 = None
        h4 = None

        if slices_raw:
            h1 = ROOT.TH1D(hID(), '',  7, 1.90, 2.04)
            h1.Sumw2()
            h2 = ROOT.TH1D(hID(), '',  7, 1.90, 2.04)
            h2.Sumw2()

        if slices_cor:
            h3 = ROOT.TH1D(hID(), '',  7, 1.90, 2.04)
            h3.Sumw2()
            h4 = ROOT.TH1D(hID(), '',  7, 1.90, 2.04)
            h4.Sumw2()

        DD.__init__(
            self,
            lambda s: 1.9 <= s.m_Ds1 <= 2.04 and 1.9 <= s.m_Ds2 <= 2.04 and cuts(
                s),
            ('Ds1', 1.9,  2.04),
            ('Ds2', 1.9,  2.04),
            weight,
            h1, h2, h3, h4)


# ========================================================================
# Create&Fill  the basic dataset for RooFit
#  @date   2011-07-22
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
class SDsLc(DD):

    """
    Create and fill the basic dataset for RooFit
    """

    def __init__(self,
                 cuts,
                 weight=lambda s: 1,
                 slices_raw=False,
                 slices_cor=False):

        h1 = None
        h2 = None
        h3 = None
        h4 = None

        if slices_raw:
            h1 = ROOT.TH1D(hID(), '',  7, 1.90, 2.04)
            h1.Sumw2()
            h2 = ROOT.TH1D(hID(), '',  6, 2.24, 2.33)
            h2.Sumw2()

        if slices_cor:
            h3 = ROOT.TH1D(hID(), '',  7, 1.90, 2.04)
            h3.Sumw2()
            h4 = ROOT.TH1D(hID(), '',  6, 2.24, 2.33)
            h4.Sumw2()

        DD.__init__(
            self,
            lambda s: 1.9 <= s.m_Ds <= 2.04 and 2.24 <= s.m_Lc <= 2.33 and cuts(
                s),
            ('Ds', 1.9,  2.04),
            ('Lc', 2.24,  2.33),
            weight,
            h1, h2, h3, h4)


# ========================================================================
# Create&Fill  the basic dataset for RooFit
#  @date   2011-07-22
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
class SLcLc(DD):

    """
    Create and fill the basic dataset for RooFit
    """

    def __init__(self,
                 cuts,
                 weight=lambda s: 1,
                 slices_raw=False,
                 slices_cor=False):

        h1 = None
        h2 = None
        h3 = None
        h4 = None

        if slices_raw:
            h1 = ROOT.TH1D(hID(), '',  9, 2.24, 2.33)
            h1.Sumw2()
            h2 = ROOT.TH1D(hID(), '',  9, 2.24, 2.33)
            h2.Sumw2()

        if slices_cor:
            h3 = ROOT.TH1D(hID(), '',  9, 2.24, 2.33)
            h3.Sumw2()
            h4 = ROOT.TH1D(hID(), '',  9, 2.24, 2.33)
            h4.Sumw2()

        DD.__init__(
            self,
            lambda s: 2.24 <= s.m_Lc1 <= 2.33 and 2.24 <= s.m_Lc2 <= 2.33 and cuts(
                s),
            ('Lc1', 2.24,  2.33),
            ('Lc2', 2.24,  2.33),
            weight,
            h1, h2, h3, h4)


# =============================================================================
if '__main__' == __name__:

    print 80 * '*'
    print __doc__
    print ' Author  : ', __author__
    print ' Version : ', __version__
    print ' Date    : ', __date__
    print ' Symbols : ', __all__
    print 80 * '*'

# =============================================================================
# The END
# =============================================================================
