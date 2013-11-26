#!/usr/bin/env ipython
# ==========================================================================================#
# $Id:$
# ========================================================================
# @file  DiCharm/Analysis/Single.py
#
#  Various selectors for 1xCharm analysis
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

Various selectors for 1xCharm analysis

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
from AnalysisPython.PyRoUts import hID, SE, VE
from AnalysisPython.PySelector import Selector
from AnalysisPython.progress_bar import ProgressBar
# =============================================================================
from math import atan2, pi
import DiCharm.Ana as Ana
# =============================================================================


class Cuts1(object):

    def __init__(self,
                 part,
                 tos_1,
                 pt_min):

        first = part
        first_name = '_' + first
        first_trg = first + '_'
        self.tos1 = lambda s: Ana.isTos(s,  first_trg)
        self.good_1 = lambda s: getattr(s, 'good' + first_name)

        if tos_1:
            self.tos = lambda s: self.tos1
        else:
            self.tos = lambda s: True

        self.pt_1 = lambda s: getattr(s, 'pt' + first_name)
        self.y_1 = lambda s: getattr(s, 'y' + first_name)

        if 0 <= first_name.find('psi'):
            pt_min = 0

        self.pt_min = pt_min
        print ' CUTS: ptmin = %f GeV ' % self.pt_min

        if 0 <= first_name.find('psi'):
            self.ok1 = lambda s:      0 <= self.pt_1(
                s) <= 12 and 2 <= self.y_1(s) <= 4.0
        else:
            self.ok1 = lambda s: pt_min <= self.pt_1(
                s) <= 12 and 2 <= self.y_1(s) <= 4.0

    def __call__(self, item):

        ok = True

        ok = ok and self.ok1(item)
        ok = ok and self.good_1(item)
        ok = ok and self.tos(item)

        return ok


class Weights1(object):

    def __init__(self,
                 first,
                 tos_1):

        first_name = '_' + first
        first_trg = first + '_'

        if 0 <= first_name.find('psi'):
            self._eff_fun_1 = lambda s: Ana.eff_Jpsi(s,  first_name)
            self._eff_trg_1 = lambda s: Ana.trgEff_Jpsi(s,  first_name)
        elif 0 <= first_name.find('D0'):
            self._eff_fun_1 = lambda s: Ana.eff_D0(s,  first_name)
            self._eff_trg_1 = lambda s: Ana.trgEff_D0(s,  first_name)
        elif 0 <= first_name.find('Dp'):
            self._eff_fun_1 = lambda s: Ana.eff_Dp(s,  first_name)
            self._eff_trg_1 = lambda s: Ana.trgEff_Dp(s,  first_name)
        elif 0 <= first_name.find('Ds'):
            self._eff_fun_1 = lambda s: Ana.eff_Ds(s,  first_name)
            self._eff_trg_1 = lambda s: Ana.trgEff_Ds(s,  first_name)
        elif 0 <= first_name.find('Lc'):
            self._eff_fun_1 = lambda s: Ana.eff_Lc(s,  first_name)
            self._eff_trg_1 = lambda s: Ana.trgEff_Lc(s,  first_name)
        else:
            raise AttributeError("Invalid first  name '%s'" % first_name)

        if tos_1:
            self._eff_trg = lambda s: self._eff_trg_1(s)
        else:
            self._eff_trg = lambda s:           VE(1, 0)

        self.pt_1 = lambda s: getattr(s, 'pt' + first_name)
        self.y_1 = lambda s: getattr(s, 'y' + first_name)

    def __call__(self, item):

        ## acceptance & reconstruction & selection
        e1 = self._eff_fun_1(item)

        # PID: pions, kaons & protons
        ePi = Ana.pidEff_pions(item)
        eK = Ana.pidEff_kaons(item)
        eP = Ana.pidEff_protons(item)

        # correct for track recontruction efficiency
        eTr = Ana.recEff_tracks(item)

        # trigger efficiency
        eTrg = self._eff_trg(item)

        eff = VE(1, 0)

        eREC = e1
        ePID = ePi * eK * eP
        eTRG = eTrg
        eTRK = eTr

        #

        eff = VE(1, 0)

        eff *= eREC
        eff *= ePID
        eff *= eTRG
        eff *= eTRK

        # final result
        weight = 1.0

        if 0 < eff.value():
            weight = 1.0 / eff.value()
        else:
            weight = 0.0

        if 0 == weight or weight > 5.e+4:
            print ' ZERO weight : ', weight, \
                  (  self.pt_1( item ), self.y_1 ( item ) ), \
                  (e1, ePi, eK, eP, eTr, eTrg)

        return weight


# ========================================================================
# Create&Fill the basic dataset for RooFit
#  @date   2011-07-22
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
class C1(Selector):

    """
    Create and fill the basic dataset for Lambda_c
    """

    def __init__(self,
                 cuts,
                 dvar,
                 weight,
                 accept):

        Selector.__init__(self, None, self)  # initialize the base

        self._cuts = cuts
        self._weight = weight

        self._dvar = dvar[0]

        self.m = ROOT.RooRealVar(
            "m_" + dvar[0], "mass(C)", dvar[1], dvar[2])
        self.pt = ROOT.RooRealVar("pt", "pt(C)", 0, 12.0)
        self.y = ROOT.RooRealVar("y", "y(C)", 2.0,  4.5)
        self.lv01 = ROOT.RooRealVar("lv01", "lv01", -1.01,  1.01)
        self.chi2dtf = ROOT.RooRealVar(
            "chi2dtf", "chi2(dtf)/ndf",  0,  1.e+100)

        self.weight = ROOT.RooRealVar("weight", "weight", 0.0,  1.e+20)

        self.varset = ROOT.RooArgSet(
            #
            self.m,
            self.pt,
            self.y,
            self.lv01,
            self.chi2dtf,
            #
            # efficiency weight
            self.weight
        )

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
            "Charm",
            "Charm",
            #
            self.varset
        )

        self._events = 0
        self._counter = SE()
        self._accept = accept
        self._progress = None

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

        self._events += 1
        if 0 == self._events % 100000:
            print self._events

        if 0 != self._events % self._accept:
            return 0

        #
        # == for more convenience
        #
        bamboo = self.fChain

        y = getattr(bamboo, 'y_' + self._dvar)
        if not 2.0 <= y <= 4.0:
            return 0

        pt = getattr(bamboo, 'pt_' + self._dvar)
        if not 0.0 <= pt <= 12.0:
            return 0

        # apply cuts
        if not self . _cuts(bamboo):
            return 0

        # calculate & store the efficiency weight
        w = self._weight(bamboo)

        if 0 == w:
            return 0  # skip invalid weights

        self.weight    . setVal(w)

        self.pt        . setVal(pt)
        self.y         . setVal(y)

        mass = getattr(bamboo, 'm_' + self._dvar)
        lv01 = getattr(bamboo, 'lv01_' + self._dvar)
        dtf = getattr(bamboo, 'dtf_' + self._dvar)

        self.m         . setVal(mass)
        self.lv01      . setVal(lv01)
        self.chi2dtf   . setVal(dtf)

        # GEC
        self.nPV       . setVal(bamboo.nPV_rec)
        self.nSpd      . setVal(bamboo.nSpd_gec)
        self.nBest     . setVal(bamboo.nBest_rec)
        self.nLong     . setVal(bamboo.nLong_rec)
        self.nOT       . setVal(bamboo.nOT_gec)

        self.data .add(self.varset)

        return 1

# ========================================================================
# Create&Fill  the basic dataset for RooFit
#  @date   2011-07-22
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class Psi(C1):

    """
    Create and fill the basic dataset for RooFit
    """

    def __init__(self,
                 cuts,
                 weight=lambda s: 1,
                 accept=50000):

        C1.__init__(self,
                    lambda s: 3.0 <= s.m_psi <= 3.2 and cuts(s),
                    ('psi', 3.0, 3.2),
                    weight,
                    accept)

# ========================================================================
# Create&Fill  the basic dataset for RooFit
#  @date   2011-07-22
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class D0(C1):

    """
    Create and fill the basic dataset for RooFit
    """

    def __init__(self,
                 cuts,
                 weight=lambda s: 1,
                 accept=50000):

        C1.__init__(self,
                    lambda s: 1.8 <= s.m_D0 <= 1.92 and cuts(s),
                    ('D0', 1.8, 1.92),
                    weight,
                    accept)

# ========================================================================
# Create&Fill  the basic dataset for RooFit
#  @date   2011-07-22
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class Dp(C1):

    """
    Create and fill the basic dataset for RooFit
    """

    def __init__(self,
                 cuts,
                 weight=lambda s: 1,
                 accept=10000):

        C1.__init__(self,
                    lambda s: 1.82 <= s.m_Dp <= 1.91 and cuts(s),
                    ('Dp', 1.82, 1.91),
                    weight,
                    accept)

# ========================================================================
# Create&Fill  the basic dataset for RooFit
#  @date   2011-07-22
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class Ds(C1):

    """
    Create and fill the basic dataset for RooFit
    """

    def __init__(self,
                 cuts,
                 weight=lambda s: 1,
                 accept=5000):

        C1.__init__(self,
                    lambda s: 1.9 <= s.m_Ds <= 2.04 and cuts(s),
                    ('Ds', 1.9, 2.04),
                    weight,
                    accept)

# ========================================================================
# Create&Fill  the basic dataset for RooFit
#  @date   2011-07-22
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class Lc(C1):

    """
    Create and fill the basic dataset for RooFit
    """

    def __init__(self,
                 cuts,
                 weight=lambda s: 1,
                 accept=50):

        C1.__init__(self,
                    lambda s: 2.24 <= s.m_Lc <= 2.33 and cuts(s),
                    ('Lc', 2.24, 2.33),
                    weight,
                    accept)


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
