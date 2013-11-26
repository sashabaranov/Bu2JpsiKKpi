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
# =============================================================================
from math import atan2, pi

# =============================================================================
# Create&Fill the histograms for Gec studies
#  @date   2011-07-22
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class PsiGec(Selector):

    """
    Create and fill the basic dataset for RooFit
    """

    def __init__(self,
                 histo,
                 pv_axis=ROOT.TAxis(10, 0,    10),
                 spd_axis=ROOT.TAxis(100, 0,  1000),
                 ot_axis=ROOT.TAxis(150, 0, 15000),
                 tlong_axis=ROOT.TAxis(100, 0,   500),
                 tbest_axis=ROOT.TAxis(100, 0,  1000)):

        Selector.__init__(self, None, self)  # initialize the base

        #
        # PV
        self.pv_axis = pv_axis
        self.pv = {}
        for i in pv_axis:
            h = histo.Clone(hID())
            h.Sumw2()
            self.pv[i] = h
        #
        # SPD
        self.spd_axis = spd_axis
        self.spd = {}
        for i in spd_axis:
            h = histo.Clone(hID())
            h.Sumw2()
            self.spd[i] = h

        # OT
        self.ot_axis = ot_axis
        self.ot = {}
        for i in ot_axis:
            h = histo.Clone(hID())
            h.Sumw2()
            self.ot[i] = h

        # Long
        self.tlong_axis = tlong_axis
        self.tlong = {}
        for i in tlong_axis:
            h = histo.Clone(hID())
            h.Sumw2()
            self.tlong[i] = h

        # Best
        self.tbest_axis = tbest_axis
        self.tbest = {}
        for i in tbest_axis:
            h = histo.Clone(hID())
            h.Sumw2()
            self.tbest[i] = h

        self._events = 0

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
        if 0 != self._events % 100:
            return 0

        if 0 != self._events % 20001:
            print self._events

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

        if not bamboo.psi_l0tos_1:
            return 0
        if not bamboo.psi_l1tos_1:
            return 0
        if not bamboo.psi_l2tos_1:
            return 0

        nPV = bamboo.nPV_rec
        nSpd = bamboo.nSpd_gec
        nOT = bamboo.nOT_gec
        nBest = bamboo.nBest_rec
        nLong = bamboo.nLong_rec

        m = bamboo.m_psi

        i = self.pv_axis .FindBin(nPV)
        if i in self.pv_axis:
            self.pv[i].Fill(m)

        i = self.spd_axis.FindBin(nSpd)
        if i in self.spd_axis:
            self.spd[i].Fill(m)

        i = self.ot_axis.FindBin(nOT)
        if i in self.ot_axis:
            self.ot[i].Fill(m)

        i = self.tlong_axis.FindBin(nLong)
        if i in self.tlong_axis:
            self.tlong[i].Fill(m)

        i = self.tbest_axis.FindBin(nBest)
        if i in self.tbest_axis:
            self.tbest[i].Fill(m)

        return 1


# =============================================================================
# Create&Fill the histograms for Gec studies
#  @date   2011-07-22
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
class DGec(PsiGec):

    """
    Create and fill the basic dataset for RooFit
    """

    def __init__(self,
                 particle,
                 histo,
                 pv_axis=ROOT.TAxis(10, 0,    10),
                 spd_axis=ROOT.TAxis(100, 0,  1000),
                 ot_axis=ROOT.TAxis(150, 0, 15000),
                 tlong_axis=ROOT.TAxis(100, 0,   500),
                 tbest_axis=ROOT.TAxis(100, 0,  1000)):

        PsiGec.__init__(self,
                        histo,
                        pv_axis,
                        spd_axis,
                        ot_axis,
                        tlong_axis,
                        tbest_axis)

        self.particle = particle
        #
        self._pt = 'pt_' + particle
        self._y = 'pt_' + particle
        self._m = 'm_' + particle
        self._good = 'good_' + particle
        self._l0tos = particle + '_l0tos_1'
        self._l1tos = particle + '_l1tos_1'
        self._l2tos = particle + '_l2tos_1'

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
        # if 0 != self._events % 10   : return 0

        if 0 != self._events % 200001:
            print self._events

        #
        # == for more convenience
        #
        bamboo = self.fChain

        pt = getattr(bamboo, self._pt)
        y = getattr(bamboo, self._y)

        #
        # acceptance
        if not 12 > pt:
            return 0
        if not 2.0 < y <= 4.0:
            return 0
        #
        # "good"
        good = getattr(bamboo, self._good)
        if not 1 == good:
            return 0  # Good J/psi

        l0tos = getattr(bamboo, self._l0tos)
        l1tos = getattr(bamboo, self._l1tos)
        l2tos = getattr(bamboo, self._l2tos)

        if not l0tos:
            return 0
        if not l1tos:
            return 0
        if not l2tos:
            return 0

        nPV = bamboo.nPV_rec
        nSpd = bamboo.nSpd_gec
        nOT = bamboo.nOT_gec
        nBest = bamboo.nBest_rec
        nLong = bamboo.nLong_rec

        m = getattr(bamboo, self._m)

        i = self.pv_axis .FindBin(nPV)
        if i in self.pv_axis:
            self.pv[i].Fill(m)

        i = self.spd_axis.FindBin(nSpd)
        if i in self.spd_axis:
            self.spd[i].Fill(m)

        i = self.ot_axis.FindBin(nOT)
        if i in self.ot_axis:
            self.ot[i].Fill(m)

        i = self.tlong_axis.FindBin(nLong)
        if i in self.tlong_axis:
            self.tlong[i].Fill(m)

        i = self.tbest_axis.FindBin(nBest)
        if i in self.tbest_axis:
            self.tbest[i].Fill(m)

        return 1


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
