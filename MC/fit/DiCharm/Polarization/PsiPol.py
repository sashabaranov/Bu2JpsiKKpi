#!/usr/bin/env ipython
# =============================================================================
# $Id:$
# =============================================================================
# @file  DiCharm/Polarization/PsiPol.py
#
#  Study for J/psi polarization
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
#  @date   2011-10-13
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#
#                    $Revision$
#  Last modification $Date$
#                 by $Author$
# =============================================================================
"""

Study for J/psi polarization 

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
__date__ = "2011-10-13"
__version__ = '$Revision$'
# =============================================================================
__all__ = (
    'PolSelector'  # py-selector to study J/psi polarization
)
# =============================================================================
import ROOT
from AnalysisPython.PyRoUts import hID, axis_bins
from AnalysisPython.PySelector import Selector
from AnalysisPython.progress_bar import ProgressBar

import DiCharm.Efficiency as Eff
# =============================================================================
# @class PolSelector
#  ROOT selector for getting J/psi polarization
#  @date   2011-10-13
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#


class PolSelector(Selector):

    """
    ROOT selector for getting prompt J/psi polarization  
    """
    # constructor from cuts and histogram bins
    #  @param var_    the prefix for trgger variable in n-tuple
    #  @param histo   the mass histogram template
    #  @param pt_axis the pt-axis binnig
    #  @param y_axis  the y-axis binnig
    #  @param l0      perform tistos'ing of L0
    #  @param l1      perform tistos'ing of L1
    #  @param l2      perform tistos'ing of L2
    #  @param cuts    the cuts (in any)

    def __init__(self,
                 # binning according to Yanxi ZHANG
                 pt_axis=axis_bins([2, 3, 4, 5, 7, 10, 12]),
                 ## pt_axis   = axis_bins ( [4,5,7,10,12] ) ,
                 y_axis=ROOT.TAxis(4, 2, 4),
                 histo=ROOT.TH1D(hID(), '', 200, 3.0, 3.2),
                 lv01_axis=ROOT.TAxis(20, 0, 1)):

        Selector.__init__(self, None, self)  # initialize the base

        self._progress = None
        self.  pt_axis = pt_axis
        self.   y_axis = y_axis
        self.lv01_axis = lv01_axis

        self.histos = {}
        self._events = 0
        for iPt in self.pt_axis:
            self.histos[iPt] = {}
            for iY in self.y_axis:
                self.histos[iPt][iY] = {}
                for iLV in self.lv01_axis:
                    self.histos[iPt][iY][iLV] = histo.Clone(hID())

    # the only one important method
    def Process(self, entry):
        """
        Fills the histograms from a tree
        """
        #
        # == getting the next entry from the tree
        if self.GetEntry(entry) <= 0:
            return 0  # RETURN
        #
        if not self._progress:
            self._progress = ProgressBar(0,
                                         self.fChain.GetEntries(),
                                         77,
                                         mode='fixed')

        self._events += 1
        if 1 == self._events % 1000:
            self._progress.increment_amount(1000)
            print self._progress, '\r',

        #
        # == for more convenience
        #
        bamboo = self.fChain

        #
        # apply trivial "acceptance" cuts
        #
        pt = bamboo.pt_psi
        y = bamboo.y_psi
        m = bamboo.m_psi

        if not 2 <= y < 4:
            return 0  # RETURN
        if not pt < 12.0:
            return 0  # RETURN
        if not 3 <= m < 3.2:
            return 0  # RETURN

        #
        # Finally: find the histogram and  fill it !
        #
        lv01 = abs(bamboo.lv01_psi)
        #
        pt_bin = self.pt_axis  . FindBin(pt)
        if not pt_bin in self.  pt_axis:
            return 0

        y_bin = self. y_axis  . FindBin(y)
        if not y_bin in self.   y_axis:
            return 0

        lv_bin = self.lv01_axis. FindBin(lv01)
        if not lv_bin in self.lv01_axis:
            return 0

        h = self.histos[pt_bin][y_bin][lv_bin]
        e = Eff.selEff_Jpsi_3D(pt, y, lv01)

        if 0 < e.value():
            w = 1.0 / e.value()
            h.Fill(m, w)

        return 1

# =============================================================================


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
