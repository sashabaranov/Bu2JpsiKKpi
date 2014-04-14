#!/usr/bin/env ipython
# ==========================================================================================#
# $Id:$
# ========================================================================
# @file  DiCharm/Trigger/GetTrgEff.py
#
#  Calculate trigger efficiencies  for charm particles
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
#  @date   2011-07-04
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#
#                    $Revision$
#  Last modification $Date$
#                 by $Author$
# =============================================================================
"""

Calculate trigger efficiencies  for charm particles 

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
__date__ = "2011-07-04"
__version__ = '$Revision$'
# =============================================================================
__all__ = (
    'MuPidSelector',  # py-selector to study mu-PID effiicency for J/psi
    'MuPidPsi',  # py-selector to study mu-PID effiicency for J/psi
)
# =============================================================================
import ROOT
from AnalysisPython.PyRoUts import hID
from AnalysisPython.PySelector import Selector
# =============================================================================
# @class MuPidSelector
#  ROOT selector for getting MuPid information from prompt charm particle
#  @date   2011-07-04
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#


class MuPidSelector(Selector):

    """
    ROOT selector for getting mu-PID information from prompt charm particle  
    """
    # constructor from cuts and histogram bins
    #  @param var_    the prefix for trgger variable in n-tuple
    #  @param histo   the mass histogram template
    #  @param pt_axis the pt-axis binnig
    #  @param y_axis  the y-axis binnig
    #  @param cuts    the cuts (in any)

    def __init__(self,
                 histo,
                 pt_axis=ROOT.TAxis(20, 2.0, 12.0),
                 y_axis=ROOT.TAxis(10, 2.0,  4.5),
                 cuts=lambda s: True):

        Selector.__init__(self, None, self)  # initialize the base

        self._cuts = cuts
        self._pt_axis = pt_axis
        self. _y_axis = y_axis

        self._histos = {}
        self._events = 0L
        #
        for iPt in pt_axis:

            for iY in y_axis:

                hTotal = histo.Clone(hID())
                hAccept = histo.Clone(hID())
                hReject = histo.Clone(hID())

                hTotal  . SetTitle('Total')
                hAccept . SetTitle('Accept')
                hReject . SetTitle('Reject')

                self._histos[(iPt, iY)] = [
                    hTotal,
                    hAccept,
                    hReject,
                ]

    def histos(self):
        return self._histos, self._pt_axis, self._y_axis

    # the only one important method
    def Process(self, entry):
        """
        Fills the histograms from a tree
        """
        #
        # == getting the next entry from the tree
        #
        if self.GetEntry(entry) <= 0:
            return 0  # RETURN

        self._events += 1
        if 0 == self._events % 500000:
            print self._events

        #
        # == for more convenience
        #
        bamboo = self.fChain

        #
        # apply trivial "acceptance" cuts
        #
        if not 2 <= bamboo.y <= 4.5:
            return 0  # RETURN
        if not bamboo.pt < 12.0:
            return 0  # RETURN

        #
        # check cuts
        #
        if not self._cuts(bamboo):
            return 0  # RETURN

        #
        # Finally: find the histogram and  fill it !
        #

        ptBin = self._pt_axis . FindBin(bamboo.pt)
        yBin = self. _y_axis . FindBin(bamboo. y)

        key = (ptBin, yBin)
        if not key in self._histos:
            return 0
        histos = self._histos[key]

        #
        # extract muPID
        #

        mu1 = bamboo.DLLmu1
        mu2 = bamboo.DLLmu2

        mass = bamboo.mass
        #
        # All:
        if True:
            histos[0] . Fill(mass)  # ALL
        #
        if min(mu1, mu2) > 0:
            histos[1] . Fill(mass)  # Accept
        else:
            histos[2] . Fill(mass)  # Reject
        #
        return 1

# =============================================================================

pt_axis_psi = ROOT.TAxis(24, 0.0, 12)
y_axis_psi = ROOT.TAxis(10, 2.0, 4.5)

# =============================================================================
# the specific selector for J/psi -> mu+ mu-


class MuPidPsi(MuPidSelector):

    """
    The specific selector for J/psi -> mu+ mu- 
    """
    # constructor from cuts and histogram bins
    #  @param pt_axis the pt-axis binnig
    #  @param y_axis  the y-axis binnig
    #  @param cuts    the cuts (in any)

    def __init__(self,
                 pt_axis=pt_axis_psi,
                 y_axis=y_axis_psi,
                 cuts=lambda s: 3.0 <= s.mass <= 3.2):

        MuPidSelector.__init__(self,
                               ROOT.TH1F(
                                   hID(), 'Psi+-template', 200, 3.0, 3.2),
                               pt_axis,
                               y_axis,
                               cuts)

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
