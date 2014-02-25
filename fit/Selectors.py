#!/usr/bin/env ipython
# ==========================================================================================#
# $Id:$
# ========================================================================
# @file  Selectors.py
#
#  Various selectors
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

Various selectors 

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
from AnalysisPython.PyRoUts import cpp, hID, SE
from AnalysisPython.PySelector import Selector, SelectorWithCuts
from AnalysisPython.progress_bar import ProgressBar
from AnalysisPython.Logger import getLogger

from variables import varlist

# =============================================================================
TisTosTob = cpp.ITisTos.TisTosTob

_name_ = __name__
if not _name_ or '__main__' == _name_:
    _name_ = 'BcT/Selectors'
logger = getLogger(_name_)


# ========================================================================
# Create&Fill  the basic dataset for RooFit
#  @date   2013-05-07
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
# class SZ0(Selector) :
class SBT(SelectorWithCuts):

    """
    Create and fill the basic dataset for RooFit
    """
    def add_variable(self, name, description, low, high):
        self.variables[name] = ROOT.RooRealVar(name, description, low, high)
        self.varset.add(self.variables[name])

    def __init__(self,
                 mass,  # mass-variable
                 selection,  # Tree-selection
                 cuts=lambda s: True,
                 weight=lambda s: 1,
                 name='',
                 tmva_weights='',
                 short=False):

        if not name:
            name = hID()

        SelectorWithCuts.__init__(self, selection)  # initialize the base

        self._cuts = cuts
        self._weight = weight
        self._short = short

        self.variables = {}

        self.mass = mass
        self.varset = ROOT.RooArgSet(self.mass)


        for v in varlist:
            self.add_variable(v[0], v[4], v[2] - 0.1, v[3] + 0.1,)

        self.data = ROOT.RooDataSet(
            #
            name,
            "Bu -> J/psi KK pi: " + self.mass.GetName(),
            #
            self.varset
        )

        ROOT.SetOwnership(self.data, False)

        self.mmin = self.mass.getMin()
        self.mmax = self.mass.getMax()

        self._events = 0
        self._counter = SE()
        self._progress = None
        self._total = 1

        self._tmva_weights = tmva_weights

    def __del__(self):

        self.data.Clear()
        self.data.reset()
        del self.data

    def dataset(self):
        return self.data

    def Terminate(self):
        print ''
        logger.info('Terminate:  %d/%d %s ' %
                    (self._events, self._total, self.cuts()))
        self.data.Print('v')

    def Init(self, chain):

        self._mlp = lambda: 0.5

        if chain and self._tmva_weights:

            from tmva_reader import Cutter

            self._mlp = Cutter(chain, self._tmva_weights)

        return SelectorWithCuts.Init(self, chain)

    # the only one important method
    def Process(self, entry):
        """
        Fills data set 
        """
        #
        # == getting the next entry from the tree
        #
        if self.GetEntry(entry) < 0:
            return 0  # RETURN
        #

        if not self._progress:
            self._total = self.fChain.GetEntries()
            logger.info('TChain entries: %d' % self._total)
            self._progress = ProgressBar(
                0,
                self._total,
                80,
                mode='fixed')

        if 0 == self._events % 1000:
            self._progress.increment_amount(1000)

        self._events += 1
        #
        # == for more convenience
        #
        bamboo = self.fChain

        if not self.mmin <= bamboo.DTFm_b < self.mmax:
            return 0
        #
        # apply cuts
        if not self . _cuts(bamboo):
            return 0

        # calculate & store the efficiency weight
        w = self._weight(bamboo)
        self._counter += 0 != w

        vv = self._mlp()
        # if vv < 0.75 : return 0

        self.mass   . setVal(bamboo.DTFm_b)
        for k, v in self.variables.items():
            v.setVal(getattr(bamboo, k))

        
        self.data .add(self.varset)

        return 1


# =============================================================================
# The END
# =============================================================================
