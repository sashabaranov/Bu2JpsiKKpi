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

        self.mass = mass

        # cuts_ =  " DTFchi2ndof > 0"
        # cuts_ += "&& DTFchi2ndof < 4"
        # cuts_ += "&& DTFctau > 0.4"
        # cuts_ += "&& c2ip_b    < 4"
        # cuts_ += "&& m_jpsi    > 3.020 && m_jpsi    < 3.135"
        # cuts_ += "&& minann_K  > 0.2   && minann_pi > 0.2"
        # cuts_ += "&& minann_mu > 0.2"

        self.DTFm_b = ROOT.RooRealVar("DTFm_b", "", 5.0, 5.5)
        self.DTFchi2ndof = ROOT.RooRealVar("DTFchi2ndof", "", -0.5, 5.0)
        self.DTFctau = ROOT.RooRealVar("DTFctau", "", 0.0, 10.0)
        self.c2ip_b = ROOT.RooRealVar("c2ip_b", "", -0.5, 20.5)
        self.m_jpsi = ROOT.RooRealVar("m_jpsi", "", 3.0, 3.2)
        self.minann_K = ROOT.RooRealVar("minann_K", "", -0.5, 1.5)
        self.minann_pi = ROOT.RooRealVar("minann_pi", "", -0.5, 1.5)
        self.minann_mu = ROOT.RooRealVar("minann_mu", "", -0.5, 1.5)

        self.DTFm_jpsikk = ROOT.RooRealVar("DTFm_jpsikk", "", 0.0, 10.0)
        self.m_Phi = ROOT.RooRealVar("m_Phi", "", 0.0, 10.0)
        self.DTFm_kpi = ROOT.RooRealVar("DTFm_kpi", "", 0.0, 10.0)
        self.DTFm_kk = ROOT.RooRealVar("DTFm_kk", "", 0.0, 10.0)
        self.DTFm_kkpi = ROOT.RooRealVar("DTFm_kkpi", "", 0.0, 10.0)

        self.ann_pion = ROOT.RooRealVar("ann_pion", "", -0.1, 1.1)
        self.ann_kaon = ROOT.RooRealVar("ann_kaon", "", -0.1, 1.1)



        self.varset = ROOT.RooArgSet(self.mass)

        self.varset.add(self.DTFchi2ndof)
        self.varset.add(self.DTFctau)
        self.varset.add(self.c2ip_b)
        self.varset.add(self.m_jpsi)
        self.varset.add(self.minann_K)
        self.varset.add(self.minann_pi)
        self.varset.add(self.minann_mu)

        self.varset.add(self.DTFm_jpsikk)
        self.varset.add(self.m_Phi)
        self.varset.add(self.DTFm_kpi)
        self.varset.add(self.DTFm_kkpi)
        self.varset.add(self.DTFm_b)
        self.varset.add(self.DTFm_kk)

        self.varset.add(self.ann_pion)
        self.varset.add(self.ann_kaon)

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
        self.DTFm_b. setVal(bamboo.DTFm_b)
        self.DTFchi2ndof.setVal(bamboo.DTFchi2ndof)
        self.DTFctau.setVal(bamboo.DTFctau)
        self.c2ip_b.setVal(bamboo.c2ip_b)
        self.m_jpsi.setVal(bamboo.m_jpsi)
        self.minann_K.setVal(bamboo.minann_K)
        self.minann_pi.setVal(bamboo.minann_pi)
        self.minann_mu.setVal(bamboo.minann_mu)

        self.DTFm_jpsikk.setVal(bamboo.DTFm_jpsikk)
        self.m_Phi.setVal(bamboo.DTFm_kk)
        self.DTFm_kpi.setVal(bamboo.DTFm_kpi)
        self.DTFm_kkpi.setVal(bamboo.DTFm_kkpi)
        self.DTFm_kk.setVal(bamboo.DTFm_kk)
        
        self.ann_pion.setVal(float(bamboo.ann_pion[0]))
        self.ann_kaon.setVal(float(bamboo.ann_kaon[1]))
        
        self.data .add(self.varset)

        return 1


# =============================================================================
# The END
# =============================================================================
