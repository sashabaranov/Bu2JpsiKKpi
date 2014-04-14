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

        self.DTFchi2ndof = ROOT.RooRealVar("DTFchi2ndof", "", -0.5, 5.0)
        self.DTFctau = ROOT.RooRealVar("DTFctau", "", 0.0, 10.0)
        self.c2ip_b = ROOT.RooRealVar("c2ip_b", "", -0.5, 20.5)
        self.m_jpsi = ROOT.RooRealVar("m_jpsi", "", 3.0, 3.2)
        

        self.minann_K = ROOT.RooRealVar("minann_K", "", -0.5, 1.5)
        self.ann_kaon1 = ROOT.RooRealVar("ann_kaon1", "", -0.5, 1.5)
        self.ann_kaon2 = ROOT.RooRealVar("ann_kaon2", "", -0.5, 1.5)
        self.ann_kaon3 = ROOT.RooRealVar("ann_kaon3", "", -0.5, 1.5)

        

        self.minann_mu = ROOT.RooRealVar("minann_mu", "", -0.5, 1.5)

        # Interesting variables
        self.pt_kaon = ROOT.RooRealVar("pt_kaon", "", -0.5, 6.0)
        self.eta_kaon = ROOT.RooRealVar("eta_kaon", "", -0.5, 6.0)

        self.y_jpsi = ROOT.RooRealVar("y_jpsi", "", 0.0, 5.0)
        self.eta_jpsi = ROOT.RooRealVar("eta_jpsi", "", -0.5, 10.0)
        self.pt_jpsi = ROOT.RooRealVar("pt_jpsi", "", -0.5, 20.0)

        self.y_b = ROOT.RooRealVar("y_b", "", 0.0, 5.0)
        self.eta_b = ROOT.RooRealVar("eta_b", "", -0.5, 10.0)
        self.pt_b = ROOT.RooRealVar("pt_b", "", -0.5, 40.0)

        self.DTFm_KKK = ROOT.RooRealVar("DTFm_KKK", "", -0.5, 40.0)

        self.minPt_track = ROOT.RooRealVar("minPt_track", "", -0.5, 20.0)
        self.m_b_misid = ROOT.RooRealVar("m_b_misid", "", 0.0, 10.0)


        self.varset = ROOT.RooArgSet(self.mass)

        self.varset.add(self.DTFchi2ndof)
        self.varset.add(self.DTFctau)
        self.varset.add(self.c2ip_b)
        self.varset.add(self.m_jpsi)
        self.varset.add(self.minann_K)
        self.varset.add(self.minann_mu)

        self.varset.add(self.pt_kaon)
        self.varset.add(self.eta_kaon)
        self.varset.add(self.y_jpsi)
        self.varset.add(self.eta_jpsi)
        self.varset.add(self.pt_jpsi)
        self.varset.add(self.y_b)
        self.varset.add(self.eta_b)
        self.varset.add(self.pt_b)
        self.varset.add(self.minPt_track)
        self.varset.add(self.DTFm_KKK)
        self.varset.add(self.m_b_misid)


        self.varset.add(self.ann_kaon1)
        self.varset.add(self.ann_kaon2)
        self.varset.add(self.ann_kaon3)



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
        self.DTFchi2ndof.setVal(bamboo.DTFchi2ndof)
        self.DTFctau.setVal(bamboo.DTFctau)
        self.c2ip_b.setVal(bamboo.c2ip_b)
        self.m_jpsi.setVal(bamboo.m_jpsi)
        self.minann_K.setVal(bamboo.minann_K)
        self.minann_mu.setVal(bamboo.minann_mu)

        self.pt_kaon.setVal(float(bamboo.pt_kaon[0]))
        self.eta_kaon.setVal(float(bamboo.eta_kaon[0]))
        self.y_jpsi.setVal(bamboo.y_jpsi)
        self.eta_jpsi.setVal(bamboo.eta_jpsi)
        self.pt_jpsi.setVal(bamboo.pt_jpsi)
        
        self.y_b.setVal(bamboo.y_b)
        self.eta_b.setVal(bamboo.eta_b)
        self.pt_b.setVal(bamboo.pt_b)
        self.minPt_track.setVal(bamboo.minPt_track)

        self.DTFm_KKK.setVal(bamboo.DTFm_KKK)
        self.m_b_misid.setVal(bamboo.m_b_misid)

        self.ann_kaon1.setVal(float(bamboo.ann_kaon[0]))
        self.ann_kaon2.setVal(float(bamboo.ann_kaon[1]))
        self.ann_kaon3.setVal(float(bamboo.ann_kaon[2]))

        self.data.add(self.varset)

        return 1


# =============================================================================
# The END
# =============================================================================
