#!/usr/bin/env ipython
# =============================================================================
# $Id:$
# =============================================================================
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
    'TrgD0',  # py-selector to study trigger efficiency for D0
    'TrgDp',  # py-selector to study trigger efficiency for D+
    'TrgDs',  # py-selector to study trigger efficiency for Ds+
    'TrgLc',  # py-selector to study trigger efficiency for Lambda_c+
    'TrgPsi',  # py-selector to study trigger efficiency for J/psi
    'TrgPsiP',  # py-selector to study trigger efficiency for   psi'
)
# =============================================================================
import ROOT
from AnalysisPython.PyRoUts import hID, axis_bins
from AnalysisPython.PySelector import Selector
from AnalysisPython.progress_bar import ProgressBar

# =============================================================================
# @class TrgSelector
#  ROOT selector for getting Tis/Tos information from prompt charm particle
#  @date   2011-07-04
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#


class TrgSelector(Selector):

    """
    ROOT selector for getting Tis/Tos information from prompt charm particle  
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
                 var_,
                 histo,
                 pt_axis,
                 y_axis,
                 l0=True,
                 l1=True,
                 l2=True,
                 cuts=lambda s: True):

        Selector.__init__(self, None, self)  # initialize the base

        self._l0 = l0
        self._l1 = l1
        self._l2 = l2

        trg = self._l0 or self._l1 or self._l2
        if not trg:
            raise TypeError, 'Trigger is not specified '

        self._cuts = cuts
        self._pt_axis = pt_axis
        self. _y_axis = y_axis

        self._histos = {}
        self._events = 0L
        #
        self.pt_ = 'pt_' + var_
        self.y_ = 'y_' + var_
        self.m_ = 'm_' + var_
        #
        self.l0_tos = var_ + '_l0tos_1'
        self.l0_tis = var_ + '_l0tis_2'
        #
        self.l1_tos = var_ + '_l1tos_1'
        self.l1_tis = var_ + '_l1tis_2'
        #
        self.l2_tos = var_ + '_l2tos_1'
        self.l2_tis = var_ + '_l2tis_2'

        self._var = var_

        for iPt in pt_axis:

            for iY in y_axis:

                hTotal = histo.Clone(hID())
                hExTos = histo.Clone(hID())
                hTisTos = histo.Clone(hID())
                hExTis = histo.Clone(hID())
                hTob = histo.Clone(hID())

                hTotal  . SetTitle('Total')
                hExTis  . SetTitle('ExTis')
                hExTos  . SetTitle('ExTos')
                hTisTos . SetTitle('Tis&Tos')
                hTob    . SetTitle('Tob')

                self._histos[(iPt, iY)] = [
                    hTotal,
                    hExTos,
                    hTisTos,
                    hExTis,
                    hTob
                ]

        #
        ## global bin
        #
        hTotal = histo.Clone(hID())
        hExTos = histo.Clone(hID())
        hTisTos = histo.Clone(hID())
        hExTis = histo.Clone(hID())
        hTob = histo.Clone(hID())

        hTotal  . SetTitle('Total')
        hExTis  . SetTitle('ExTis')
        hExTos  . SetTitle('ExTos')
        hTisTos . SetTitle('Tis&Tos')
        hTob    . SetTitle('Tob')

        self._histos0 = [
            hTotal,
            hExTos,
            hTisTos,
            hExTis,
            hTob
        ]

        #
        # progress bar
        #
        self._progress = None
        self._events = 0

    def histos(self):
        print ' '
        return self._histos, self._pt_axis, self._y_axis, self._histos0

    # start processing
    def Begin(self):
        """
        Start processing
        """
        print 80 * '*'
        print 'Begin     %s    L0 %s  L1 %s  L2 %s ' % (self._var,
                                                        self._l0,
                                                        self._l1,
                                                        self._l2)
    # end processing

    def Terminate(self):
        """
        End processing
        """
        print '\nTerminate %s    L0 %s  L1 %s  L2 %s ' % (self._var,
                                                          self._l0,
                                                          self._l1,
                                                          self._l2)

        n_total = self._histos0[0] . accumulate()
        n_extos = self._histos0[1] . accumulate()
        n_tistos = self._histos0[2] . accumulate()
        n_extis = self._histos0[3] . accumulate()
        n_tob = self._histos0[4] . accumulate()

        print ' #tot %.0f #extos %.0f #tistos %.0f #extis %.f #tob %.0f ' % (n_total  . value(),
                                                                             n_extos  . value(
                                                                             ),
                                                                             n_tistos . value(
                                                                             ),
                                                                             n_extis  . value(
                                                                             ),
                                                                             n_tob    . value())
        eTOS = 1 / (1 + n_extis / n_tistos)
        eTIS = 1 / (1 + n_extos / n_tistos)

        print ' eTOS=%s%%  eTIS=%s%% ' % (eTOS * 100, eTIS * 100)

        print 80 * '*'

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

        if not self._progress:
            print ' Trg-Chain entries: %d' % self.fChain.GetEntries()
            self._progress = ProgressBar(
                0,
                self.fChain.GetEntries(),
                80,
                mode='fixed')

        self._events += 1
        if 0 == self._events % 1000:
            self._progress.increment_amount(1000)
            print self._progress, '\r',

        #
        # == for more convenience
        #
        bamboo = self.fChain

        #
        # apply trivial "acceptance" cuts
        #
        pt = getattr(bamboo, self.pt_)
        y = getattr(bamboo, self.y_)

        if not 2 <= y <= 4.5:
            return 0  # RETURN
        if not pt <= 12.0:
            return 0  # RETURN

        #
        # check cuts
        #
        if not self._cuts(bamboo):
            return 0  # RETURN

        #
        # Finally: find the histogram and  fill it !
        #

        ptBin = self._pt_axis . FindBin(pt)
        yBin = self. _y_axis . FindBin(y)

        key = (ptBin, yBin)
        if not key in self._histos:
            return 0
        histos = self._histos[key]
        histos0 = self._histos0

        #
        # extract TisTos information:
        #

        # ALL:
        tos = True
        tis = True

        #
        # L0
        if self._l0 and (tis or tos):
            #
            tos_l0 = 1 == getattr(bamboo, self.l0_tos)
            tis_l0 = 1 == getattr(bamboo, self.l0_tis)
            #
            tos = tos and tos_l0
            tis = tis and tis_l0

        #
        # HLT1
        if self._l1 and (tis or tos):
            #
            tos_l1 = 1 == getattr(bamboo, self.l1_tos)
            tis_l1 = 1 == getattr(bamboo, self.l1_tis)

            tos = tos and tos_l1
            tis = tis and tis_l1

        #
        # HLT2
        if self._l2 and (tis or tos):
            tos_l2 = 1 == getattr(bamboo, self.l2_tos)
            tis_l2 = 1 == getattr(bamboo, self.l2_tis)

            tos = tos and tos_l2
            tis = tis and tis_l2

        #
        # TIS/TOS categories
        exTos = tos and not tis
        exTis = tis and not tos
        tisTos = tis and tos

        mass = getattr(bamboo, self.m_)
        #
        # All:
        if True:
            histos[0] . Fill(mass)  # ALL
            histos0[0] . Fill(mass)  # ALL
        #
        # start TisTos'ing:
        if exTos:
            histos[1] . Fill(mass)  # ExTos
            histos0[1] . Fill(mass)  # ExTos
        elif tisTos:
            histos[2] . Fill(mass)  # Tis&Tos
            histos0[2] . Fill(mass)  # Tis&Tos
        elif exTis:
            histos[3] . Fill(mass)  # ExTis
            histos0[3] . Fill(mass)  # ExTis
        else:
            histos[4] . Fill(mass)  # Tob
            histos0[4] . Fill(mass)  # Tob
        #
        return 1


pt_bins_psi = [0, 0.5,
               1, 1.5,
               2, 2.5,
               3, 3.5,
               4, 4.5,
               5, 6, 7, 8, 10, 12]

pt_bins_had = [2, 2.5,
               3, 3.5,
               4, 4.5,
               5, 5.5,
               6, 7, 8, 10, 12]

pt_bins_lc = [2.0,
              2.5,
              3.0, 4, 5, 6, 7, 8, 10, 12]


y_bins_all = [2, 2.5, 3.0, 3.5, 4.0, 4.5]


pt_axis_psi = axis_bins(pt_bins_psi)
pt_axis_had = axis_bins(pt_bins_had)
pt_axis_lc = axis_bins(pt_bins_lc)
pt_axis_ds = axis_bins(pt_bins_lc)
y_axis_all = axis_bins(y_bins_all)

pt_axis_psip = axis_bins([3, 5, 8, 12])
y_axis_psip = axis_bins([2, 3, 4,  4.5])

pt_axis_ups = axis_bins([0, 1, 2, 3, 4, 5, 7, 10, 15])
y_axis_ups = y_axis_all

# =============================================================================
# the specific selector for D0 -> K- pi+


class TrgD0(TrgSelector):

    """
    The specific selector for D0 -> K- pi+
    """
    # constructor from cuts and histogram bins
    #  @param pt_axis the pt-axis binnig
    #  @param y_axis  the y-axis binnig
    #  @param l0      perform tistos'ing of L0
    #  @param l1      perform tistos'ing of L1
    #  @param l2      perform tistos'ing of L2
    #  @param cuts    the cuts (in any)

    def __init__(self,
                 pt_axis=pt_axis_had,
                 y_axis=y_axis_all,
                 l0=True,
                 l1=True,
                 l2=True,
                 cuts=lambda s: 1.82 <= s.m_D0 <= 1.92):

        TrgSelector.__init__(
            self,
            'D0',  # see TrgEff
            ROOT.TH1F(hID(), 'D0-template', 100, 1.82, 1.92),
            pt_axis,
            y_axis,
            l0,
            l1,
            l2,
            cuts)

# =============================================================================
# the specific selector for D+ -> K- pi+ pi+


class TrgDp (TrgSelector):

    """
    The specific selector for D+ -> K- pi+ pi+
    """
    # constructor from cuts and histogram bins
    #  @param pt_axis the pt-axis binnig
    #  @param y_axis  the y-axis binnig
    #  @param l0      perform tistos'ing of L0
    #  @param l1      perform tistos'ing of L1
    #  @param l2      perform tistos'ing of L2
    #  @param cuts    the cuts (in any)

    def __init__(self,
                 pt_axis=pt_axis_had,
                 y_axis=y_axis_all,
                 l0=True,
                 l1=True,
                 l2=True,
                 cuts=lambda s: 1.82 <= s.m_Dp <= 1.91):

        TrgSelector.__init__(
            self,
            'Dp',  # see TrgEff
            ROOT.TH1F(hID(), 'D+-template', 90, 1.82, 1.91),
            pt_axis,
            y_axis,
            l0,
            l1,
            l2,
            cuts)


# =============================================================================
# the specific selector for Lambda_c+ -> p K- pi+
class TrgLc(TrgSelector):

    """
    The specific selector for Lambda_c+ -> p K- pi+ 
    """
    # constructor from cuts and histogram bins
    #  @param pt_axis the pt-axis binnig
    #  @param y_axis  the y-axis binnig
    #  @param l0      perform tistos'ing of L0
    #  @param l1      perform tistos'ing of L1
    #  @param l2      perform tistos'ing of L2
    #  @param cuts    the cuts (in any)

    def __init__(self,
                 pt_axis=pt_axis_lc,
                 y_axis=y_axis_all,
                 l0=True,
                 l1=True,
                 l2=True,
                 cuts=lambda s: 2.24 <= s.m_Lc <= 2.33):

        TrgSelector.__init__(
            self,
            'Lc',  # see TrgEff
            ROOT.TH1F(hID(), 'c+-template', 90, 2.24, 2.33),
            pt_axis,
            y_axis,
            l0,
            l1,
            l2,
            cuts)

# =============================================================================
# the specific selector for Ds+ -> (K- K+)pi+


class TrgDs(TrgSelector):

    """
    The specific selector for Ds+ -> (K- K+)pi+
    """
    # constructor from cuts and histogram bins
    #  @param pt_axis the pt-axis binnig
    #  @param y_axis  the y-axis binnig
    #  @param l0      perform tistos'ing of L0
    #  @param l1      perform tistos'ing of L1
    #  @param l2      perform tistos'ing of L2
    #  @param cuts    the cuts (in any)

    def __init__(self,
                 pt_axis=pt_axis_ds,
                 y_axis=y_axis_all,
                 l0=True,
                 l1=True,
                 l2=True,
                 cuts=lambda s: 1.9 <= s.m_Ds <= 2.04):

        TrgSelector.__init__(
            self,
            'Ds',  # see TrgEff
            ROOT.TH1F(hID(), 'Ds+-template', 140, 1.9, 2.04),
            pt_axis,
            y_axis,
            l0,
            l1,
            l2,
            cuts)


# =============================================================================
# the specific selector for J/psi -> mu+ mu-
class TrgPsi(TrgSelector):

    """
    The specific selector for J/psi -> mu+ mu-
    """
    # constructor from cuts and histogram bins
    #  @param pt_axis the pt-axis binnig
    #  @param y_axis  the y-axis binnig
    #  @param l0      perform tistos'ing of L0
    #  @param l1      perform tistos'ing of L1
    #  @param l2      perform tistos'ing of L2
    #  @param cuts    the cuts (in any)

    def __init__(self,
                 pt_axis=pt_axis_psi,
                 y_axis=y_axis_all,
                 l0=True,
                 l1=True,
                 l2=True,
                 cuts=lambda s: 3.0 <= s.m_psi <= 3.2):

        TrgSelector.__init__(
            self,
            'psi',  # see TrgEff
            ROOT.TH1F(hID(), 'Psi+-template', 200, 3.0, 3.2),
            pt_axis,
            y_axis,
            l0,
            l1,
            l2,
            cuts)

# =============================================================================
# @class TrgUps
#  ROOT selector for getting Tis/Tos information for Y -> mu+ mu-
#  @date   2012-05-13
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class TrgUps(TrgSelector):

    """
    The specific selector for Y -> mu+ mu-
    """
    # constructor from cuts and histogram bins
    #  @param pt_axis the pt-axis binnig
    #  @param y_axis  the y-axis binnig
    #  @param l0      perform tistos'ing of L0
    #  @param l1      perform tistos'ing of L1
    #  @param l2      perform tistos'ing of L2
    #  @param cuts    the cuts (in any)

    def __init__(self,
                 pt_axis=pt_axis_ups,
                 y_axis=y_axis_ups,
                 l0=True,
                 l1=True,
                 l2=True,
                 cuts=lambda s: 8.5 <= s.m_ups <= 11):

        TrgSelector.__init__(
            self,
            'ups',  # see TrgEff
            ROOT.TH1F(hID(), "Y-template", 250, 8.5, 11),
            pt_axis,
            y_axis,
            l0,
            l1,
            l2,
            cuts)

# =============================================================================
# @class TrgPsiP
#  ROOT selector for getting Tis/Tos information for psi' -> mu+ mu-
#  @date   2012-05-13
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class TrgPsiP(TrgSelector):

    """
    The specific selector for psi' -> mu+ mu-
    """
    # constructor from cuts and histogram bins
    #  @param pt_axis the pt-axis binnig
    #  @param y_axis  the y-axis binnig
    #  @param l0      perform tistos'ing of L0
    #  @param l1      perform tistos'ing of L1
    #  @param l2      perform tistos'ing of L2
    #  @param cuts    the cuts (in any)

    def __init__(self,
                 pt_axis=pt_axis_psip,
                 y_axis=y_axis_psip,
                 l0=True,
                 l1=True,
                 l2=True,
                 cuts=lambda s: 3.61 <= s.m_psip <= 3.77):

        TrgSelector.__init__(
            self,
            'psip',  # see TrgEff
            ROOT.TH1F(hID(), "Psi'-template", 160, 3.61, 3.77),
            pt_axis,
            y_axis,
            l0,
            l1,
            l2,
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
