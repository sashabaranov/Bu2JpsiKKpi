#!/usr/bin/env ipython
# ========================================================================
# @file DiCharm/Alg
#
#  Simple template algorithm for Charm+X
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
#  @date   2011-08-18
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#
#                    $Revision$
#  Last modification $Date$
#                 by $Author$
# =============================================================================
"""

Simple template algorithm for 2xCharm

This file is a part of BENDER project:
    ``Python-based Interactive Environment for Smart and Friendly Physics Analysis''

The project has been designed with the kind help from Pere MATO and Andrey TSAREGORODTSEV. 

And it is based on the LoKi project:
    ``C++ ToolKit for Smart and Friendly Physics Analysis''

By usage of this code one clearly states the disagreement 
with the smear campaign of Dr.O.Callot et al.: 
    ``No Vanya's lines are allowed in LHCb/Gaudi software.''

"""
# =============================================================================
__author__ = 'Vanya BELYAEV  Ivan.Belyaev@cern.ch'
__date__ = "2011-07-18"
__version__ = '$Revision$'
# =============================================================================
## needed to produce/visualize the histograms
import ROOT
from Bender.Main import *  # import all bender goodies
## easy access to various geometry routines
import LHCbMath.Types
from Gaudi.Configuration import *  # needed for job configuration
# =============================================================================
from GaudiKernel.SystemOfUnits import GeV, MeV, mm
from GaudiKernel.PhysicalConstants import c_light
# =============================================================================
# logging
# =============================================================================
from Bender.Logger import getLogger
logger = getLogger(__name__)
# =============================================================================
import BenderTools.TisTos  # add methods for TisTos
import DiCharm.Pids  # add methods for Pid/Track information
# =============================================================================
from DiCharm.Alg import DiCharmAlg
# =============================================================================

# =============================================================================
# @class CharmY
#  Simple template algorithm to study Charm & Y
#  @date   2011-05-27
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class CharmY (DiCharmAlg):

    """
    Simple algorithm to study Charm & Z production
    """

    # the only one esential method:
    def analyse(self):
        """
        The only one essential method
        """

        primaries = self.vselect('PVs', ISPRIMARY)
        if primaries.empty():
            return self.Warning('No primary vertices are found', SUCCESS)

        #
        # rec summary
        rc_summary = self.get('/Event/Rec/Summary').summaryData()
        #

        # D0&Y
        tup = self.nTuple('D0Y')
        D0_Y = self.select(
            'D0Y', '[ Meson -> ( D0 -> K- pi+) ( Meson -> mu+ mu- ) ]CC ')
        for pair in D0_Y:
            if M1(pair) < 8 * GeV:
                continue
            self.addItem(tup,
                         pair,
                         'Y', 'Upsilon', self.goodY,
                         'D0', 'D0', self.goodD0)
            #
            # add some reco-summary information
            self.addRecSummary(tup, rc_summary)
            #
            tup.write()

        # D+&Y
        tup = self.nTuple('DpY')
        Dp_Y = self.select(
            'DpY', '[ Meson -> ( D+ -> K- pi+ pi+) ( Meson -> mu+ mu- ) ]CC ')
        for pair in Dp_Y:
            if M1(pair) < 8 * GeV:
                continue
            self.addItem(tup,
                         pair,
                         'Y', 'Upsilon', self.goodY,
                         'Dp', 'D+', self.goodDp)
            #
            # add some reco-summary information
            self.addRecSummary(tup, rc_summary)
            #
            tup.write()

        # Ds+&Z
        tup = self.nTuple('DsY')
        Ds_Y = self.select(
            'DsY', '[ Meson -> ( D_s+ -> K- K+ pi+) ( Meson -> mu+ mu- ) ]CC ')
        for pair in Ds_Y:
            if M1(pair) < 8 * GeV:
                continue
            self.addItem(tup,
                         pair,
                         'Y', 'Upsilon', self.goodY,
                         'Ds', 'Ds+', self.goodDs)
            #
            # add some reco-summary information
            self.addRecSummary(tup, rc_summary)
            #
            tup.write()

        # Lc+&Y
        tup = self.nTuple('LcY')
        Lc_Y = self.select(
            'LcY', '[ Meson -> ( Lambda_c+ -> p+ K- pi+) ( Meson -> mu+ mu- ) ]CC ')
        for pair in Lc_Y:
            if M1(pair) < 8:
                continue
            self.addItem(tup,
                         pair,
                         'Y', 'Upsilon', self.goodY,
                         'Lc', 'Lc+', self.goodLc)
            #
            # add some reco-summary information
            self.addRecSummary(tup, rc_summary)
            #
            tup.write()

        return SUCCESS


# =============================================================================
# configure the job
def configure(datafiles,
              catalogs=[],
              castor=False):
    """
    Job configuration 
    """

    from Configurables import DaVinci  # needed for job configuration

    the_year = '2011'
    rootInTES = '/Event/Charm2Mu/Charm'

    davinci = DaVinci(
        DataType=the_year,
        InputType='MDST',
        Simulation=False,
        PrintFreq=1000,
        EvtMax=-1,
        #
        HistogramFile='CharmY_Histos.root',
        TupleFile='CharmY.root',
        #
        Lumi=True,
        #
    )

    from Configurables import CondDB
    CondDB(LatestGlobalTagByDataType=the_year)

    # from Configurables import Gaudi__IODataManager as IODataManager
    # IODataManager().AgeLimit = 2

    # ------- decoding set-up start ----------
    from BenderTools.MicroDST import uDstConf
    uDstConf(rootInTES)
    # ------- decoding set-up end  -----------

    # come back to Bender
    setData(datafiles, catalogs, castor)

    #
    # start Gaudi
    #
    gaudi = appMgr()

    #
    # Set properties of the TisTosTools
    #
    for t in (gaudi.tool('CharmY.L0TriggerTisTos'),
              gaudi.tool('CharmY.TriggerTisTos')):

        # t . UseParticle2LHCbIDsMap = 2
        t . PropertiesPrint = True

    algY = CharmY(
        'CharmY',  # Algorithm name ,
        RootInTES=rootInTES,
        Inputs=['Phys/DiMuonAndCharmForPromptCharm/Particles']
    )
    #

    mainSeq = gaudi.algorithm('GaudiSequencer/DaVinciUserSequence', True)
    mainSeq.Members += [algY.name()]

    return SUCCESS


# =============================================================================
# The actual job steering
if '__main__' == __name__:

    # import DiCharm.Samples as Samples
    # from   DiCharm.Utils     import lfn2disk

    #data = Samples.DATA_2011_1fb_samples

    # dicharm        = lfn2disk (
    # Samples.DATA_samples [ 'DiCharm/Down'     ] +
    #   Samples.DATA_samples [ 'DiCharm/Up'       ]  )

    # psi_and_charm  = lfn2disk (
    #    Samples.DATA_samples [ 'DiMuon&Charm/Down'] +
    #    Samples.DATA_samples [ 'DiMuon&Charm/Up'  ]
    #    )

    ## lfns = Samples.DATA_samples [ 'DiMuon&Charm/Down']
    ## + Samples.DATA_samples [ 'DiMuon&Charm/Up'  ]

    # psi_and_psi     = data [ '2xDiMuon/Up'     ]
    # psi_and_charm   = data [ 'DiMuon&Charm/Up' ]
    # charm_and_charm = data [ 'DiCharm/Up'      ]

    import shelve
    db = shelve.open('/afs/cern.ch/user/i/ibelyaev/public/LFNs/lfns.db')

    # charm_ew  = db['charm/EW-down-iso']
    # charm_ew += db['charm/EW-up-iso'  ]

    charm_2mu = db['Charm&2mu:17/Down']
    charm_2mu += db['Charm&2mu:17/Up']

    input = charm_2mu

    configure(input, castor=True)

    run(1000)

    gaudi = appMgr()

# =============================================================================
# The END
# =============================================================================
