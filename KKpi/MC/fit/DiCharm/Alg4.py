#!/usr/bin/env ipython
# ========================================================================
# @file DiCharm/Alg4.py
#
#  Simple template algorithm for Y+(2mu)
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

Simple template algorithm for Y+(2mu)

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
# @class Ydimu
#  Simple template algorithm to study Charm & Y
#  @date   2011-05-27
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class Ydimu (DiCharmAlg):

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
            self.Warning('No primary vertices are found', SUCCESS)

        #
        # rec summary
        rc_summary = self.get('/Event/Rec/Summary').summaryData()
        #

        y2mu = self.select(
            'e2mu', 'Meson -> ( Meson -> mu+ mu- ) ( Meson -> mu+ mu- ) ')

        #  2xY
        tup = self.nTuple('YY')
        for pair in y2mu:
            #
            self.addItem(tup,
                         pair,
                         'Y1', 'Upsilon', self.goodY,
                         'Y2', 'Upsilon', self.goodY)

            # add some reco-summary information
            self.addRecSummary(tup, rc_summary)
            #
            tup.write()

        #  Y&J/psi
        tup = self.nTuple('Ypsi')
        for pair in y2mu:
            #
            d1 = pair.child(1)
            d2 = pair.child(2)

            y = d1
            psi = d2

            m1 = M(d1)
            m2 = M(d2)

            if m1 < 7 * GeV and m2 > 7 * GeV:
                y = d2
                psi = d1
            elif m1 > 7 * GeV and m2 > 7 * GeV:
                continue
            elif m1 < 7 * GeV and m2 < 7 * GeV:
                continue

            self.addItem2(tup,
                          pair,
                          y, 'Y', 'Upsilon', self.goodY,
                          psi, 'psi', 'psi', self.goodJpsi)

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

    ## the_year  = '2011'
    the_year = '2012'

    rootInTES = '/Event/Y2Mu'
    davinci = DaVinci(
        DataType=the_year,
        InputType='MDST',
        Simulation=False,
        PrintFreq=1000,
        EvtMax=-1,
        #
        HistogramFile='Ydimu_Histos.root',
        TupleFile='Ydimu.root',
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
    for t in (gaudi.tool('Ydimu.L0TriggerTisTos'),
              gaudi.tool('Ydimu.TriggerTisTos')):

        # t . UseParticle2LHCbIDsMap = 2
        t . PropertiesPrint = True

    alg = Ydimu(
        'Ydimu',  # Algorithm name ,
        RootInTES=rootInTES,
        Inputs=['Phys/Ymu2/Particles']
    )
    #
    mainSeq = gaudi.algorithm('GaudiSequencer/DaVinciUserSequence', True)
    mainSeq.Members += [alg.name()]

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

    # charm_2mu  = db [ 'Charm&2mu:17/Down' ]
    # charm_2mu += db [ 'Charm&2mu:17/Up'   ]

    y_2mu = db['DiMuon.Y2mu.Stripping17;Up']
    y_2mu += db['DiMuon.Y2mu.Stripping17;Down']

    input = y_2mu

    db.close()

    configure(input, castor=True)

    run(1000)

    gaudi = appMgr()

# =============================================================================
# The END
# =============================================================================
