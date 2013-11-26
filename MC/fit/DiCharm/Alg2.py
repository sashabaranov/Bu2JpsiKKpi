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
# @class CharmW
#  Simple template algorithm to study Charm & W
#  @date   2011-05-27
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class CharmW (DiCharmAlg):

    """
    Simple algorithm to study Charm & W production
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

        # D0&W
        tup = self.nTuple('D0W')
        D0_W = self.select('D0W', '[ H_10 -> ( D0 -> K- pi+) [mu+]cc ]CC ')
        for pair in D0_W:
            self.addItem(tup,
                         pair,
                         'D0', 'D0', self.goodD0,
                         'W', 'W+', self.goodW)
            #
            # add some reco-summary information
            self.addRecSummary(tup, rc_summary)
            #
            tup.write()

        # D+&W
        tup = self.nTuple('DpW')
        Dp_W = self.select(
            'DpW', '[ H_10 -> ( D+ -> K- pi+ pi+) [mu+]cc ]CC ')
        for pair in Dp_W:
            self.addItem(tup,
                         pair,
                         'Dp', 'D+', self.goodDp,
                         'W', 'W+', self.goodW)
            #
            # add some reco-summary information
            self.addRecSummary(tup, rc_summary)
            #
            tup.write()

        # Ds+&W
        tup = self.nTuple('DsW')
        Ds_W = self.select(
            'DsW', '[ H_10 -> ( D_s+ -> K- K+ pi+) [mu+]cc ]CC ')
        for pair in Ds_W:
            self.addItem(tup,
                         pair,
                         'Ds', 'Ds+', self.goodDs,
                         'W', 'W+', self.goodW)
            #
            # add some reco-summary information
            self.addRecSummary(tup, rc_summary)
            #
            tup.write()

        # Lc+&W
        tup = self.nTuple('LcW')
        Lc_W = self.select(
            'LcW', '[ H_10 -> ( Lambda_c+ -> p+ K- pi+) [mu+]cc ]CC ')
        for pair in Lc_W:
            self.addItem(tup,
                         pair,
                         'Lc', 'Lc+', self.goodLc,
                         'W', 'W+', self.goodW)
            #
            # add some reco-summary information
            self.addRecSummary(tup, rc_summary)
            #
            tup.write()

        return SUCCESS

# =============================================================================
# @class CharmZ
#  Simple template algorithm to study Charm & Z
#  @date   2011-05-27
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class CharmZ (CharmW):

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

        # D0&Z
        tup = self.nTuple('D0Z')
        D0_Z = self.select(
            'D0Z', '[ H_10 -> ( D0 -> K- pi+) ( Z0 -> mu+ mu- ) ]CC ')
        for pair in D0_Z:
            self.addItem(tup,
                         pair,
                         'D0', 'D0', self.goodD0,
                         'Z', 'Z0', self.goodZ)
            #
#            print pair.daughters()[0]
#            import pdb; pdb.set_trace()
            D0 = pair.daughters()[0]
            tup.column('hasRich1', CHILDCUT(HASRICH,  1)(D0))
            tup.column('hasRich2', CHILDCUT(HASRICH,  2)(D0))
            tup.column('KL1 ', CHILDFUN(TINFO(101, 1e+09), 1)(D0))
            tup.column('KL2 ', CHILDFUN(TINFO(101, 1e+09), 2)(D0))
            tup.column('BPVLTIME', BPVLTIME(
                'LoKi::LifetimeFitter/lifetime:PUBLIC', 9)(D0))
            tup.column('MINTREE_MIPCHI2DV', MINTREE(
                MIPCHI2DV(''), (ISBASIC & HASTRACK))(D0))
            tup.column(
                'piPID',  MINTREE((PPINFO(602, 0, -1000) - PPINFO(603, 0, -1000)), (ABSID == 211))(D0))
            tup.column(
                'KPID',  MINTREE((PPINFO(603, 0, -1000) - PPINFO(602, 0, -1000)), (ABSID == 321))(D0))

            # add some reco-summary information
            self.addRecSummary(tup, rc_summary)
            #
            tup.write()

        # D+&Z
        tup = self.nTuple('DpZ')
        Dp_Z = self.select(
            'DpZ', '[ H_10 -> ( D+ -> K- pi+ pi+) ( Z0 -> mu+ mu- ) ]CC ')
        for pair in Dp_Z:
            self.addItem(tup,
                         pair,
                         'Dp', 'D+', self.goodDp,
                         'Z', 'Z0', self.goodZ)
            #
            # add some reco-summary information
            self.addRecSummary(tup, rc_summary)
            #
            tup.write()

        # Ds+&Z
        tup = self.nTuple('DsZ')
        Ds_Z = self.select(
            'DsZ', '[ H_10 -> ( D_s+ -> K- K+ pi+) ( Z0 -> mu+ mu- ) ]CC ')
        for pair in Ds_Z:
            self.addItem(tup,
                         pair,
                         'Ds', 'Ds+', self.goodDs,
                         'Z', 'Z0', self.goodZ)
            #
            # add some reco-summary information
            self.addRecSummary(tup, rc_summary)
            #
            tup.write()

        # Lc+&W
        tup = self.nTuple('LcZ')
        Lc_Z = self.select(
            'LcZ', '[ H_10 -> ( Lambda_c+ -> p+ K- pi+) ( Z0 -> mu+ mu- ) ]CC ')
        for pair in Lc_Z:
            self.addItem(tup,
                         pair,
                         'Lc', 'Lc+', self.goodLc,
                         'Z', 'Z0', self.goodZ)
            #
            # add some reco-summary information
            self.addRecSummary(tup, rc_summary)
            #
            tup.write()

        return SUCCESS

    def addItem2(self,
                 tup,
                 pair,
                 d1, name1, trg1, fun1,
                 d2, name2, trg2, fun2):
        CharmW.addItem2(self,
                        tup,
                        pair,
                        d1, name1, trg1, fun1,
                        d2, name2, trg2, fun2)

#        print d2
#        import pdb; pdb.set_trace()
        mup = d2.child(1)
        mum = d2.child(2)
        self.tisTos(mup, tup, 'muplus_', self.lines[trg2])
        self.tisTos(mum, tup, 'muminus_', self.lines[trg2])


# =============================================================================
# @class DiMuZ
#  Simple template algorithm to study 2mu+Z
#  @date   2011-05-27
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
class DiMuZ (DiCharmAlg):

    """
    Simple algorithm to study 2mu+W production
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

        Z2mu = self.select(
            'z2mu', ' H_10 ->  ( J/psi(1S) -> mu+ mu-) ( Z0 -> mu+ mu- ) ')

        # J/psi&Z
        tup = self.nTuple('Zpsi')
        for pair in Z2mu:
            if 3.3 * GeV < M1(pair):
                continue
            self.addItem(tup,
                         pair,
                         'psi', 'psi', self.goodJpsi,
                         'Z', 'Z0', self.goodZ)
            #
            # add some reco-summary information
            self.addRecSummary(tup, rc_summary)
            #
            tup.write()

        # Y+Z
        tup = self.nTuple('ZY')
        for pair in Z2mu:
            if not 7.5 * GeV < M1(pair) < 15 * GeV:
                continue
            self.addItem(tup,
                         pair,
                         'Y', 'Upsilon', self.goodY,
                         'Z', 'Z0', self.goodZ)
            #
            # add some reco-summary information
            self.addRecSummary(tup, rc_summary)
            #
            tup.write()

        # Z+Z
        tup = self.nTuple('ZZ')
        for pair in Z2mu:
            if M1(pair) < 15 * GeV:
                continue
            self.addItem(tup,
                         pair,
                         'Z1', 'Z0', self.goodZ,
                         'Z2', 'Z0', self.goodZ)
            #
            # add some reco-summary information
            self.addRecSummary(tup, rc_summary)
            #
            tup.write()

        return SUCCESS

# =============================================================================
# @class DiMuW
#  Simple template algorithm to study 2mu+W
#  @date   2011-05-27
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class DiMuW (DiCharmAlg):

    """
    Simple algorithm to study 2mu+W production
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

        W2mu = self.select(
            'w2mu', ' [H_10 -> ( J/psi(1S) -> mu+ mu- )  mu+]CC')

        # J/psi&W
        tup = self.nTuple('Wpsi')
        for pair in W2mu:
            if 3.3 * GeV < M1(pair):
                continue
            self.addItem(tup,
                         pair,
                         'psi', 'psi', self.goodJpsi,
                         'W', 'W+', self.goodW)
            #
            # add some reco-summary information
            self.addRecSummary(tup, rc_summary)
            #
            tup.write()

        # Y+W
        tup = self.nTuple('WY')
        for pair in W2mu:
            if not 7.5 * GeV < M1(pair) < 15 * GeV:
                continue
            self.addItem(tup,

                         pair,
                         'Y', 'Upsilon', self.goodY,
                         'W', 'W+', self.goodW)
            #
            # add some reco-summary information
            self.addRecSummary(tup, rc_summary)
            #
            tup.write()

        # Z+W
        tup = self.nTuple('WZ')
        for pair in W2mu:
            if M1(pair) < 15 * GeV:
                continue
            self.addItem(tup,
                         pair,
                         'Z', 'Z0', self.goodZ,
                         'W', 'W+', self.goodW)
            #
            # add some reco-summary information
            self.addRecSummary(tup, rc_summary)
            #
            tup.write()

        return SUCCESS

# =============================================================================

# =============================================================================
# configure the job


def configure_EW(rootInTES,
                 datafiles,
                 catalogs=[],
                 castor=False):
    """
    Job configuration 
    """

    from Configurables import DaVinci  # needed for job configuration
    from Configurables import EventSelector  # needed for job configuration
    from Configurables import FileCatalog  # needed for job configuration
    from Configurables import NTupleSvc  # needed for job configuration

    from PhysConf.Filters import LoKi_Filters

    ## the_year = '2011'
    the_year = '2012'

    davinci = DaVinci(
        DataType=the_year,
        InputType='MDST',
        Simulation=False,
        PrintFreq=1000,
        ## EventPreFilters = filters ,
        EvtMax=-1,
        #
        HistogramFile='CharmEW_Histos.root',
        TupleFile='CharmEW.root',
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

    return SUCCESS


# =============================================================================
# configure the job
def configure_CW(datafiles, catalogs=[], castor=False):
    """
    Job configuration 
    """
    #
    # common configuration
    #
    rootInTES = '/Event/CharmEW'
    configure_EW(
        rootInTES,
        datafiles,
        catalogs,
        castor)

    #
    # start Gaudi
    #
    gaudi = appMgr()

    #
    # Set properties of the TisTosTools
    #
    for t in (gaudi.tool('CharmW.L0TriggerTisTos'),
              gaudi.tool('CharmW.TriggerTisTos'),
              gaudi.tool('CharmZ.L0TriggerTisTos'),
              gaudi.tool('CharmZ.TriggerTisTos')):

        # t . UseParticle2LHCbIDsMap = 2
        t . PropertiesPrint = True

    algW = CharmW(
        'CharmW',  # Algorithm name ,
        RootInTES=rootInTES,
        Inputs=['Phys/CharmAndW/Particles']
    )
    algZ = CharmZ(
        'CharmZ',  # Algorithm name ,
        RootInTES=rootInTES,
        Inputs=['Phys/CharmAndZ/Particles']
    )
    #
    algWdimu = DiMuW(
        'DiMuW',  # Algorithm name ,
        RootInTES=rootInTES,
        Inputs=['Phys/DiMuAndW/Particles']
    )
    algZdimu = DiMuZ(
        'DiMuZ',  # Algorithm name ,
        RootInTES=rootInTES,
        Inputs=['Phys/DiMuAndZ/Particles']
    )

    mainSeq = gaudi.algorithm('GaudiSequencer/DaVinciUserSequence', True)
    mainSeq.Members += [algW     . name(),
                        algZ     . name(),
                        algWdimu . name(),
                        algZdimu . name()]

    return SUCCESS


# =============================================================================
# configure the job
def configure(datafiles, catalogs=[], castor=False):
    """
    Configure the job
    """
    return configure_CW(datafiles, catalogs, castor)

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

    from socket import gethostname as hostname
    runningAtUZH = 'physik.uzh.ch' in hostname()
    if not runningAtUZH:
        import shelve
        db = shelve.open('/afs/cern.ch/user/i/ibelyaev/public/LFNs/lfns.db')

        charm_ew = db['charm/EW-down-iso']
        charm_ew += db['charm/EW-up-iso']
    else:
        from glob import glob
        charm_ew = glob(
            '/home/hep/bursche/wd/data/R12S17mDSTS/v/*/*/EW.CharmEW.mdst')

    configure(charm_ew, [], True)

    run(-1)

    gaudi = appMgr()

# =============================================================================
# The END
# =============================================================================
