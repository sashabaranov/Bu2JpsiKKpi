#!/usr/bin/env ipython
# ========================================================================
# @file DiCharm/BcT.py
#
#  Simple algorithm for Bc
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

Simple algorithm for Bc

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
import BenderTools.Fill  # add methods for TisTos

# =============================================================================
# @class B2PsiH
#  Simple algorithh for B -> J/psi H
#  @date   2011-05-27
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class B2PsiH (Algo):

    """
    Simple algorithm 
    """

    def initialize(self):

        sc = Algo.       initialize(self)
        if sc.isFailure():
            return sc

        triggers = {}
        triggers['psi'] = {}
        from DiCharm.TisTos import lines

        sc = self.tisTos_initialize(triggers, lines)
        if sc.isFailure():
            return sc

        sc = self.  fill_initialize()
        if sc.isFailure():
            return sc

        self.f_c2dtf = DTF_CHI2NDOF(
            True, strings(['J/psi(1S)', 'psi(2S)']))
        self.f_mfit = DTF_FUN(
            M, True, strings(['J/psi(1S)', 'psi(2S)']))
        self.f_ctau = DTF_CTAU(
            0, True, strings(['J/psi(1S)', 'psi(2S)']))

        return SUCCESS

    def finalize(self):

        self.f_c2dtf = None
        self.f_mfit = None
        self.f_ctau = None

        self.  fill_finalize()
        self.tisTos_finalize()

        return Algo.finalize(self)

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

        bees = self.select(
            'b', '[ Beauty -> ( J/psi(1S) -> mu+ mu- ) (pi+|K+) ]CC')

        tup = self.nTuple('B')
        for b in bees:

            m = M(b) / GeV
            if m > 7.1:
                continue
            #
            psi = b(1)
            h = b(2)

            m_psi = M(psi) / GeV
            if m_psi < 2.7:
                continue
            elif m_psi > 4.0:
                continue
            elif m_psi > 3.3:
                psi.setParticleID(LHCb.ParticleID(100443))

            c2 = self.f_c2dtf(b)
            if not -1 <= c2 < 10000:
                continue

            ct = self.f_ctau(b)
            if not -10 <= ct < 100:
                continue

            tup.column_float('mass', self.f_mfit(b))
            tup.column_float('c2dtf', c2)
            tup.column_float('ctau', ct)

            self.treatKine(tup, b, '_b')
            self.treatKine(tup, psi, '_psi')
            self.treatKine(tup, h, '_h')

            # add the information needed for TisTos
            self.tisTos(psi, tup, 'psi_',
                        self.lines['psi'], self.l0tistos, self.tistos)

            # add the information for Pid efficiency correction
            self.treatPions(tup, b)
            self.treatKaons(tup, b)
            self.treatMuons(tup, b)
            self.treatTracks(tup, b)

            # add some reco-summary information
            self.addRecSummary(tup, rc_summary)
            #
            tup.write()

        return SUCCESS

# =============================================================================
# configure the job


def configure(datafiles,
              catalogs=[],
              castor=False,
              params={}):
    """
    Job configuration 
    """

    ## needed for job configuration
    from Configurables import DaVinci

    from BenderTools.Parser import hasInFile

    the_year = "2011"

    if params:
        the_year = params['Year']
        logger.info('Year is set from params to be %s ' % the_year)
    else:
        if hasInFile(datafiles, 'Collision11'):
            the_year = '2011'
        elif hasInFile(datafiles, 'Collision12'):
            the_year = '2012'
        elif hasInFile(datafiles, 'Stripping17'):
            the_year = '2011'
        elif hasInFile(datafiles, 'Stripping13'):
            the_year = '2011'
        elif hasInFile(datafiles, 'Stripping15'):
            the_year = '2011'
        elif hasInFile(datafiles, 'Stripping17'):
            the_year = '2011'
        elif hasInFile(datafiles, 'Stripping19'):
            the_year = '2012'
        logger.info('Year is set from files  to be %s ' % the_year)

    #
    # check
    #
    if '2011' == the_year and hasInFile(datafiles, 'Collision12'):
        raise AttributeError, 'Invalid Year %s ' % the_year
    if '2012' == the_year and hasInFile(datafiles, 'Collision11'):
        raise AttributeError, 'Invalid Year %s ' % the_year

    logger.info('Use the Year = %s ' % the_year)

    rootInTES = '/Event/PSIX'

    # Read only fired events to speed up
    from PhysConf.Filters import LoKi_Filters
    fltrs = LoKi_Filters(
        VOID_Code="""
        ( CONTAINS('/Event/PSIX/Phys/SelPsiKForPsiX/Particles'  ) > 0.5 ) |
        ( CONTAINS('/Event/PSIX/Phys/SelPsiPiForPsiX/Particles' ) > 0.5 ) 
        """
    )

    davinci = DaVinci(
        EventPreFilters=fltrs.filters('Filters'),
        DataType=the_year,
        InputType='MDST',
        RootInTES=rootInTES,
        Simulation=False,
        PrintFreq=1000,
        EvtMax=-1,
        #
        HistogramFile='Bc_Histos.root',
        TupleFile='Bc.root',
        #
        Lumi=True,
        #
    )

    from BenderTools.Utils import silence
    silence()

    from Configurables import TrackScaleState
    state_scale = TrackScaleState('StateScale')
    davinci.UserAlgorithms = [
        state_scale,
        'B2K',
        'B2P'
    ]

    from Configurables import CondDB
    CondDB(LatestGlobalTagByDataType=the_year)

    #
    # come back to Bender
    #
    setData(datafiles, catalogs, castor)

    #
    # start Gaudi
    #
    gaudi = appMgr()

    #
    alg1 = B2PsiH(
        'B2K',  # Algorithm name ,
        RootInTES=rootInTES,
        Inputs=['Phys/SelPsiKForPsiX/Particles']
    )  # Phys/SelPsiKForPsiX/Particles

    alg2 = B2PsiH(
        'B2P',  # Algorithm name ,
        RootInTES=rootInTES,
        Inputs=['Phys/SelPsiPiForPsiX/Particles']
    )

    return SUCCESS

# =============================================================================
# The actual job steering
if '__main__' == __name__:

    input = [
        '/lhcb/LHCb/Collision12/PSIX.MDST/00024671/0000/00024671_00000001_1.psix.mdst'
    ]

    configure(input, castor=True)

    run(1000)

    gaudi = appMgr()


# =============================================================================
# The END
# =============================================================================
