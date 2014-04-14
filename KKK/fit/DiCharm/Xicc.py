#!/usr/bin/env ipython
# ========================================================================
# @file DiCharm/Xicc.py
#
#  Simple template algorithm for Xicc
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

Simple template algorithm for Xicc

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

# =============================================================================
# @class Xicc
#  Simple template algorithm to search for Xicc
#  @date   2011-05-27
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class Xicc(Algo):

    """
    Simple algorithm 
    """
    # standard constructor

    def __init__(self, name, **args):
        """
        Standard constructor
        """
        Algo.__init__(self, name, **args)
        #
        # modes
        #
        self.mode = " [ Xi_cc+ -> Lambda_c+ (K-|K+) (pi+|pi-) ]CC"

    def initialize(self):

        sc = Algo.initialize(self)
        if sc.isFailure():
            return sc

        sc = self.  fill_initialize()
        if sc.isFailure():
            return sc

        self.c2dtf = DTF_CHI2NDOF(True, strings(['Lambda_c+']))
        self.mfit = DTF_FUN(
            M,                     True, strings(['Lambda_c+']))
        self.ctau = DTF_CTAU(
            0,                     True, strings(['Lambda_c+']))
        self.ctauC = DTF_CTAU(
            'Lambda_c+' == ABSID, True, strings(['Lambda_c+']))
        self.ctauCs = DTF_CTAUSIGNIFICANCE(
            'Lambda_c+' == ABSID, True, strings(['Lambda_c+']))

        return SUCCESS

    def finalize(self):

        self.c2dtf = None
        self.mfit = None
        self.ctau = None
        self.ctauC = None
        self.ctauCs = None

        self.  fill_finalize()

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

        xicc = self.select('xic', self.mode)

        tup = self.nTuple('Xic')
        for xi in xicc:
            #

            lc = xi(1)
            k = xi(2)
            pi = xi(3)
            pi2 = xi(4)

            self.treatKine(tup, xi, '_xi')
            self.treatKine(tup, lc, '_lc')

            self.treatKine(tup, k, '_k')
            self.treatKine(tup, pi, '_pi')

            if pi2:
                self.treatKine(tup, pi2, '_pi2')

            tup.column('mass', self.mfit(xi) / GeV)
            tup.column('c2dtf', self.c2dtf(xi))
            tup.column('ctau', self.ctau(xi))
            tup.column('ctauC', self.ctauC(xi))
            tup.column('ctauCs', self.ctauCs(xi))

            #
            # full combinatoric of masses
            #
            self.fillMasses(tup, xi)

            #
            # add the information for Pid efficiency correction
            #
            self.treatPions(tup, xi)
            self.treatKaons(tup, xi)
            self.treatMuons(tup, xi)
            self.treatProtons(tup, xi)
            self.treatTracks(tup, xi)

            #
            # add some reco-summary information
            #
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

    rootInTES = '/Event/Charm'

    from PhysConf.Filters import LoKi_Filters
    fltrs = LoKi_Filters(
        STRIP_Code="""
        HLT_PASS_RE ( 'Stripping.*Xicc.*ForPromptCharm.*' ) 
        """
    )

    davinci = DaVinci(
        DataType=the_year,
        EventPreFilters=fltrs.filters('Filters'),
        InputType='MDST',
        RootInTES=rootInTES,
        Simulation=False,
        PrintFreq=1000,
        EvtMax=-1,
        #
        HistogramFile='Xic_Histos.root',
        TupleFile='Xic.root',
        #
        Lumi=True,
        #
    )

    from BenderTools.Utils import silence
    silence()

    from Configurables import TrackScaleState
    state_scale = TrackScaleState(
        'StateScale',
    )

    davinci.UserAlgorithms = [state_scale, 'Xi1', 'Xi2']

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
    alg1 = Xicc(
        'Xi1',  # Algorithm name ,
        RootInTES=rootInTES,
        Inputs=['Phys/Xicc+ForPromptCharm/Particles']
    )

    alg2 = Xicc(
        'Xi2',  # Algorithm name ,
        RootInTES=rootInTES,
        Inputs=['Phys/Xicc++ForPromptCharm/Particles']
    )

    alg2.mode = " [ Xi_cc++ -> Lambda_c+ (K-|K+) (pi+|pi-) (pi+|pi-)]CC"

    return SUCCESS

# =============================================================================
# The actual job steering
if '__main__' == __name__:

    input = [
        '/lhcb/LHCb/Collision12/CHARM.MDST/00024181/0000/00024181_0000%04d_1.charm.mdst' % i for i in range(1, 999)
    ]

    configure(input, castor=True)

    run(5000)

    gaudi = appMgr()


# =============================================================================
# The END
# =============================================================================
