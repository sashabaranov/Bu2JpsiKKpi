#!/usr/bin/env python
# ========================================================================
# @file DiCharm/BB.py
#
#  Seek for new peaks in B&Q/WG-v4
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

Seek for new peaks in B&Q/WG-v4

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
import BenderTools.Fill  # add methods to fill n-tuple info
# =============================================================================

# =============================================================================
# @class BB
#  Simple template algorithm to study 2xB
#  @date   2011-05-27
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class BB (Algo):

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

    def initialize(self):

        sc = Algo.       initialize(self)
        if sc.isFailure():
            return sc

        triggers = {}
        from DiCharm.TisTos import lines

        sc = self.tisTos_initialize(triggers, lines)
        if sc.isFailure():
            return sc

        sc = self.  fill_initialize()
        if sc.isFailure():
            return sc

        self.mfit = DTF_FUN(M, True, strings(['J/psi(1S)', 'psi(2S)']))
        self.mfit1 = DTF_FUN(M1, True, strings(['J/psi(1S)', 'psi(2S)']))
        self.mfit2 = DTF_FUN(M2, True, strings(['J/psi(1S)', 'psi(2S)']))

        self.c2dtf = DTF_CHI2NDOF(True, strings(['J/psi(1S)', 'psi(2S)']))
        self.ctau1 = DTF_CTAU(
            1, True, strings(['J/psi(1S)', 'psi(2S)']))
        self.ctau2 = DTF_CTAU(
            1, True, strings(['J/psi(1S)', 'psi(2S)']))

        self.jpsi = LoKi.Child.Selector(
            (M < 3.3 * GeV) & DECTREE('Meson -> mu+ mu-')
        )

        self.nPi = NINTREE('pi+' == ABSID)
        self.nK = NINTREE('K+' == ABSID)

        return SUCCESS

    # finalize the algorithm
    def finalize(self):
        """
        Finalize the action 
        """
        #
        self.mfit = None
        self.mfit1 = None
        self.mfit2 = None
        #
        self.c2dtf = None
        self.ctau1 = None
        self.ctau2 = None

        self.jpsi = None

        self.nPi = None
        self.nK = None

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
            return self.Warning('No primary vertices are found', SUCCESS)

        #
        # rec summary
        rc_summary = self.get('/Event/Rec/Summary').summaryData()
        #

        bb = self.select('bb', 'Meson -> Beauty Beauty')
        if bb.empty():
            all = self.select('all', PALL)
            for a in all:
                print 'COMPONENT:', a.decay()

        tup = self.nTuple('B')

        for b in bb:
            #
            print 'B', b
            psis = LHCb.Particle.ConstVector()
            b.children(self.jpsi, psis)
            if 2 != len(psis):
                continue

            b1 = b(1)
            b2 = b(2)
            #

            self.treatKine(tup, b, '_2b')
            self.treatKine(tup, b1, '_b1')
            self.treatKine(tup, b2, '_b2')

            tup.column_float('mass', self.mfit(b) / GeV)
            tup.column_float('m1', self.mfit(b) / GeV)
            tup.column_float('m2', self.mfit(b) / GeV)
            tup.column_float('c2dtf', self.c2dtf(b))
            tup.column_float('ctau1', self.ctau1(b))
            tup.column_float('ctau2', self.ctau2(b))

            self.column_int('npi1', int(self.nPi(b1)))
            self.column_int('k1', int(self.nK(b1)))
            self.column_int('npi2', int(self.nPi(b2)))
            self.column_int('k2', int(self.nK(b2)))

            # add the information for Pid efficiency correction
            self.treatPions(tup, b)
            self.treatMuons(tup, b)
            self.treatKaons(tup, b)
            # add the infomation needed for track efficiency correction
            self.treatTracks(tup, b)

            i = 0
            for j in psis:
                # add the information needed for TisTos
                i += 1
                self.tisTos(j, tup, 'psi%d_' % i,
                            self.lines['psi'],
                            self.l0tistos,
                            self.tistos)

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

    the_year = "2011"

    from BenderTools.Parser import hasInFile

    if params:
        the_year = params['Year']
        logger.info('Year is set from params to be %s ' % the_year)
    else:
        if hasInFile(datafiles, 'Collision11'):
            the_year = '2011'
        elif hasInFile(datafiles, 'Collision12'):
            the_year = '2012'
        elif hasInFile(datafiles, 'Collision13'):
            the_year = '2013'
        elif hasInFile(datafiles, 'Stripping17'):
            the_year = '2011'
        elif hasInFile(datafiles, 'Stripping13'):
            the_year = '2011'
        elif hasInFile(datafiles, 'Stripping15'):
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
        0.5 < CONTAINS('/Event/PSIX/Phys/SelB&B/Particles' )
        """
    )

    dv = DaVinci(
        DataType=the_year,
        InputType='MDST',
        RootInTES=rootInTES,
        Simulation=False,
        PrintFreq=1000,
        EvtMax=-1,
        #
        HistogramFile='BB_Histos.root',
        TupleFile='BB.root',
        #
        Lumi=True,
        #
    )

    dv.UserAlgorithms = ['BB']

    from Configurables import CondDB
    CondDB(LatestGlobalTagByDataType=the_year)

    from BenderTools.Utils import silence
    silence()

    #
    # come back to Bender
    #
    setData(datafiles, catalogs, castor)

    #
    # start Gaudi
    #
    gaudi = appMgr()
    alg = BB(
        'BB',
        RootInTES=rootInTES,
        Inputs=['Phys/SelB&B/Particles']
    )

    return SUCCESS

# =============================================================================
# The actual job steering
if '__main__' == __name__:

    input = [
        '/lhcb/LHCb/Collision12/PSIX.MDST/00024671/0000/00024671_00000001_1.psix.mdst'
    ]

    configure(input, castor=True)

    run(10000)

    gaudi = appMgr()

# =============================================================================
# The END
# =============================================================================
