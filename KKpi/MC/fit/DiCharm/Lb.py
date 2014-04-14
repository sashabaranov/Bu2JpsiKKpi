#!/usr/bin/env python
# ========================================================================
# @file DiCharm/Lb.py
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
# @class LB
#  Simple template algorithm to study 2xB
#  @date   2011-05-27
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class LB (Algo):

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

        constraint_Lc = ['Lambda_c+']
        constraint_Ds = ['D+', 'D_s+']

        constraints = strings(constraint_Lc + constraint_Ds)

        #
        self.mfit = DTF_FUN(M, True, constraints)
        #
        self.c2dtf = DTF_CHI2NDOF(True, constraints)
        self.c2dtf0 = DTF_CHI2NDOF(True)
        self.c2dtf1 = DTF_CHI2NDOF(True, strings(constraint_Ds))
        self.c2dtf2 = DTF_CHI2NDOF(True, 'Lambda_c+')

        self.ctau = DTF_CTAU(0, True, constraints)
        self.ctau1 = DTF_CTAU(
            'Lambda_c+' == ABSID, True, constraints)
        self.ctau2 = DTF_CTAU(
            ('D_s+' == ABSID) | ('D+' == ABSID), True, constraints)
        self.ctau1s = DTF_CTAUSIGNIFICANCE(
            'Lambda_c+' == ABSID, True, constraints)
        self.ctau2s = DTF_CTAUSIGNIFICANCE(
            ('D_s+' == ABSID) | ('D+' == ABSID), True, constraints)

        self.m1 = DTF_FUN(
            M1, True, strings(constraint_Ds))
        self.m2 = DTF_FUN(
            M2, True, strings(constraint_Lc))

        return SUCCESS

    # finalize the algorithm
    def finalize(self):
        """
        Finalize the action 
        """
        #
        self.mfit = None

        self.ctau = None

        self.ctau1 = None
        self.ctau2 = None

        self.ctau1s = None
        self.ctau2s = None
        self.m1 = None
        self.m2 = None
        #
        self.c2dtf = None
        self.c2dtf0 = None
        self.c2dtf1 = None
        self.c2dtf2 = None
        #
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

        Lb = self.select(
            'lb', ' [ Lambda_b0 -> Lambda_c+ ( D_s- | D_s+ ) ]CC')

        tup = self.nTuple('B')

        for l in Lb:
            #
            lc = l(1)
            ds = l(2)
            #
            if M(ds) < 1.9 * GeV:
                pid = ds.particleID().pid()
                if 431 == pid:
                    ds.setParticleID(LHCb.ParticleID(411))
                elif -431 == pid:
                    ds.setParticleID(LHCb.ParticleID(-411))

            self.treatKine(tup, l, '_Lb')
            self.treatKine(tup, lc, '_Lc')
            self.treatKine(tup, ds, '_D')

            tup.column_float('mass', self.mfit(l) / GeV)

            tup.column_float('c2dtf', self.c2dtf(l))
            tup.column_float('c2dtf0', self.c2dtf0(l))
            tup.column_float('c2dtf1', self.c2dtf1(l))
            tup.column_float('c2dtf2', self.c2dtf2(l))

            tup.column_float('ctau', self.ctau(l))

            tup.column_float('ctau1', self.ctau1(l))
            tup.column_float('ctau2', self.ctau2(l))

            tup.column_float('ctau1s', self.ctau1s(l))
            tup.column_float('ctau2s', self.ctau2s(l))

            tup.column_float('m1', self.m1(l) / GeV)
            tup.column_float('m2', self.m2(l) / GeV)

            # add the information for Pid efficiency correction
            self.treatPions(tup, l)
            self.treatMuons(tup, l)
            self.treatKaons(tup, l)
            self.treatProtons(tup, l)

            # add the infomation needed for track efficiency correction
            self.treatTracks(tup, l)

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

    the_year = "2012"

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

    rootInTES = '/Event/Charm'

    # Read only fired events to speed up
    from PhysConf.Filters import LoKi_Filters
    fltrs = LoKi_Filters(
        STRIP_Code="""
        HLT_PASS_RE( 'Stripping.*DsLamCForPromptCharm.*' ) 
        """
    )

    dv = DaVinci(
        EventPreFilters=fltrs.filters('Filters'),
        DataType=the_year,
        InputType='MDST',
        RootInTES=rootInTES,
        Simulation=False,
        PrintFreq=1000,
        EvtMax=-1,
        #
        HistogramFile='LB_Histos.root',
        TupleFile='LB.root',
        #
        Lumi=True,
        #
    )

    dv.UserAlgorithms = ['LB']

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
    alg = LB(
        'LB',
        RootInTES=rootInTES,
        Inputs=['Phys/DsLamCForPromptCharm/Particles']
    )
    return SUCCESS

# =============================================================================
# The actual job steering
if '__main__' == __name__:

    input = [
        '/lhcb/LHCb/Collision12/CHARM.MDST/00024181/0000/00024181_0000%04d_1.charm.mdst' % i for i in range(1, 999)
    ]

    configure(input, castor=True)

    run(10000)

    gaudi = appMgr()

# =============================================================================
# The END
# =============================================================================
