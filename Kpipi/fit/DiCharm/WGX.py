#!/usr/bin/env ipython
# ========================================================================
# @file WGX.py
#
#  check for B+ -> ( J/psi pi+ pi-) K+ statistics in Matt's mDSTs
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
#  @date   2013-06-11
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#
#                    $Revision$
#  Last modification $Date$
#                 by $Author$
# =============================================================================
"""

Check for B+ -> ( J/psi pi+ pi-) K+ statistics in Matt's mDSTs 

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
# @class  B2X
#  Simple  template algorithm for B -> psi X
#  @date   2013-06-11
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class B2X (Algo):

    """
    Simple algorithm 
    """

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

        #
        masses = ['J/psi(1S)', 'psi(2S)']
        constraints = strings(masses)

        self.mfit = DTF_FUN(M, True, constraints)
        self.c2dtf = DTF_CHI2NDOF(True, constraints)
        self.ctau = DTF_CTAU(0, True, constraints)
        self.m123 = DTF_FUN(MASS(1, 2, 3), True, constraints)

        return SUCCESS

    def finalize(self):

        self.mfit = None
        self.c2dtf = None
        self.ctau = None
        self.m123 = None

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

        B = self.select(
            'b', '[Beauty -> ( J/psi(1S) -> mu+ mu- ) pi+ pi- K+]CC')

        tup = self.nTuple('B')
        for b in B:
            #

            m = M(b) / GeV
            if not 5.0 < m < 5.6:
                continue
            #

            psi = b(1)
            pi1 = b(2)
            pi2 = b(3)
            k = b(4)

            mpsi = M(psi)
            if 4.5 * GeV < mpsi:
                continue
            elif 3.3 * GeV < mpsi:
                psi.setParticleID(LHCb.ParticleID(100443))

            tup.column('m', M(b) / GeV)
            tup.column('mass', self.mfit(b) / GeV)
            tup.column('c2dtf', self.c2dtf(b))
            tup.column('ctau', self.ctau(b))
            tup.column('mX', self.m123(b))

            self.treatKine(tup, b, '_B')
            self.treatKine(tup, psi, '_psi')
            self.treatKine(tup, pi1, '_pi1')
            self.treatKine(tup, pi2, '_pi2')
            self.treatKine(tup, k, '_K')

            #
            # full combinatoric of masses
            #
            self.fillMasses(tup, b)

            # add the information needed for TisTos
            self.tisTos(psi, tup, 'psi_',
                        self.lines['psi'],
                        self.l0tistos,
                        self.tistos)

            # add the information for Pid efficiency correction
            self.treatPions(tup, b)
            self.treatMuons(tup, b)
            self.treatKaons(tup, b)

            # add the infomation needed for track efficiency correction
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

    rootInTES = '/Event/B'

    davinci = DaVinci(
        DataType=the_year,
        InputType='MDST',
        Simulation=False,
        PrintFreq=1000,
        EvtMax=-1,
        #
        HistogramFile='B2X_Histos.root',
        TupleFile='B2X.root',
        #
        Lumi=True,
        #
    )

    from Configurables import CondDB
    CondDB(LatestGlobalTagByDataType=the_year)

    # ------- decoding set-up start ----------
    from BenderTools.MicroDST import uDstConf
    uDstConf(rootInTES)
    # ------- decoding set-up end  -----------

    # suppress some unnesessary prints
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

    alg = B2X(
        'B2X',  # Algorithm name ,
        RootInTES=rootInTES,
        Inputs=['Phys/Bplus/Particles']
    )

    mainSeq = gaudi.algorithm('GaudiSequencer/DaVinciUserSequence', True)
    mainSeq.Members += [alg . name()]

    return SUCCESS

# =============================================================================
# The actual job steering
if '__main__' == __name__:

    inputs = []

    import shelve
    db = shelve.open('/afs/cern.ch/user/i/ibelyaev/public/LFNs/lfns.db', 'r')

    inputs += db['B+ -> X(3872)K+, 2k+11, MagnetDown']
    # inputs += db ['B+ -> X(3872)K+, 2k+11, MagnetUp']
    # inputs += db ['B+ -> X(3872)K+, 2k+12, MagnetDown']
    # inputs += db ['B+ -> X(3872)K+, 2k+12, MagnetUp']

    # inputs = [
    #    '/lhcb/LHCb/Collision12/PSIX0.MDST/00023621/0000/00023621_00000001_1.psix0.mdst'
    #    ]

    configure(inputs[:10], castor=True)

    run(5000)

    gaudi = appMgr()

# =============================================================================
# The END
# =============================================================================
