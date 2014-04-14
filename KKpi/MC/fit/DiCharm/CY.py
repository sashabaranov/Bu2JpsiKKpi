#!/usr/bin/env python
# ========================================================================
# @file DiCharm/CY.py
#
#  Associative production open charm + Y
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
#  @date   2013-07-11
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#
#                    $Revision$
#  Last modification $Date$
#                 by $Author$
# =============================================================================
"""

Associative production open charm + Y 

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
# @class CY
#  Simple template algorithm to study open charm + Y
#  @date   2011-05-27
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class CY (Algo):

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
        self._mode = '[Meson -> ( Meson -> mu+ mu- ) Charm ]CC'

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
        stable_charm = ['D0', 'D+', 'D_s+', 'Lambda_c+']
        excited_charm = []
        upsilons = ['Upsilon(1S)', 'Upsilon(2S)', 'Upsilon(3S)']

        all_constraints = strings(stable_charm + upsilons)

        self.mfit = DTF_FUN(M, True, all_constraints) / GeV
        self.c2dtf = DTF_CHI2NDOF(True, all_constraints)

        self.mfit_0 = DTF_FUN(M, True) / GeV
        self.c2dtf_0 = DTF_CHI2NDOF(True)

        self.mfit_1 = DTF_FUN(M, True, strings(upsilons)) / GeV
        self.c2dtf_1 = DTF_CHI2NDOF(True, strings(upsilons))

        self.mfit_2 = DTF_FUN(M, True, strings(stable_charm)) / GeV
        self.c2dtf_2 = DTF_CHI2NDOF(True, strings(stable_charm))

        self.ctauC = DTF_CTAU(strings(stable_charm) == ID, True)

        _pid = LoKi.Particles.pidFromName
        self.dimu_pids = [
            (ADMASS('J/psi(1S)'),  _pid('J/psi(1S)')),
            (ADMASS('psi(2S)'),  _pid('psi(2S)')),
            (ADMASS('Upsilon(1S)'),  _pid('Upsilon(1S)')),
            (ADMASS('Upsilon(2S)'),  _pid('Upsilon(2S)')),
            (ADMASS('Upsilon(3S)'),  _pid('Upsilon(3S)')),
        ]

        self._PID_Dstar = 'D*(2010)+' == ABSID
        self._PID_Sigmac_z = 'Sigma_c0' == ABSID
        self._PID_Sigmac_pp = 'Sigma_c++' == ABSID
        self._PID_Lc_2593 = 'Lambda_c(2593)+' == ABSID
        self._PID_Lc_2625 = 'Lambda_c(2625)+' == ABSID

        return SUCCESS

    # finalize the algorithm
    def finalize(self):
        """
        Finalize the action 
        """

        self.mfit = None
        self.c2dtf = None

        self.mfit_0 = None
        self.c2dtf_0 = None

        self.mfit_1 = None
        self.c2dtf_1 = None

        self.mfit_2 = None
        self.c2dtf_2 = None

        self.ctauC = None

        self.dimu_pids = None

        self.  fill_finalize()
        self.tisTos_finalize()

        return Algo.finalize(self)

     # set new particle ID for dimuons,
    #  accoring to nearest resonance
    def dimuon_ID(self, dimu):
        """
        Set new particle ID for dimuons,
        accoring to nearest resonance     
        """
        pid = LHCb.ParticleID()
        m_ = -1
        for i in self.dimu_pids:
            dm = i[0](dimu)
            if dm < m_ or m_ < 0:
                m_ = dm
                pid = i[1]

        dimu.setParticleID(pid)

        return pid

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

        yc = self.select('CY', self._mode)

        tup = self.nTuple('CY')

        for b in yc:
            #
            Y = b(1)
            C = b(2)
            #
            # rdefine dimuon pid
            self.dimuon_ID(Y)
            #
            #
            # skip fit errors
            #
            c0 = self.c2dtf_0(b)
            if c0 <= -1:
                continue
            #
            tup.column_float('c2dtf', self.c2dtf(b))
            tup.column_float('c2dtf0', self.c2dtf_0(b))
            tup.column_float('c2dtf1', self.c2dtf_1(b))
            tup.column_float('c2dtf2', self.c2dtf_2(b))

            tup.column_float('c2dtfY', self.c2dtf(Y))
            tup.column_float('c2dtfY0', self.c2dtf_0(Y))

            tup.column_float('c2dtfC', self.c2dtf(C))
            tup.column_float('c2dtfC0', self.c2dtf_0(C))

            tup.column_float('ctauC', self.ctauC(b))

            self.treatKine(tup, b, '_CY')
            self.treatKine(tup, Y, '_Y')
            self.treatKine(tup, C, '_C')

            if self._PID_Dstar(C):
                self.treatKine(tup, C(1), '_c2')
            elif self._PID_Sigmac_z(C):
                self.treatKine(tup, C(1), '_c2')
            elif self._PID_Sigmac_pp(C):
                self.treatKine(tup, C(1), '_c2')
            elif self._PID_Lc_2593(C):
                self.treatKine(tup, C(1), '_c2')
            elif self._PID_Lc_2625(C):
                self.treatKine(tup, C(1), '_c2')

            # add the information needed for TisTos
            self.tisTos(Y, tup, 'ups_',
                        self.lines['Upsilon'],
                        self.l0tistos,
                        self.tistos)

            # add the information needed for TisTos
            self.tisTos(Y, tup, 'ups1_',
                        self.lines['Upsilon1'],
                        self.l0tistos,
                        self.tistos)

            # add the information for Pid efficiency correction
            self.treatPions(tup, b)
            self.treatMuons(tup, b)
            self.treatKaons(tup, b)
            self.treatProtons(tup, b)
            # add the infomation needed for track efficiency correction
            self.treatTracks(tup, b)

            # add some reco-summary information
            self.addRecSummary(tup, rc_summary)
            #
            tup.write()

        return SUCCESS

# =============================================================================
# @class C1
#  Simple template algorithm to study open charm + Y
#  @date   2011-05-27
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class C1 (CY):

    """
    Simple algorithm 
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

        c1 = self.select('C1', self._mode)

        tup = self.nTuple('C1')

        for b in c1:
            #
            c = b
            #
            #
            # skip fit errors
            #
            c0 = self.c2dtf_0(b)
            if c0 <= -1:
                continue
            #
            #
            tup.column_float('c2dtf', self.c2dtf(c))
            tup.column_float('c2dtf0', self.c2dtf_0(c))
            tup.column_float('mass', self.mfit_0(c))

            self.treatKine(tup, c, '_c1')

            # add the information for Pid efficiency correction
            self.treatPions(tup, b)
            self.treatKaons(tup, b)
            self.treatProtons(tup, b)
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

    rootInTES = '/Event/YC'

    davinci = DaVinci(
        DataType=the_year,
        InputType='DST',
        Simulation=False,
        PrintFreq=1000,
        EvtMax=-1,
        #
        HistogramFile='CY_Histos.root',
        TupleFile='CY.root',
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

    alg1 = CY(
        'YD0', RootInTES=rootInTES, Inputs=['Phys/Y&D0/Particles'])
    alg2 = CY(
        'YD+', RootInTES=rootInTES, Inputs=['Phys/Y&D+/Particles'])
    alg3 = CY(
        'YDs', RootInTES=rootInTES, Inputs=['Phys/Y&Ds+/Particles'])
    alg4 = CY(
        'YLc', RootInTES=rootInTES, Inputs=['Phys/Y&Lc+/Particles'])
    alg5 = CY(
        'YDst', RootInTES=rootInTES, Inputs=['Phys/Y&Dstar+/Particles'])
    alg6 = CY(
        'YSc', RootInTES=rootInTES, Inputs=['Phys/Y&Sigmac/Particles'])
    alg7 = CY(
        'YLcst', RootInTES=rootInTES, Inputs=['Phys/Y&Lcstar+/Particles'])

    # single charm
    alg8 = C1(
        'D01', RootInTES=rootInTES, Inputs=['Phys/SelD02HHForPromptCharm/Particles '])
    alg9 = C1(
        'Dp1', RootInTES=rootInTES, Inputs=['Phys/SelDForPromptCharm/Particles'])
    alg10 = C1(
        'Ds1', RootInTES=rootInTES, Inputs=['Phys/SelDsForPromptCharm/Particles'])
    alg11 = C1(
        'Lc1', RootInTES=rootInTES, Inputs=['Phys/SelLambdaCForPromptCharm/Particles'])

    alg8 ._mode = "[ D0         -> K- pi+    ]CC"
    alg9 ._mode = "[ D+        --> K- pi+ pi+]CC"
    alg10._mode = "[ D_s+      --> K- K+  pi+]CC"
    alg11._mode = "[ Lambda_c+ --> p+ K-  pi+]CC"

    mainSeq = gaudi.algorithm('GaudiSequencer/DaVinciUserSequence', True)
    mainSeq.Members += [alg1  . name(),
                        alg2  . name(),
                        alg3  . name(),
                        alg4  . name(),
                        alg5  . name(),
                        # alg6  . name () ,
                        # alg7  . name () ,
                        #
                        alg8  . name(),
                        alg9  . name(),
                        alg10 . name(),
                        alg11 . name()
                        ]

    return SUCCESS

# =============================================================================
# The actual job steering
if '__main__' == __name__:

    import shelve

    db = shelve.open('/afs/cern.ch/user/i/ibelyaev/public/LFNs/lfns.db', 'r')

    input = db['Y&Charm,mdst:2k+12,down']
    input += db['Y&Charm,mdst:2k+12,up']
    # input += db['Y&Charm,mdst:2k+11,down']
    # input += db['Y&Charm,mdst:2k+11,up'  ]

    db.close()

    configure(input, castor=True)

    run(1000)

    gaudi = appMgr()

# =============================================================================
# The END
# =============================================================================
