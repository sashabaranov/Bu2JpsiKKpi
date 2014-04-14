#!/usr/bin/env python
# ========================================================================
# @file DiCharm/WG3.py
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

Seek for new peaks in B&Q/WG-v3 

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
# @class B2Keta
#  Simple template algorithm to study B -> J/psi K (K) ( eta -> gamma gamma )
#  @date   2011-05-27
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class B2Keta (Algo):

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
        self.mode1K = '[ Beauty -> ( Meson -> mu+ mu- ) K+     eta  ]CC'
        self.mode2K = '  Beauty -> ( Meson -> mu+ mu- ) K+ K-  eta     '

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
        masses = ['J/psi(1S)', 'psi(2S)', 'eta', 'eta_prime', 'pi0']
        constraints = strings(masses)

        self.mfit = DTF_FUN(M, True, constraints)
        self.c2dtf = DTF_CHI2NDOF(True, constraints)
        self.ctau = DTF_CTAU(0, True, constraints)

        self.Print('Mode1K: %s ' % self.mode1K)
        self.Print('Mode2K: %s ' % self.mode2K)

        return SUCCESS

    # finalize the algorithm
    def finalize(self):
        """
        Finalize the action 
        """
        self.mfit = None
        self.c2dtf = None
        self.ctau = None

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

        b1K = self.select('b1K', self.mode1K)
        b2K = self.select('b2K', self.mode2K)

        tup1K = self.nTuple('B1K')
        tup2K = self.nTuple('B2K')

        for b in b1K:
            #
            m = M(b) / GeV
            if not 4.9 < m < 5.6:
                continue
            #
            psi = b(1)
            k = b(2)
            eta = b(3)

            mpsi = M(psi)
            if 4.5 * GeV < mpsi:
                continue
            elif 3.3 * GeV < mpsi:
                psi.setParticleID(LHCb.ParticleID(100443))

            chi2_dtf = self.c2dtf(b)
            if not 0 <= chi2_dtf < 1000:
                continue

            mass = self.mfit(b) / GeV
            if not 5.1 <= mass < 5.5:
                continue

            tup1K.column_float('c2dtf', chi2_dtf)
            tup1K.column_float('mass', mass)

            self.treatKine(tup1K, b, '_b')
            self.treatKine(tup1K, psi, '_psi')
            self.treatKine(tup1K, k, '_k')
            self.treatKine(tup1K, eta, '_x0')

            tup1K.column_float('ctau', self.ctau(b))

            #
            # fictive K**
            #
            p = LoKi.LorentzVector()
            p += k.momentum()
            p += eta.momentum()

            tup1K.column_float('m_kst', p.M() / GeV)
            tup1K.column_float('p_kst', p.P() / GeV)
            tup1K.column_float('pt_kst', p.Pt() / GeV)

            #
            # full combinatoric of masses
            #
            self.fillMasses(tup1K, b)

            # add the information needed for TisTos
            self.tisTos(psi, tup1K, 'psi_',
                        self.lines['psi'],
                        self.l0tistos,
                        self.tistos)

            # add the information for Pid efficiency correction
            self.treatPions(tup1K, b)
            self.treatMuons(tup1K, b)
            self.treatKaons(tup1K, b)
            self.treatPhotons(tup1K, b)
            self.treatDiGamma(tup1K, b)
            # add the infomation needed for track efficiency correction
            self.treatTracks(tup1K, b)

            # add some reco-summary information
            self.addRecSummary(tup1K, rc_summary)
            #
            tup1K.write()

        # check the decays with two kaons
        for b in b2K:
            #
            m = M(b) / GeV
            if not 4.9 < m < 5.6:
                continue
            #
            psi = b(1)
            k1 = b(2)
            k2 = b(3)
            eta = b(4)

            mpsi = M(psi)
            if 4.5 * GeV < mpsi:
                continue
            elif 3.3 * GeV < mpsi:
                psi.setParticleID(LHCb.ParticleID(100443))

            chi2_dtf = self.c2dtf(b)
            if not 0 <= chi2_dtf < 1000:
                continue

            mass = self.mfit(b) / GeV
            if not 5.1 <= mass < 5.5:
                continue

            tup2K.column_float('c2dtf', chi2_dtf)
            tup2K.column_float('mass', mass)

            self.treatKine(tup2K, b, '_b')
            self.treatKine(tup2K, psi, '_psi')
            self.treatKine(tup2K, k1, '_k1')
            self.treatKine(tup2K, k2, '_k2')
            self.treatKine(tup2K, eta, '_x0')

            tup2K.column_float('ctau', self.ctau(b))

            #
            # fictive phi
            #
            p = LoKi.LorentzVector()
            p += k1.momentum()
            p += k2.momentum()

            tup2K.column_float('m_KK', p.M() / GeV)
            tup2K.column_float('p_KK', p.P() / GeV)
            tup2K.column_float('pt_KK', p.Pt() / GeV)

            p += eta.momentum()

            tup2K.column_float('m_kst', p.M() / GeV)
            tup2K.column_float('p_kst', p.P() / GeV)
            tup2K.column_float('pt_kst', p.Pt() / GeV)

            #
            # full combinatoric of masses
            #
            self.fillMasses(tup2K, b)

            # add the information needed for TisTos
            self.tisTos(psi, tup2K, 'psi_',
                        self.lines['psi'],
                        self.l0tistos,
                        self.tistos)

            # add the information for Pid efficiency correction
            self.treatPions(tup2K, b)
            self.treatMuons(tup2K, b)
            self.treatKaons(tup2K, b)
            self.treatPhotons(tup2K, b)
            self.treatDiGamma(tup2K, b)
            # add the infomation needed for track efficiency correction
            self.treatTracks(tup2K, b)

            # add some reco-summary information
            self.addRecSummary(tup2K, rc_summary)
            #
            tup2K.write()

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

    rootInTES = '/Event/PSIX0'

    davinci = DaVinci(
        DataType=the_year,
        InputType='DST',
        Simulation=False,
        PrintFreq=1000,
        EvtMax=-1,
        #
        HistogramFile='WGv4_Histos.root',
        TupleFile='WGv4.root',
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

    #
    # come back to Bender
    #
    setData(datafiles, catalogs, castor)

    #
    # start Gaudi
    #
    gaudi = appMgr()

    alg1 = B2Keta(
        'Keta2g',
        RootInTES=rootInTES,
        Inputs=['Phys/SelB2PsiKEta2ggForPsiX0/Particles']
    )

    alg2 = B2Keta(
        'Keta3pi',
        RootInTES=rootInTES,
        Inputs=['Phys/SelB2PsiKEta23piForPsiX0/Particles']
    )

    alg3 = B2Keta(
        'Ketap2rg',
        RootInTES=rootInTES,
        Inputs=['Phys/SelB2PsiKEtap2rhogForPsiX0/Particles']
    )

    alg4 = B2Keta(
        'Ketap2ppe',
        RootInTES=rootInTES,
        Inputs=['Phys/SelB2PsiKEtap2pipietaForPsiX0/Particles']
    )

    alg5 = B2Keta(
        'Komega',
        RootInTES=rootInTES,
        Inputs=['Phys/SelB2PsiKOmegaForPsiX0/Particles']
    )

    #
    # kaons
    #

    alg3.mode1K = '[ Beauty -> ( Meson -> mu+ mu- ) K+     eta_prime  ]CC'
    alg3.mode2K = '  Beauty -> ( Meson -> mu+ mu- ) K+ K-  eta_prime     '

    alg4.mode1K = '[ Beauty -> ( Meson -> mu+ mu- ) K+     eta_prime  ]CC'
    alg4.mode2K = '  Beauty -> ( Meson -> mu+ mu- ) K+ K-  eta_prime     '

    alg5.mode1K = '[ Beauty -> ( Meson -> mu+ mu- ) K+     omega(782) ]CC'
    alg5.mode2K = '  Beauty -> ( Meson -> mu+ mu- ) K+ K-  omega(782)    '

    alg6 = B2Keta(
        'Pieta2g',
        RootInTES=rootInTES,
        Inputs=['Phys/SelB2PsiPiEta2ggForPsiX0/Particles']
    )

    alg7 = B2Keta(
        'Pieta3pi',
        RootInTES=rootInTES,
        Inputs=['Phys/SelB2PsiPiEta23piForPsiX0/Particles']
    )

    alg8 = B2Keta(
        'Pietap2rg',
        RootInTES=rootInTES,
        Inputs=['Phys/SelB2PsiPiEtap2rhogForPsiX0/Particles']
    )

    alg9 = B2Keta(
        'Pietap2ppe',
        RootInTES=rootInTES,
        Inputs=['Phys/SelB2PsiPiEtap2pipietaForPsiX0/Particles']
    )

    alg10 = B2Keta(
        'Piomega',
        RootInTES=rootInTES,
        Inputs=['Phys/SelB2PsiPiOmegaForPsiX0/Particles']
    )

    alg6 .mode1K = '[ Beauty -> ( Meson -> mu+ mu- ) pi+     eta        ]CC'
    alg6 .mode2K = '  Beauty -> ( Meson -> mu+ mu- ) pi+ pi- eta           '

    alg7 .mode1K = '[ Beauty -> ( Meson -> mu+ mu- ) pi+     eta        ]CC'
    alg7 .mode2K = '  Beauty -> ( Meson -> mu+ mu- ) pi+ pi- eta           '

    alg8 .mode1K = '[ Beauty -> ( Meson -> mu+ mu- ) pi+     eta_prime  ]CC'
    alg8 .mode2K = '  Beauty -> ( Meson -> mu+ mu- ) pi+ pi- eta_prime     '

    alg9 .mode1K = '[ Beauty -> ( Meson -> mu+ mu- ) pi+     eta_prime  ]CC'
    alg9 .mode2K = '  Beauty -> ( Meson -> mu+ mu- ) pi+ pi- eta_prime     '

    alg10.mode1K = '[ Beauty -> ( Meson -> mu+ mu- ) pi+     omega(782) ]CC'
    alg10.mode2K = '  Beauty -> ( Meson -> mu+ mu- ) pi+ pi- omega(782)    '

    mainSeq = gaudi.algorithm('GaudiSequencer/DaVinciUserSequence', True)
    mainSeq.Members += [alg1  . name(),
                        alg2  . name(),
                        alg3  . name(),
                        alg4  . name(),
                        alg5  . name(),
                        #
                        alg5  . name(),
                        alg6  . name(),
                        alg7  . name(),
                        alg8  . name(),
                        alg9  . name(),
                        alg10 . name()]

    return SUCCESS

# =============================================================================
# The actual job steering
if '__main__' == __name__:

    input = [
        '/lhcb/LHCb/Collision12/PSIX0.MDST/00024672/0000/00024672_00000002_1.psix0.mdst',
        '/lhcb/LHCb/Collision12/PSIX0.MDST/00024672/0000/00024672_00000003_1.psix0.mdst',
        '/lhcb/LHCb/Collision12/PSIX0.MDST/00024672/0000/00024672_00000006_1.psix0.mdst',
        '/lhcb/LHCb/Collision12/PSIX0.MDST/00024672/0000/00024672_00000009_1.psix0.mdst',
        '/lhcb/LHCb/Collision12/PSIX0.MDST/00024672/0000/00024672_00000045_1.psix0.mdst',
        '/lhcb/LHCb/Collision12/PSIX0.MDST/00024672/0000/00024672_00000047_1.psix0.mdst',
        '/lhcb/LHCb/Collision12/PSIX0.MDST/00024672/0000/00024672_00000049_1.psix0.mdst',
        '/lhcb/LHCb/Collision12/PSIX0.MDST/00024672/0000/00024672_00000052_1.psix0.mdst',
        '/lhcb/LHCb/Collision12/PSIX0.MDST/00024672/0000/00024672_00000053_1.psix0.mdst',
        '/lhcb/LHCb/Collision12/PSIX0.MDST/00024672/0000/00024672_00000055_1.psix0.mdst',
        '/lhcb/LHCb/Collision12/PSIX0.MDST/00024672/0000/00024672_00000058_1.psix0.mdst',
        '/lhcb/LHCb/Collision12/PSIX0.MDST/00024672/0000/00024672_00000059_1.psix0.mdst',
        '/lhcb/LHCb/Collision12/PSIX0.MDST/00024672/0000/00024672_00000060_1.psix0.mdst',
        '/lhcb/LHCb/Collision12/PSIX0.MDST/00024672/0000/00024672_00000061_1.psix0.mdst',
        '/lhcb/LHCb/Collision12/PSIX0.MDST/00024672/0000/00024672_00000063_1.psix0.mdst',
        '/lhcb/LHCb/Collision12/PSIX0.MDST/00024672/0000/00024672_00000064_1.psix0.mdst',
        '/lhcb/LHCb/Collision12/PSIX0.MDST/00024672/0000/00024672_00000065_1.psix0.mdst',
        '/lhcb/LHCb/Collision12/PSIX0.MDST/00024672/0000/00024672_00000067_1.psix0.mdst',
        '/lhcb/LHCb/Collision12/PSIX0.MDST/00024672/0000/00024672_00000068_1.psix0.mdst',
        '/lhcb/LHCb/Collision12/PSIX0.MDST/00024672/0000/00024672_00000069_1.psix0.mdst',
        '/lhcb/LHCb/Collision12/PSIX0.MDST/00024672/0000/00024672_00000071_1.psix0.mdst',
        '/lhcb/LHCb/Collision12/PSIX0.MDST/00024672/0000/00024672_00000073_1.psix0.mdst',
        '/lhcb/LHCb/Collision12/PSIX0.MDST/00024672/0000/00024672_00000074_1.psix0.mdst',
        '/lhcb/LHCb/Collision12/PSIX0.MDST/00024672/0000/00024672_00000075_1.psix0.mdst',
        '/lhcb/LHCb/Collision12/PSIX0.MDST/00024672/0000/00024672_00000076_1.psix0.mdst',
        '/lhcb/LHCb/Collision12/PSIX0.MDST/00024672/0000/00024672_00000078_1.psix0.mdst',
        '/lhcb/LHCb/Collision12/PSIX0.MDST/00024672/0000/00024672_00000079_1.psix0.mdst'
    ]

    configure(input, castor=True)

    run(5000)

    gaudi = appMgr()

# =============================================================================
# The END
# =============================================================================
