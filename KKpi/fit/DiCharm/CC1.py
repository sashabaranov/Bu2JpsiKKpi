# ========================================================================
# @file DiCharm/Bc.py
#
#  Simple template algorithm for Bc
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

Simple template algorithm for Bc

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
# @class B2Kpp
#  Simple template algorithm to study Bc -> J/psi K pi pi
#  @date   2011-05-27
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class B2Kpp (DiCharmAlg):

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

            self.Warning('No primary vertices are found', SUCCESS)
            mips = -1000 * PONE()

        else:

            mips = MIPCHI2(self.geo(), primaries)

        #
        # rec summary
        rc_summary = self.get('/Event/Rec/Summary').summaryData()
        #
        bc = self.select(
            'bc', '[ Beauty -> ( J/psi(1S) -> mu+ mu- ) K+ pi+ pi-]CC')

        c2dtf = DTF_CHI2NDOF(True, 'J/psi(1S)')
        mfit = DTF_FUN(M, True, 'J/psi(1S)')
        ctau = DTF_CTAU(0, True, strings(['J/psi(1S)', 'D_s+']))

        miniph = MINTREE(('pi+' == ABSID) | ('K+' == ABSID), BPVIPCHI2())

        tup = self.nTuple('B2Kpp')
        for b in bc:
            #
            m = M(b) / GeV
            if not 4.9 < m < 5.6:
                continue
            #
            psi = b(1)
            k1 = b(2)
            pi2 = b(3)
            pi3 = b(4)

            pi1 = k1

            self.treatKine(tup, b, '_bc')
            self.treatKine(tup, psi, '_psi')
            self.treatKine(tup, k1, '_pi1')
            self.treatKine(tup, pi2, '_pi2')
            self.treatKine(tup, pi3, '_pi3')

            tup.column('mass', mfit(b) / GeV)
            tup.column('c2dtf', c2dtf(b))
            tup.column('ctau', ctau(b))
            tup.column('minipH', miniph(b))

            #
            # fictive Ds
            p = LoKi.LorentzVector()
            p += k1.momentum()
            p += pi2.momentum()
            p += pi3.momentum()

            tup.column('m_ds', p.M() / GeV)
            tup.column('p_ds', p.P() / GeV)
            tup.column('pt_ds', p.Pt() / GeV)

            #
            #
            #

            # add the information needed for TisTos
            self.tisTos(psi, tup, 'psi_',
                        self.lines['psi'],
                        self.l0tistos,
                        self.tistos)

            # add the information for Pid efficiency correction
            self.treatPions(tup, b)
            self.treatMuons(tup, b)
            # add the infomation needed for track efficiency correction
            self.treatTracks(tup, b)
            self.treatKaons(tup, b)

            # add some reco-summary information
            self.addRecSummary(tup, rc_summary)
            #
            tup.write()

        return SUCCESS

# =============================================================================
# @class B2K
#  Simple template algorithm to study B -> J/psi K
#  @date   2011-05-27
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class B2K (DiCharmAlg):

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
            self.Warning('No primary vertices are found', SUCCESS)

        #
        # rec summary
        rc_summary = self.get('/Event/Rec/Summary').summaryData()
        #

        bc = self.select(
            'bc', '[ Beauty -> ( J/psi(1S) -> mu+ mu- ) K+ ]CC')

        c2dtf = DTF_CHI2NDOF(True, 'J/psi(1S)')
        mfit = DTF_FUN(M, True, 'J/psi(1S)')
        ctau = DTF_CTAU(0, True, strings(['J/psi(1S)', 'D_s+']))

        miniph = MINTREE(('pi+' == ABSID) | ('K+' == ABSID), BPVIPCHI2())

        doca = DOCA(1, 2, self.distanceCalculator())

        tup = self.nTuple('B2K')
        for b in bc:
            #
            m = M(b) / GeV
            if not 4.9 < m < 5.6:
                continue
            #
            psi = b(1)
            k1 = b(2)

            pi1 = k1

            self.treatKine(tup, b, '_bc')
            self.treatKine(tup, psi, '_psi')
            self.treatKine(tup, k1, '_pi1')

            tup.column('mass', mfit(b) / GeV)
            tup.column('c2dtf', c2dtf(b))
            tup.column('ctau', ctau(b))
            tup.column('minipH', miniph(b))

            #
            # fictive Ds
            #
            tup.column('m_ds', M(pi1) / GeV)
            tup.column('p_ds', P(pi1) / GeV)
            tup.column('pt_ds', PT(pi1) / GeV)

            #
            #
            #
            basic = LHCb.Particle.ConstVector()

            b.children(LoKi.Child.Selector(ISBASIC & HASTRACK), basic)

            tup.column('doca', doca.docachi2max(basic))
            #
            #
            #

            # add the information needed for TisTos
            self.tisTos(psi, tup, 'psi_',
                        self.lines['psi'],
                        self.l0tistos,
                        self.tistos)

            # add the information for Pid efficiency correction
            self.treatPions(tup, b)
            self.treatMuons(tup, b)
            # add the infomation needed for track efficiency correction
            self.treatTracks(tup, b)
            self.treatKaons(tup, b)

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
        elif hasInFile(datafiles, 'Stripping17'):
            the_year = '2011'
        elif hasInFile(datafiles, 'Stripping13'):
            the_year = '2011'
        elif hasInFile(datafiles, 'Stripping15'):
            the_year = '2011'
        elif hasInFile(datafiles, 'Stripping19'):
            the_year = '2012'
        elif hasInFile(datafiles, 'Stripping20'):
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

    #
    Jpsi_location = '/Event/Dimuon/Phys/FullDSTDiMuonJpsi2MuMuDetachedLine/Particles'
    #

    # Read only fired events to speed up
    from PhysConf.Filters import LoKi_Filters
    fltrs = LoKi_Filters(
        STRIP_Code="HLT_PASS_RE('Stripping.*DiMuonJpsi2MuMuDeta.*')",
        VOID_Code="""
        0.5 < CONTAINS('%s')
        """ % Jpsi_location
    )
    #
    # protection agains ``corrupted'' Stripping 17b DIMUON.DST
    fltrs_0 = LoKi_Filters(
        VOID_Code="""
        ( EXISTS ( '/Event/DAQ/RawEvent') | EXISTS('/Event/Trigger/RawEvent' ) ) 
        & EXISTS ( '/Event/Strip/Phys/DecReports') 
        """
    )

    davinci = DaVinci(
        EventPreFilters=fltrs_0.filters(
            'Filters0') + fltrs.filters('Filters'),
        DataType=the_year,
        InputType='DST',
        Simulation=False,
        PrintFreq=1000,
        EvtMax=-1,
        #
        HistogramFile='Bcc1_Histos.root',
        TupleFile='Bcc1.root',
        #
        Lumi=True,
        #
    )

    from Configurables import CondDB
    CondDB(LatestGlobalTagByDataType=the_year)

    # ------- decoding set-up start ----------
    ## from BenderTools.MicroDST import uDstConf
    ## uDstConf ( rootInTES )
    # ------- decoding set-up end  -----------

    #
    # dimuon locations in DIMUON.DST
    #
    from PhysSelPython.Wrappers import AutomaticData
    jpsi = AutomaticData(Location=Jpsi_location)
    #
    # get the prompt charm
    #
    from StrippingSelections.StrippingPromptCharm import StrippingPromptCharmConf as PC
    #
    # ======================================
    pc = PC('PromptCharm', {
        'TrackCuts'       : """
        ( TRCHI2DOF   < 4   ) &
        ( TRGHOSTPROB < 0.5 ) &           
        ( PT > 250 * MeV    ) &
        in_range  ( 2 , ETA , 5 ) 
        """ ,
        'KaonCuts': ' & in_range ( 3.2 * GeV , P , 100 * GeV ) & ( 2 < PIDK  - PIDpi ) ',
        'PionCuts': ' & in_range ( 3.2 * GeV , P , 100 * GeV ) & ( 0 < PIDpi - PIDK  ) ',
    }
    )

    pions = pc.pions()
    kaons = pc.kaons()

    Preambulo = [
        # shortcut for chi2 of vertex fit
        'chi2vx = VFASPF(VCHI2) ',
        # shortcut for the c*tau
        "from GaudiKernel.PhysicalConstants import c_light",
        # use the embedded cut for chi2(LifetimeFit)<9 !!!
        "ctau   = BPVLTIME ( 9 ) * c_light "  # ATTENTION, 9 is here!
    ]

    from GaudiConfUtils.ConfigurableGenerators import CombineParticles

    # ========================================================================
    # B -> J/psi + K pi pi
    # ========================================================================
    from GaudiConfUtils.ConfigurableGenerators import CombineParticles
    bc_Kpp = CombineParticles(
        DecayDescriptor='[B+ -> J/psi(1S) K+ pi+ pi-]cc',
        #
        Preambulo=Preambulo,
        DaughtersCuts={
            "J/psi(1S)": " in_range( 3.096 * GeV - 45 * MeV , M , 3.096 * GeV + 45 * MeV ) "
        },
        #
        CombinationCut="""
        in_range ( 5.0 * GeV , AM ,  5.6 * GeV ) 
        """ ,
        #
        MotherCut="""
        in_range  ( 5.1 * GeV , M , 5.5 * GeV ) &
        ( PT      > 1 * GeV          ) &
        ( chi2vx  <  49              ) &
        in_range ( 150 * micrometer , ctau , 1000 * micrometer ) 
        """ ,
        #
        ParticleCombiners={'': 'LoKi::VertexFitter'},
        ReFitPVs=True
    )
    #
    from PhysSelPython.Wrappers import Selection
    Bc_Kpp = Selection(
        'PsiKpp',
        Algorithm=bc_Kpp,
        RequiredSelections=[jpsi, pions, kaons]
    )

    # ========================================================================
    # B -> J/psi + K
    # ========================================================================
    from GaudiConfUtils.ConfigurableGenerators import CombineParticles
    bc_K = CombineParticles(
        DecayDescriptor='[ B+ -> J/psi(1S) K+ ]cc',
        #
        Preambulo=Preambulo,
        DaughtersCuts={
            "J/psi(1S)": " in_range( 3.096 * GeV - 45 * MeV , M , 3.096 * GeV + 45 * MeV ) "
        },
        #
        CombinationCut="""
        in_range ( 5.0 * GeV , AM ,  5.6 * GeV ) 
        """ ,
        #
        MotherCut="""
        in_range  ( 5.1 * GeV , M , 5.5 * GeV ) &
        ( PT      > 1 * GeV          ) &
        ( chi2vx  <  16              ) &
        in_range ( 150 * micrometer , ctau , 1000 * micrometer ) 
        """ ,
        ParticleCombiners={'': 'LoKi::VertexFitter'},
        ReFitPVs=True
    )
    #
    Bc_K = Selection(
        'PsiK',
        Algorithm=bc_K,
        RequiredSelections=[jpsi, kaons]
    )

    from PhysSelPython.Wrappers import SelectionSequence
    Seq_Kpp = SelectionSequence("PSIKPP", TopSelection=Bc_Kpp)
    Seq_K = SelectionSequence("PSIK", TopSelection=Bc_K)

    from Configurables import GaudiSequencer
    davinci.UserAlgorithms = [
        GaudiSequencer(
            'K', Members=[Seq_K  .sequence(), 'B2PsiK']),
        GaudiSequencer('KPP', Members=[Seq_Kpp.sequence(), 'B2PsiKpp'])
    ]

    #
    # come back to Bender
    #
    setData(datafiles, catalogs, castor)

    #
    # start Gaudi
    #
    gaudi = appMgr()

    #
    algKpp = B2Kpp(
        'B2PsiKpp',  # Algorithm name ,
        Inputs=[Seq_Kpp.outputLocation()]
    )
    algK = B2K(
        'B2PsiK',  # Algorithm name ,
        Inputs=[Seq_K  .outputLocation()]
    )

    return SUCCESS

# =============================================================================
# The actual job steering
if '__main__' == __name__:

    input = [
        '/lhcb/LHCb/Collision12/DIMUON.DST/00018987/0000/00018987_00000094_1.dimuon.dst',
        '/lhcb/LHCb/Collision12/DIMUON.DST/00018987/0000/00018987_00000099_1.dimuon.dst',
        '/lhcb/LHCb/Collision12/DIMUON.DST/00018987/0000/00018987_00000108_1.dimuon.dst',
        '/lhcb/LHCb/Collision12/DIMUON.DST/00018987/0000/00018987_00000109_1.dimuon.dst',
        '/lhcb/LHCb/Collision12/DIMUON.DST/00018987/0000/00018987_00000110_1.dimuon.dst',
        '/lhcb/LHCb/Collision12/DIMUON.DST/00018987/0000/00018987_00000777_1.dimuon.dst'
    ]

    configure(input, castor=True)

    run(5000)

    gaudi = appMgr()

# =============================================================================
# The END
# =============================================================================
