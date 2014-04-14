#!/usr/bin/env ipython
# ========================================================================
# @file DiCharm/MCBc.py
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
from Bender.MainMC import *  # import all bender goodies
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
import BenderTools.TisTosMC  # add methods for TisTos
import DiCharm.Pids  # add methods for Pid/Track information
# =============================================================================

# =============================================================================
# @class MCBc
#  Simple template algorithm to study Bc
#  @date   2011-05-27
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class MCBc (AlgoMC):

    """
    Simple algorithm 
    """

    def initialize(self):
        """
        Initialization
        """
        sc = AlgoMC.    initialize(self)
        if sc.isFailure():
            return sc

        triggers = {}
        triggers['J/psi(1S)'] = {}
        from DiCharm.TisTos import lines

        sc = self.tisTos_initialize(triggers, lines)
        if sc.isFailure():
            return sc

        sc = self.   pid_initialize()
        if sc.isFailure():
            return sc

        return SUCCESS

    # finalize it!
    def finalize(self):
        """
        finalize the algorithm
        """
        #
        self . tisTos_finalize()
        self .    pid_finalize()
        #
        self . dumpHistos()
        #
        return AlgoMC.finalize(self)

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

        bc = self.  select(
            'bc', '[ B_c+ -> ( J/psi(1S) -> mu+ mu- ) ( D_s+ --> K+ K- pi+ ) ]CC')

        if bc.empty():
            return SUCCESS

        mcbc = self.mcselect(
            'mcbc', '[ B_c+ => ( J/psi(1S) => mu+ mu- ) ( D_s+ ==> K+ K- pi+ ) ]CC')
        if mcbc.empty():
            self.Warning('No true MC-decays are found!', SUCCESS)
            mcbc = self.mcselect('mcbcall', MCABSID == 'B_c+')
            print mcbc
            return SUCCESS

        mcBc = MCTRUTH(self.mcTruth(), mcbc)

        c2dtf = DTF_CHI2NDOF(True, strings(['J/psi(1S)', 'D_s+']))
        mfit = DTF_FUN(M, True, strings(['J/psi(1S)', 'D_s+']))
        ctau = DTF_CTAU(0, True, strings(['J/psi(1S)', 'D_s+']))

        ctau_9 = BPVLTIME(9) * c_light
        ipchi2 = BPVIPCHI2()

        ctauC = DTF_CTAU(
            'D_s+' == ABSID, True, strings(['J/psi(1S)', 'D_s+']))
        ctauCs = DTF_CTAUSIGNIFICANCE(
            'D_s+' == ABSID, True, strings(['J/psi(1S)', 'D_s+']))

        maxTrGH = MAXTREE(ISBASIC & HASTRACK, TRGHOSTPROB)
        maxTrC2 = MAXTREE(ISBASIC & HASTRACK, TRCHI2DOF)

        minDllK = MINTREE('K+' == ABSID, PIDK - PIDpi)
        minDllmu = MINTREE('mu+' == ABSID, PIDmu - PIDpi)
        minDllpi = MINTREE('pi+' == ABSID, PIDpi - PIDK)
        minPTmu = MINTREE('mu+' == ABSID, PT)

        tup = self.nTuple('Bc')
        for b in bc:
            #
            m = M(b) / GeV
            if not m < 7.1:
                continue
            #
            psi = b(1)
            ds = b(2)

            tup.column('mc', mcBc(b))

            tup.column('mass', mfit(b) / GeV)
            tup.column('c2dtf', c2dtf(b))

            tup.column('ctau', ctau(b))

            tup.column('ctau_bc', ctau_9(b))
            tup.column('ipchi2_bc', ipchi2(b))

            tup.column('ctauC', ctauC(b))
            tup.column('ctauCs', ctauCs(b))

            tup.column('mKK', M12(ds) / GeV)

            tup.column('m_psi', M(psi) / GeV)
            tup.column('m_ds', M(ds) / GeV)
            tup.column('pt_psi', PT(psi) / GeV)
            tup.column('pt_ds', PT(ds) / GeV)
            tup.column('pt_bc', PT(b) / GeV)

            tup.column('pt_mu', minPTmu(b) / GeV)

            tup.column('maxChi2_track', maxTrC2(b))
            tup.column('maxTrGh_track', maxTrGH(b))
            tup.column('mindll_mu', minDllmu(b))
            tup.column('mindll_K', minDllK(b))
            tup.column('mindll_piK', minDllpi(b))

            # add the information needed for TisTos
            self.tisTos(psi, tup, 'psi_',
                        self.lines['psi'],
                        self.l0tistos,
                        self.tistos)
            self.tisTos(psi, tup, 'psicl_',
                        self.lines['psi_clean'],
                        self.l0tistos,
                        self.tistos)
# self.decisions ( psi                        ,
##                              self.triggers['J/psi(1S)'] ,
##                              self.l0tistos              ,
# self.  tistos              )

# add the information for Pid efficiency correction
##             self.treatPions   ( tup , b )
##             self.treatKaons   ( tup , b )
##             self.treatMuons   ( tup , b )
##             self.treatTracks  ( tup , b )

# add some reco-summary information
##             self.addRecSummary ( tup , rc_summary   )
            #
            tup.write()

        return SUCCESS


# =============================================================================
# @class MCBc1Pi
#  Simple template algorithm to study Bc -> J/psi pi
#  @date   2011-05-27
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
class MCBc1Pi (MCBc):

    """
    Simple algorithm 
    """
    # the only one esential method:

    def analyse(self):
        """
        The only one essential method
        """

##         mct = MCDECTREE( ' [ Lambda0 -> p+ pi-]CC')
##         mcl  = self.mcselect ( 'l'  , '[ Lambda0 =>  p+  pi-]CC' )
##         mcp  = self.mcselect ( 'p'  , '[ Lambda0 => ^p+  pi-]CC' )
##         mcpi = self.mcselect ( 'pi' , '[ Lambda0 =>  p+ ^pi-]CC' )
# print '# of lambdas:', len(mcl), len(mcp), len(mcpi)
# for l in mcl :
# print 'DECAY', l.decay()
# if not mct ( l ) :
# print l
# return SUCCESS
        primaries = self.vselect('PVs', ISPRIMARY)
        if primaries.empty():
            self.Warning('No primary vertices are found', SUCCESS)

        #
        # rec summary
        rc_summary = self.get('/Event/Rec/Summary').summaryData()
        #

        bc = self.  select(
            'bc', '[ B_c+ -> ( J/psi(1S) -> mu+ mu- ) pi+ ]CC')

        if bc.empty():
            return SUCCESS

        mcbc = self.mcselect(
            'mcbc', '[ B_c+ => ( J/psi(1S) => mu+ mu- ) pi+ ]CC')
        if mcbc.empty():
            self.Warning('No true MC-decays are found!', SUCCESS)
            mcbc = self.mcselect('mcbcall', MCABSID == 'B_c+')
            print mcbc
            return SUCCESS

        mcBc = MCTRUTH(self.mcTruth(), mcbc)

        c2dtf = DTF_CHI2NDOF(True, strings(['J/psi(1S)']))
        mfit = DTF_FUN(M, True, strings(['J/psi(1S)']))
        ctau = DTF_CTAU(0, True, strings(['J/psi(1S)']))

        ctau_9 = BPVLTIME(9) * c_light
        ipchi2 = BPVIPCHI2()

        maxTrGH = MAXTREE(ISBASIC & HASTRACK, TRGHOSTPROB)
        maxTrC2 = MAXTREE(ISBASIC & HASTRACK, TRCHI2DOF)

        minDllK = MINTREE('K+' == ABSID, PIDK - PIDpi)
        minDllmu = MINTREE('mu+' == ABSID, PIDmu - PIDpi)
        minDllpi = MINTREE('pi+' == ABSID, PIDpi - PIDK)
        minPTmu = MINTREE('mu+' == ABSID, PT)

        tup = self.nTuple('Bc')
        for b in bc:
            #
            m = M(b) / GeV
            if not m < 7.1:
                continue
            #
            psi = b(1)
            pion = b(2)

            tup.column('mc', mcBc(b))

            tup.column('mass', mfit(b) / GeV)
            tup.column('c2dtf', c2dtf(b))
            tup.column('ctau', ctau(b))

            tup.column('ctau_bc', ctau_9(b))
            tup.column('ipchi2_bc', ipchi2(b))

            tup.column('m_psi', M(psi) / GeV)
            tup.column('pt_psi', PT(psi) / GeV)
            tup.column('pt_bc', PT(b) / GeV)

            tup.column('m_pi', M(pion) / GeV)
            tup.column('pt_pi', PT(pion) / GeV)

            ds = pion
            tup.column('m_ds', M(ds) / GeV)
            tup.column('pt_ds', PT(ds) / GeV)

            tup.column('pt_mu', minPTmu(b) / GeV)

            tup.column('maxChi2_track', maxTrC2(b))
            tup.column('maxTrGh_track', maxTrGH(b))
            tup.column('mindll_mu', minDllmu(b))
            tup.column('mindll_K', minDllK(b))
            tup.column('mindll_piK', minDllpi(b))

            # add the information needed for TisTos
            self.tisTos(psi, tup, 'psi_',
                        self.lines['psi'],
                        self.l0tistos,
                        self.tistos)
            self.tisTos(psi, tup, 'psicl_',
                        self.lines['psi_clean'],
                        self.l0tistos,
                        self.tistos)
# self.decisions ( psi                        ,
##                              self.triggers['J/psi(1S)'] ,
##                              self.l0tistos              ,
# self.  tistos              )

# add the information for Pid efficiency correction
##             self.treatPions   ( tup , b )
##             self.treatKaons   ( tup , b )
##             self.treatMuons   ( tup , b )
##             self.treatTracks  ( tup , b )

# add some reco-summary information
##             self.addRecSummary ( tup , rc_summary   )
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

    # =======================================================================
    # J/psi
    # =======================================================================

    from PhysSelPython.Wrappers import AutomaticData
    jpsi_location = '/Event/AllStreams/Phys/FullDSTDiMuonJpsi2MuMuDetachedLine/Particles'
    jpsi = AutomaticData(Location=jpsi_location)

    from StrippingSelections.StrippingPromptCharm import StrippingPromptCharmConf as PC

    # ======================================
    # REMOVED !!!                         +
    # ( TRGHOSTPROB < 0.5 ) &             +
    # ======================================
    pc = PC('PromptCharm', {
        'TrackCuts'       : """
        ( TRCHI2DOF < 5     ) &
        ( PT > 250 * MeV    ) &
        in_range  ( 2 , ETA , 5 )
        """ ,
        'KaonCuts': ' & ( 0 < PIDK  - PIDpi ) ',
        'PionCuts': ' & ( 0 < PIDpi - PIDK  ) ',
    }
    )

    pions = pc.pions()
    kaons = pc.kaons()

    Preambulo = [
        # shortcut for chi2 of vertex fit
        'chi2vx = VFASPF(VCHI2) ',
        # shortcut for the c*tau
        "from GaudiKernel.PhysicalConstants import c_light",
        # use the embedded cut for chi2(LifetimeFit)<25
        "ctau   = BPVLTIME ( 25 ) * c_light "
    ]

    from GaudiConfUtils.ConfigurableGenerators import CombineParticles
    # ========================================================================
    # prepare Ds+
    # ========================================================================

    Ds_alg = CombineParticles(
        # the decays to be reconstructed
        DecayDescriptor="[D_s+ -> K-  K+ pi+]cc ",
        #
        Preambulo=Preambulo,
        # combination cut : wide mass-cut & PT-cut
        CombinationCut="""
        (   AM12           < 1050 * MeV ) & 
        ( ( ADAMASS('D+')  <   60 * MeV ) | ( ADAMASS('D_s+') <   60 * MeV ) ) & 
        (   APT            >  900 * MeV )
        """ ,
        # mother cut
        MotherCut="""
        (   chi2vx        < 25       ) &
        (   PT            >  1 * GeV ) &
        ( ( ADMASS('D+')  < 50 * MeV ) | ( ADMASS('D_s+') < 50 * MeV ) ) 
        """ ,
        #
        ParticleCombiners={'': 'LoKi::VertexFitter'}
    )

    from PhysSelPython.Wrappers import Selection
    # make selection
    Ds_sel = Selection(
        'DsForPsiC',
        Algorithm=Ds_alg,
        RequiredSelections=[kaons, pions]
    )

    # ========================================================================
    BDs_alg = CombineParticles(
        # the decays to be reconstructed
        DecayDescriptor="[ B_c+ -> J/psi(1S) D_s+ ]cc ",
        #
        Preambulo=Preambulo,
        # combination cut :
        CombinationCut="""
        in_range ( 5 * GeV , AM ,  7 * GeV ) 
        """ ,
        # mother cut
        MotherCut="""
        ( chi2vx  <  16              ) &
        ( ctau    > -10 * micrometer ) 
        """ ,
        #
        ParticleCombiners={'': 'LoKi::VertexFitter'},
        ReFitPVs=True
    )

    # make the selection
    BDs_sel = Selection(
        'PsiDs',
        Algorithm=BDs_alg,
        RequiredSelections=[jpsi, Ds_sel]
    )

    # ========================================================================
    # Bc -> J/psi + 3pi
    # ========================================================================
    from GaudiConfUtils.ConfigurableGenerators import CombineParticles
    bc_3pi = CombineParticles(
        DecayDescriptor='[B_c+ -> J/psi(1S) pi+ pi+ pi-]cc',
        #
        Preambulo=Preambulo,
        DaughtersCuts={
            "J/psi(1S)": " M < 4.0 * GeV "
        },
        #
        CombinationCut="""
        in_range ( 5.9 * GeV , AM ,  6.6 * GeV ) 
        """ ,
        #
        MotherCut="""
        ( chi2vx  < 49              ) &
        ( ctau    > 50 * micrometer ) 
        """ ,
        #
        ParticleCombiners={'': 'LoKi::VertexFitter'},
        ReFitPVs=True
    )
    #
    Bc_3p = Selection(
        'Psi3pi',
        Algorithm=bc_3pi,
        RequiredSelections=[jpsi, pions]
    )

    # ========================================================================
    # Bc -> J/psi + pi
    # ========================================================================
    from GaudiConfUtils.ConfigurableGenerators import CombineParticles
    bc_pi = CombineParticles(
        DecayDescriptor='[B_c+ -> J/psi(1S) pi+ ]cc',
        #
        Preambulo=Preambulo,
        DaughtersCuts={
            "J/psi(1S)": " M < 4.0 * GeV "
        },
        #
        CombinationCut="""
        in_range ( 5.9 * GeV , AM ,  6.6 * GeV ) 
        """ ,
        #
        MotherCut="""
        ( chi2vx  < 16              ) &
        ( ctau    > 50 * micrometer ) 
        """ ,
        ParticleCombiners={'': 'LoKi::VertexFitter'},
        ReFitPVs=True
    )
    #
    Bc_1p = Selection(
        'Psi1Pi',
        Algorithm=bc_pi,
        RequiredSelections=[jpsi, pions]
    )

    from PhysSelPython.Wrappers import SelectionSequence

    Bc_1PI = SelectionSequence("PSIPi", TopSelection=Bc_1p)
    Bc_3PI = SelectionSequence("PSI3Pi", TopSelection=Bc_3p)
    Bc_Ds = SelectionSequence("PSIDs", TopSelection=BDs_sel)

    # Read only fired events to speed up
    from PhysConf.Filters import LoKi_Filters
    fltrs = LoKi_Filters(
        STRIP_Code="HLT_PASS_RE('Stripping.*DiMuonJpsi2MuMuDeta.*')",
        VOID_Code="""
        0.5 < CONTAINS('%s')
        """ % jpsi_location
    )

    # =====================================================================
    from Configurables import DaVinci  # needed for job configuration

    the_year = '2011'

    davinci = DaVinci(
        EventPreFilters=fltrs.filters('FilterMC'),
        DataType=the_year,
        InputType='DST',
        Simulation=True,
        PrintFreq=1000,
        EvtMax=-1,
        #
        HistogramFile='BcMC_Histos.root',
        TupleFile='BcMC.root',
        ## HistogramFile = 'BcMC1pi_Histos.root' ,
        ## TupleFile     = 'BcMC1pi.root'        ,
        #
        # for Ds
        #DDDBtag   = 'MC2011-20120727'          ,
        #CondDBtag = 'MC2011-20120727-vc-mu100' ,
        #
        # for  pi
        # DDDBtag   =  'MC11-20111102' ,
        # CondDBtag =  'sim-20111111-vc-md100' ,
        #
        ## Lumi          = True ,
        #
    )

    from Configurables import GaudiSequencer
    seqDs = GaudiSequencer(
        'DS',
        Members=[Bc_Ds   . sequence(),
                 "MCBc2Ds"]
    )
    from Configurables import GaudiSequencer
    seq1pi = GaudiSequencer(
        'Pi1',
        Members=[Bc_1PI . sequence(),
                 "MCBc1Pi"]
    )

    from Configurables import TrackSmearState
    state_smear = TrackSmearState(
        'StateSmear',
        ## RootInTES  = rootInTES  ,
    )

    davinci.UserAlgorithms = [state_smear, seqDs]

    ## davinci.UserAlgorithms = [ seq1pi ]

    # from Configurables import CondDB
    # CondDB ( LatestGlobalTagByDataType = the_year )

    # ------- decoding set-up start ----------
    # from BenderTools.MicroDST import uDstConf
    # uDstConf ( rootInTES )
    # ------- decoding set-up end  -----------

    # come back to Bender
    setData(datafiles, catalogs, castor)

    #
    # start Gaudi
    #
    gaudi = appMgr()

    alg1 = MCBc(
        'MCBc2Ds',  # Algorithm name ,
        Inputs=[Bc_Ds.outputLocation()],
        PP2MCs=['Relations/Rec/ProtoP/Charged']
    )

    alg2 = MCBc1Pi(
        'MCBc1Pi',  # Algorithm name ,
        Inputs=[Bc_1PI.outputLocation()],
        PP2MCs=['Relations/Rec/ProtoP/Charged']
    )

    return SUCCESS

# =============================================================================
# The actual job steering
if '__main__' == __name__:

    # Bc -> J/psi Ds
    input = [
        '/lhcb/MC/MC11a/ALLSTREAMS.DST/00020276/0000/00020276_00000%03d_1.allstreams.dst' % i for i in range(1, 50)
    ]

    # Bc -> J/psi pi
    input = [
        '/lhcb/MC/MC11a/ALLSTREAMS.DST/00014469/0000/00014469_00000%03d_1.allstreams.dst' % i for i in range(1, 50)
    ]

    configure(input, castor=True)

    run(100)

    import DetCond.HistoCond

    gaudi = appMgr()

    det = gaudi.detSvc()
    cond = det['/dd/Conditions/Calibration/LHCb/MomentumScale']
    cpp.ISolid
    DetDesc = cpp.DetDesc

    print det.ls('/dd/Conditions/Calibration/LHCb')

    #h1      = DetDesc.Params.paramAsHisto2D ( cond , 'IdpPlus'  )
    #h2      = DetDesc.Params.paramAsHisto2D ( cond , 'IdpMinus' )
    #offsets = DetDesc.Params.paramAsHisto1D ( cond , 'Offsets'  )

    #h1      = cond .paramAsHisto2D ( 'IdpPlus'  )
    #h2      = cond .paramAsHisto2D ( 'IdpMinus' )
    #offsets = cond .paramAsHisto1D ( 'Offsets'  )
    #delta   = cond .paramAsDouble  ( 'Delta' )

# =============================================================================
# The END
# =============================================================================
