#!/usr/bin/env ipython
# ========================================================================
# @file DiCharm/BcTMC.py
#
#  Simple algorithm for MC-Bc
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
from Bender.MainMC import *  # import all bender goodies
## easy access to various geometry routines
import LHCbMath.Types
from Gaudi.Configuration import *  # needed for job configuration
# =============================================================================
from GaudiKernel.SystemOfUnits import GeV, MeV, mm, meter
from GaudiKernel.PhysicalConstants import c_light
# =============================================================================
# logging
# =============================================================================
from Bender.Logger import getLogger
logger = getLogger(__name__)
# =============================================================================
import BenderTools.TisTosMC  # add methods for TisTos
import BenderTools.Fill  # add methods for TisTos

# =============================================================================
# @class MCB2PsiH
#  Simple algorithh for B -> J/psi H
#  @date   2011-05-27
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class MCB2PsiH (AlgoMC):

    """
    Simple algorithm 
    """

    def __init__(self, name, **kwargs):

        AlgoMC.__init__(self, name, **kwargs)
        self._decay = '[ Beauty -> ( J/psi(1S) -> mu+ mu- ) (pi+|K+) ]CC'

    def initialize(self):

        sc = AlgoMC.       initialize(self)
        if sc.isFailure():
            return sc

        triggers = {}
        triggers['psi'] = {}
        triggers['psi1'] = {}
        triggers['psi2'] = {}
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

        bees = self.select('b', self._decay)
        if bees.empty():
            print bees
            return self.Warning('No good RC-decays are found!' + self._decay, SUCCESS)

        mcbc = self.mcselect(
            'mcbc', '[ B_c+   => ( J/psi(1S) => mu+ mu- )                  pi+   ]CC')
        mcbp = self.mcselect(
            'mcbp', '[ B+     => ( J/psi(1S) => mu+ mu- )               K+       ]CC')
        mcbd = self.mcselect(
            'mcbp', '[ Beauty => ( J/psi(1S) => mu+ mu- ) ( K*(892)0 => K+ pi- ) ]CC')

        if mcbc.empty() and mcbp.empty() and mcbd.empty():
            self.Warning('No true MC-decays are found!', SUCCESS)
            mcbc = self.mcselect('mcbcall',
                                 (MCABSID == 'B+') |
                                 (MCABSID == 'B0') |
                                 (MCABSID == 'B_c+'))
            print mcbc
            return SUCCESS

        mcpsi = self.mcselect(
            'mcpsi', '[ Beauty => ^( J/psi(1S) => mu+ mu- )  ( pi+ | K+ | K*(892)0 )  ]CC')
        mch = self.mcselect(
            'mch', '[ Beauty =>  ( J/psi(1S) => mu+ mu- ) ^( pi+ | K+ | K*(892)0 )  ]CC')
        mck2 = self.mcselect(
            'mck2', '[ Beauty =>  ( J/psi(1S) => mu+ mu- )  ( K*(892)0 => ^K+  pi- ) ]CC')
        mcp2 = self.mcselect(
            'mcp2', '[ Beauty =>  ( J/psi(1S) => mu+ mu- )  ( K*(892)0 =>  K+ ^pi- ) ]CC')

        mcBc = NONE if mcbc .empty() else MCTRUTH(self.mcTruth(), mcbc)
        mcBp = NONE if mcbp .empty() else MCTRUTH(self.mcTruth(), mcbp)
        mcBd = NONE if mcbd .empty() else MCTRUTH(self.mcTruth(), mcbd)
        mcPsi = NONE if mcpsi.empty() else MCTRUTH(self.mcTruth(), mcpsi)
        mcH = NONE if mch  .empty() else MCTRUTH(self.mcTruth(), mch)
        mcK2 = NONE if mck2 .empty() else MCTRUTH(self.mcTruth(), mck2)
        mcP2 = NONE if mcp2 .empty() else MCTRUTH(self.mcTruth(), mcp2)

        tup = self.nTuple('B')
        for b in bees:

            m = M(b) / GeV
            if m > 7.1:
                continue
            #
            psi = b(1)
            h = b(2)
            h2 = b(3)

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

            _mcbc = True if mcBc(b) else False
            _mcbp = True if mcBp(b) else False
            _mcbd = True if mcBd(b) else False
            _mcpsi = True if mcPsi(psi) else False

            _mch = True if mcH(h) else False

            tup.column_bool('mcbc', _mcbc)
            tup.column_bool('mcbp', _mcbp)
            tup.column_bool('mcbd', _mcbd)
            tup.column_bool('mcpsi', _mcpsi)
            tup.column_bool('mch', _mch)

            tup.column_float('mass', self.f_mfit(b))
            tup.column_float('c2dtf', c2)
            tup.column_float('ctau', ct)

            mcct_c = -1 * meter
            mcct_p = -1 * meter
            mcct_0 = -1 * meter
            if not mcbc.empty():
                mcct_c = MCCTAU(mcbc(0))
            if not mcbp.empty():
                mcct_p = MCCTAU(mcbp(0))
            if not mcbd.empty():
                mcct_0 = MCCTAU(mcbd(0))

            tup.column_float('mcct_c', mcct_c)
            tup.column_float('mcct_p', mcct_p)
            tup.column_float('mcct_0', mcct_0)

            self.treatKine(tup, b, '_b')
            self.treatKine(tup, psi, '_psi')
            self.treatKine(tup, h, '_h')

            if h2:

                _mch2 = True if mcH(h2) else False
                _mck2 = True if mcK2(h) else False
                _mcp2 = True if mcP2(h2) else False
                tup.column_bool('mch2', _mch)
                tup.column_bool('mck2', _mck2)
                tup.column_bool('mcp2', _mcp2)
                self.treatKine(tup, h2, '_h2')

            # add the information needed for TisTos
            self.tisTos(psi, tup, 'psi_',
                        self.lines['psi'], self.l0tistos, self.tistos)
            self.tisTos(psi, tup, 'psi1_',
                        self.lines['psi1'], self.l0tistos, self.tistos)
            self.tisTos(psi, tup, 'psi2_',
                        self.lines['psi2'], self.l0tistos, self.tistos)

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

    #
    # pseudo stripping & WG-production
    #
    # from PhysSelPython.Wrappers import AutomaticData
    # jpsi_location = '/Event/AllStreams/Phys/FullDSTDiMuonJpsi2MuMuDetachedLine/Particles'
    # jpsi = AutomaticData ( Location = jpsi_location )

# if   '2012' == the_year :
##         import StrippingArchive.Stripping20.StrippingDiMuonNew              as DiMuon
##         import StrippingSettings.Stripping20.LineConfigDictionaries_BandQ   as LineSettings
# elif '2011' == the_year :
##         import StrippingArchive.Stripping20r1.StrippingDiMuonNew            as DiMuon
##         import StrippingSettings.Stripping20r1.LineConfigDictionaries_BandQ as LineSettings
##     config  = LineSettings.FullDSTDiMuon['CONFIG']
##     name    = 'FullDST'
##     builder = DiMuon.DiMuonConf ( name , config )
# selection
##     jpsi  = builder.SelJpsi2MuMuDetached
    import CommonParticles.StdLooseJpsi2MuMu

    # from StandardParticles import StdLooseJpsi2MuMu
    # jpsi  = StdLooseJpsi2MuMu

    from PhysSelPython.Wrappers import AutomaticData
    jpsi = AutomaticData('/Event/Phys/StdLooseJpsi2MuMu/Particles')

    from GaudiConfUtils.ConfigurableGenerators import FilterDesktop
    # pions :
    alg_pions = FilterDesktop(
        #
        Code="""
        ( PT > 500  * MeV      ) & 
        ( CLONEDIST   > 5000   ) & 
        ( TRGHOSTPROB < 0.5    ) &
        ( TRCHI2DOF   < 4      ) & 
        in_range ( 2 , ETA , 5 ) &
        HASRICH                  &
        ( PIDpi - PIDK > -2    ) 
        """ ,
    )

    from PhysSelPython.Wrappers import Selection
    from StandardParticles import StdAllLoosePions
    pions = Selection(
        "ThePions",
        Algorithm=alg_pions,
        RequiredSelections=[StdAllLoosePions]
    )

    alg_kaons = FilterDesktop(
        #
        Code="""
        ( PT > 500  * MeV      ) & 
        ( CLONEDIST   > 5000   ) & 
        ( TRCHI2DOF   < 4      ) & 
        ( TRGHOSTPROB < 0.5    ) & 
        in_range ( 2 , ETA , 5 ) &
        HASRICH                  &
        ( PIDK - PIDpi > -2    ) 
        """
    )

    from StandardParticles import StdAllLooseKaons
    kaons = Selection(
        "TheKaons",
        Algorithm=alg_kaons,
        RequiredSelections=[StdAllLooseKaons]
    )

    #
    Preambulo = [
        # shortcut for chi2 of vertex fit
        "chi2vx   = VFASPF(VCHI2)",
        # shortcut for the c*tau
        "from GaudiKernel.PhysicalConstants import c_light",
        # use the embedded cut for chi2(LifetimeFit)<16
        "ctau_25   = BPVLTIME ( 25 ) * c_light ",
        "mbc_acut  = in_range ( 6.050 * GeV , AM , 6.550 * GeV ) ",
        "mbp_acut  = in_range ( 5.100 * GeV , AM , 5.550 * GeV ) ",
        # mass-cut for beauty particles
        "mbc_cut   = in_range ( 6.100 * GeV ,  M , 6.500 * GeV ) ",
        "mbp_cut   = in_range ( 5.150 * GeV ,  M , 5.500 * GeV ) ",
    ]

    # Bc :
    from GaudiConfUtils.ConfigurableGenerators import CombineParticles
    alg_bc = CombineParticles(
        DecayDescriptor="[B_c+ -> J/psi(1S) pi+ ]cc",
        Preambulo=Preambulo,
        CombinationCut=" mbc_acut ",
        MotherCut="""
        ( chi2vx < 20  ) & mbc_cut &  ( ctau_25 > %s )
        """ %  ( 40 * micrometer )
    )

    sel_bc = Selection(
        "TheBc",
        Algorithm=alg_bc,
        RequiredSelections=[jpsi, pions]
    )

    # B+ :
    alg_bu = CombineParticles(
        DecayDescriptor="[B+ -> J/psi(1S) K+ ]cc",
        Preambulo=Preambulo,
        CombinationCut=" mbp_acut ",
        MotherCut="""
        ( chi2vx < 20  ) & mbp_cut &  ( ctau_25 > %s )
        """ %  ( 40 * micrometer )
    )

    sel_bu = Selection(
        "TheB",
        Algorithm=alg_bu,
        RequiredSelections=[jpsi, kaons]
    )

    # B0 :
    alg_bd = CombineParticles(
        DecayDescriptor="[B0 -> J/psi(1S) K+ pi-]cc",
        Preambulo=Preambulo,
        CombinationCut="""
        mbp_acut & in_range ( 600 * MeV , AM23  , 1.2 * GeV ) 
        """ ,
        MotherCut="""
        ( chi2vx < 50  ) & mbp_cut &  ( ctau_25 > %s )
        """ %  ( 40 * micrometer )
    )

    sel_bd = Selection(
        "TheB0",
        Algorithm=alg_bd,
        RequiredSelections=[jpsi, kaons, pions]
    )

    from PhysSelPython.Wrappers import SelectionSequence
    selseq_pi = SelectionSequence("Bc", TopSelection=sel_bc)
    selseq_k = SelectionSequence("Bu", TopSelection=sel_bu)
    selseq_kst = SelectionSequence("Bd", TopSelection=sel_bd)

    # Read only fired events to speed up
    from PhysConf.Filters import LoKi_Filters
    fltrs = LoKi_Filters(
        #
        # Require at least one primary vertex
        #
        VOID_Code="""
        ( RECSUMMARY (  0 , -1 ) >  0.5 ) 
        """
    )

    #
    # make the final sequencers:
    #

    from Configurables import GaudiSequencer
    seq_pi = GaudiSequencer(
        'PION',
        Members=[selseq_pi  . sequence(), "MCB2PI"]
    )
    seq_k = GaudiSequencer(
        'KAON',
        Members=[selseq_k   . sequence(), "MCB2K"]
    )
    seq_kst = GaudiSequencer(
        'KSTAR',
        Members=[selseq_kst . sequence(), "MCB2KST"]
    )

    #
    # make the final choice
    #

    mode = params['Mode']
    mode = mode.upper()
    if 0 <= mode.find('BC') or 0 <= mode.find('PI'):
        seq = seq_pi
    elif  0 <= mode.find( 'B+' ) or 0 <= mode.find('BP') or \
            0 <= mode.find('BU') or 0 <= mode.find('K+'):
        seq = seq_k
    elif 0 <= mode.find('B0') or 0 <= mode.find('BZ') or 0 <= mode.find('K*'):
        seq = seq_kst

    #
    # finally: DaVinci
    #
    ## needed for job configuration
    from Configurables import DaVinci
    davinci = DaVinci(
        EventPreFilters=fltrs.filters('Filters'),
        DataType=the_year,
        InputType='MDST',
        Simulation=True,
        PrintFreq=1000,
        EvtMax=-1,
        #
        HistogramFile='MCBc_Histos.root',
        TupleFile='MCBc.root',
        #
        DDDBtag=params['DDDB'],
        CondDBtag=params['SIMCOND'],
        #
        Lumi=False  # True ,
        #
    )

    from BenderTools.Utils import silence
    silence()

    #
    # finally inform Davinci about algorithsm
    #
    davinci.UserAlgorithms = [
        seq
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

    #
    alg1 = MCB2PsiH(
        'MCB2PI',  # Algorithm name ,
        Inputs=[selseq_pi.outputLocation()],
        PP2MCs=['Relations/Rec/ProtoP/Charged'],
        ReFitPVs=True
    )

    alg2 = MCB2PsiH(
        'MCB2K',  # Algorithm name ,
        Inputs=[selseq_k.outputLocation()],
        PP2MCs=['Relations/Rec/ProtoP/Charged'],
        ReFitPVs=True
    )

    alg3 = MCB2PsiH(
        'MCB2KST',  # Algorithm name ,
        Inputs=[selseq_kst.outputLocation()],
        PP2MCs=['Relations/Rec/ProtoP/Charged'],
        ReFitPVs=True
    )

    alg3._decay = '[ Beauty -> ( J/psi(1S) -> mu+ mu- ) K+ pi- ]CC'

    return SUCCESS

# =============================================================================
# The actual job steering
if '__main__' == __name__:

    #
    # Bc->J/psi pi
    #
    # path = '/MC/MC11a/Beam3500GeV-2011-MagDown-Nu2-EmNoCuts/Sim05/Trig0x40760037Flagged/Reco12a/Stripping17NoPrescalingFlagged/42122001/ALLSTREAMS.DST'
    inputs_Bc = [
        '/lhcb/MC/MC11a/ALLSTREAMS.DST/00014472/0000/00014472_00000006_1.allstreams.dst',
        '/lhcb/MC/MC11a/ALLSTREAMS.DST/00014472/0000/00014472_00000015_1.allstreams.dst',
        '/lhcb/MC/MC11a/ALLSTREAMS.DST/00014472/0000/00014472_00000053_1.allstreams.dst',
        '/lhcb/MC/MC11a/ALLSTREAMS.DST/00014472/0000/00014472_00000011_1.allstreams.dst'
    ]
    params_Bc = {
        'Mode': 'Bc -> J/psi pi',
        'Year': '2011',
        'DDDB': 'MC11-20111102',
        'SIMCOND': 'sim-20111111-vc-md100'
    }

    #
    # B+ -> J/psi K
    #
    # path = '/MC/MC11a/Beam3500GeV-2011-MagDown-Nu2-EmNoCuts/Sim05/Trig0x40760037Flagged/Reco12/Stripping17NoPrescalingFlagged/12143001/ALLSTREAMS.DST'
    inputs_Bu = [
        '/lhcb/MC/MC11a/ALLSTREAMS.DST/00013192/0000/00013192_00000001_1.allstreams.dst',
        '/lhcb/MC/MC11a/ALLSTREAMS.DST/00013416/0000/00013416_00000022_1.allstreams.dst',
        '/lhcb/MC/MC11a/ALLSTREAMS.DST/00013192/0000/00013192_00000010_1.allstreams.dst',
        '/lhcb/MC/MC11a/ALLSTREAMS.DST/00013416/0000/00013416_00000017_1.allstreams.dst',
        '/lhcb/MC/MC11a/ALLSTREAMS.DST/00013416/0000/00013416_00000031_1.allstreams.dst'
    ]
    params_Bu = {
        'Mode': 'B+ -> J/psi K+',
        'Year': '2011',
        'DDDB': 'MC11-20111102',
        'SIMCOND': 'sim-20111111-vc-md100'
    }
    #

    #
    # B0 -> J/psi K*
    #
    # path = '/MC/MC11a/Beam3500GeV-2011-MagDown-Nu2-EmNoCuts/Sim05/Trig0x40760037Flagged/Reco12a/Stripping17NoPrescalingFlagged/11144001/ALLSTREAMS.DST'
    inputs_Bd = [
        '/lhcb/MC/MC11a/ALLSTREAMS.DST/00016191/0000/00016191_00000176_1.allstreams.dst',
        '/lhcb/MC/MC11a/ALLSTREAMS.DST/00016191/0000/00016191_00000250_1.allstreams.dst',
        '/lhcb/MC/MC11a/ALLSTREAMS.DST/00015672/0000/00015672_00000066_1.allstreams.dst',
        '/lhcb/MC/MC11a/ALLSTREAMS.DST/00015672/0000/00015672_00000037_1.allstreams.dst',
        '/lhcb/MC/MC11a/ALLSTREAMS.DST/00016191/0000/00016191_00000296_1.allstreams.dst',
        '/lhcb/MC/MC11a/ALLSTREAMS.DST/00016191/0000/00016191_00000115_1.allstreams.dst',
        '/lhcb/MC/MC11a/ALLSTREAMS.DST/00015672/0000/00015672_00000073_1.allstreams.dst'
    ]
    params_Bd = {
        'Mode': 'B0 -> J/psi K*',
        'Year': '2011',
        'DDDB': 'MC11-20111102',
        'SIMCOND': 'sim-20111111-vc-md100'
    }

    configure(
        inputs_Bd,
        castor=True,
        params=params_Bd
    )

    run(1000)

    gaudi = appMgr()


# =============================================================================
# The END
# =============================================================================
