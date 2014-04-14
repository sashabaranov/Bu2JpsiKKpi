#!/usr/bin/env gaudirun.py
# =============================================================================
# $Id:$
# =============================================================================
# @file WriteCharmY.py
#
# Helper script to write muDST with Charm+Y
#
# @author Vanya BELYAEV Ivan.Belyaev@itep.ru
# @date 2013-09-01
#
#                   $Revision:$
# Last modification $Date$
#                by $Author$
# =============================================================================
"""
Helper script to write MicroDST with Charm+Y selection 
"""
# =============================================================================
__author__ = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__ = "2013-09-01"
__version__ = "$Revision:$"
# =============================================================================
from GaudiKernel.SystemOfUnits import MeV, GeV
# =============================================================================
# logging
# =============================================================================
try:
    from AnalysisPython.Logger import getLogger
    _name = __name__
    if 0 == _name.find('__'):
        _name = 'C&Y'
    logger = getLogger(_name)
except:
    from Gaudi.Configuration import log as logger
# =============================================================================


# =============================================================================
# 0. get upsilons from DIMUON.DST
# =============================================================================

ups_name = 'FullDSTDiMuonDiMuonHighMass'

from PhysSelPython.Wrappers import AutomaticData
ups = AutomaticData('/Event/Dimuon/Phys/%sLine/Particles' % ups_name)

# =============================================================================
# 1. filter them to be in agreement with BOTTOM.MDST
# =============================================================================
# Upsilon -> mumu, cuts by  Giulia Manca
# =============================================================================
from GaudiConfUtils.ConfigurableGenerators import FilterDesktop
alg_ups = FilterDesktop(
    #
    Code="""
    ( M > 7 * GeV ) & 
    DECTREE ( 'Meson -> mu+ mu-'   )                      &
    CHILDCUT( 1 , HASMUON & ISMUON )                      &
    CHILDCUT( 2 , HASMUON & ISMUON )                      & 
    ( MINTREE ( 'mu+' == ABSID , PT ) > 1 * GeV         ) &
    ( MAXTREE ( ISBASIC & HASTRACK , TRCHI2DOF ) < 4    ) & 
    ( MINTREE ( ISBASIC & HASTRACK , CLONEDIST ) > 5000 ) & 
    ( VFASPF  ( VPCHI2 ) > 0.5/100 ) 
    & ( abs ( BPV ( VZ    ) ) <  0.5 * meter     ) 
    & (       BPV ( vrho2 )   < ( 10 * mm ) ** 2 ) 
    """ ,
    #
    Preambulo=[
        "vrho2 = VX**2 + VY**2"
    ],
    #
    ReFitPVs=True
    #
)

from PhysSelPython.Wrappers import Selection
sel_ups = Selection(
    'Upsilon',
    Algorithm=alg_ups,
    RequiredSelections=[ups]
)


# =============================================================================
# 2. get prompt charm selection
# =============================================================================
from StrippingSelections.StrippingPromptCharm import StrippingPromptCharmConf as PC
from StrippingSelections.StrippingPromptCharm import FITTER

from GaudiConfUtils.ConfigurableGenerators import FilterDesktop
from GaudiConfUtils.ConfigurableGenerators import CombineParticles

# =============================================================================
# unify the pion& kaon  selections
# =============================================================================
_PionCut_  = """
( CLONEDIST   > 5000   ) & 
( TRCHI2DOF   < 4      ) &
( TRGHOSTPROB < 0.4    ) &
( PT          > 250 * MeV               ) & 
in_range ( 2          , ETA , 4.9       ) &
in_range ( 3.2 * GeV  , P   , 150 * GeV ) &
HASRICH                  &
( PROBNNpi     > 0.15  ) &
( MIPCHI2DV()  > 9.    )
"""
_KaonCut_  = """
( CLONEDIST   > 5000   ) & 
( TRCHI2DOF   < 4      ) & 
( TRGHOSTPROB < 0.4    ) & 
( PT          > 200 * MeV               ) & 
in_range ( 2          , ETA , 4.9       ) &
in_range ( 3.2 * GeV  , P   , 150 * GeV ) &
HASRICH                  &
( PROBNNk      > 0.15  ) &
( MIPCHI2DV()  > 9.    ) 
"""
_ProtonCut_ = """
( CLONEDIST   > 5000    ) & 
( TRCHI2DOF   < 4       ) & 
( TRGHOSTPROB < 0.4     ) & 
( PT          > 200 * MeV              ) & 
in_range (  2        , ETA , 4.9       ) &
in_range ( 10 * GeV  , P   , 150 * GeV ) &
HASRICH                   &
( PROBNNp      > 0.15   ) &
( MIPCHI2DV()  > 9.     ) 
"""

from GaudiConfUtils.ConfigurableGenerators import FilterDesktop
_alg_pi = FilterDesktop(
    #
    Code=_PionCut_,
    #
)

from PhysSelPython.Wrappers import Selection
from StandardParticles import StdAllNoPIDsPions as input_pions
pions = Selection(
    "SelPiForCEW",
    Algorithm=_alg_pi,
    RequiredSelections=[input_pions]
)

from GaudiConfUtils.ConfigurableGenerators import FilterDesktop
_alg_k = FilterDesktop(
    #
    Code=_KaonCut_,
    #
)

from PhysSelPython.Wrappers import Selection
from StandardParticles import StdAllNoPIDsKaons as input_kaons
kaons = Selection(
    "SelKForCEW",
    Algorithm=_alg_k,
    RequiredSelections=[input_kaons]
)

from GaudiConfUtils.ConfigurableGenerators import FilterDesktop
_alg_p = FilterDesktop(
    #
    Code=_ProtonCut_,
    #
)

from PhysSelPython.Wrappers import Selection
from StandardParticles import StdAllNoPIDsProtons as input_protons
protons = Selection(
    "SelPForCEW",
    Algorithm=_alg_p,
    RequiredSelections=[input_protons]
)


def _protons_(self):
    return protons


def _kaons_(self):
    return kaons


def _pions_(self):
    return pions

#
# get the selections
#
PC.pions = _pions_
PC.kaons = _kaons_
PC.protons = _protons_
#
logger.warning("Redefine PromptCharm.pions   ")
logger.warning("Redefine PromptCharm.kaons   ")
logger.warning("Redefine PromptCharm.protons ")
#

# =============================================================================
# Lambda_C* -> Lambda_C pi pi selection
# =============================================================================


def _LamCstar_(self):

    sel = self._selection('LambdaCstarForPromptCharm_Selection')
    if sel:
        return sel

    cmb = CombineParticles(
        #
        Monitor=self['Monitor'],
        HistoProduce=self['Monitor'],
        #
        DecayDescriptors=[
            " [ Lambda_c(2625)+ -> Lambda_c+ pi+ pi-]cc"
        ],
        #
        DaughtersCuts={
            'Lambda_c+': ' INTREE ( good_proton ) ',
            'pi+':  self.slowPionCuts()
        },
        #
        Preambulo=self.preambulo() + [
            "hDmLamC = Gaudi.Histo1DDef ( 'dm(Lambda_c)' , 200 , 700 , 200 )"
        ],
        #
        CombinationCut="""
        AM - AM1 < 650 * MeV
        """ ,
        MotherCut="""
        ( chi2vx  < 25 )
        """ ,
        #
        MotherMonitor="""
        process ( monitor ( M-M1 , hDmLamC , 'dm(Lambda_c)' ) ) >> ~EMPTY
        """ ,
        #
        # make the selection faster
        #
        ParticleCombiners={'': FITTER},
        #
    )

    from StandardParticles import StdAllNoPIDsPions
    # convert it to selection
    sel = Selection(
        "SelLambdaCstarFor" + self.name(),
        Algorithm=cmb,
        RequiredSelections=[self.LamC(),
                            StdAllNoPIDsPions]  # slow prompt pion!
    )

    return self._add_selection('LambdaCstarForPromptCharm_Selection',  sel)

PC.LamCstar = _LamCstar_
logger.warning("Redefine PromptCharm.LamcStar ")

# ========================================================================
pc = PC('PromptCharm', {
    #
    'pT(D0)':  0.50 * GeV,  # pt-cut for  prompt   D0
    'pT(D0) for D*+':  0.50 * GeV,  # pt-cut for  D0 from  D*+
    'pT(D+)':  0.50 * GeV,  # pt-cut for  prompt   D+
    'pT(Ds+)':  0.50 * GeV,  # pt-cut for  prompt   Ds+
    'pT(Lc+)':  0.50 * GeV,  # pt-cut for  prompt   Lc+
    #
})

# =============================================================================
# 3. combine upsilon + charms
# =============================================================================
from GaudiConfUtils.ConfigurableGenerators import CombineParticles
# =============================================================================
alg_D0 = CombineParticles(
    #
    DecayDescriptor="[ chi_b1(1P) -> J/psi(1S) D0 ]cc ",
    #
    CombinationCut=" AALL ",
    MotherCut=" PALL",
    #
    ParticleCombiners={'': 'LoKi::VertexFitter:PUBLIC'},
    #
    ReFitPVs=True
    #
)
# =============================================================================
alg_Dp = CombineParticles(
    #
    DecayDescriptor="[ chi_b1(1P) -> J/psi(1S) D+ ]cc ",
    #
    CombinationCut=" AALL ",
    MotherCut=" PALL",
    #
    ParticleCombiners={'': 'LoKi::VertexFitter:PUBLIC'},
    #
    ReFitPVs=True
    #
)
# =============================================================================
alg_Ds = CombineParticles(
    #
    DecayDescriptor="[ chi_b1(1P) -> J/psi(1S) D_s+ ]cc ",
    #
    CombinationCut=" AALL ",
    MotherCut=" PALL",
    #
    ParticleCombiners={'': 'LoKi::VertexFitter:PUBLIC'},
    #
    ReFitPVs=True
    #
)
# =============================================================================
alg_Dst = CombineParticles(
    #
    DecayDescriptor="[ chi_b1(1P) -> J/psi(1S) D*(2010)+ ]cc ",
    #
    CombinationCut=" AALL ",
    MotherCut=" PALL",
    #
    ParticleCombiners={'': 'LoKi::VertexFitter:PUBLIC'},
    #
    ReFitPVs=True
    #
)
# ============================================================================
alg_Lc = CombineParticles(
    #
    DecayDescriptor="[ chi_b1(1P) -> J/psi(1S) Lambda_c+ ]cc ",
    #
    CombinationCut=" AALL ",
    MotherCut=" PALL",
    #
    ParticleCombiners={'': 'LoKi::VertexFitter:PUBLIC'},
    #
    ReFitPVs=True
    #
)
# ============================================================================
alg_Lcst = CombineParticles(
    #
    DecayDescriptor="[ chi_b1(1P) -> J/psi(1S) Lambda_c(2625)+ ]cc ",
    #
    CombinationCut=" AALL ",
    MotherCut=" PALL",
    #
    ParticleCombiners={'': 'LoKi::VertexFitter:PUBLIC'},
    #
    ReFitPVs=True
    #
)
# ============================================================================
alg_Sc = CombineParticles(
    #
    DecayDescriptors=[
        "[ chi_b1(1P) -> J/psi(1S) Sigma_c0  ]cc ",
        "[ chi_b1(1P) -> J/psi(1S) Sigma_c++ ]cc "
    ],
    #
    CombinationCut=" AALL ",
    MotherCut=" PALL",
    #
    ParticleCombiners={'': 'LoKi::VertexFitter:PUBLIC'},
    #
    ReFitPVs=True
    #
)

from PhysSelPython.Wrappers import Selection
sel_D0 = Selection(
    'Y&D0',
    Algorithm=alg_D0,
    RequiredSelections=[sel_ups, pc.D02HH()])
sel_Dp = Selection(
    'Y&D+',
    Algorithm=alg_Dp,
    RequiredSelections=[sel_ups, pc.Dplus()])
sel_Ds = Selection(
    'Y&Ds+',
    Algorithm=alg_Ds,
    RequiredSelections=[sel_ups, pc.Ds()])
sel_Lc = Selection(
    'Y&Lc+',
    Algorithm=alg_Lc,
    RequiredSelections=[sel_ups, pc.LamC()])
sel_Dst = Selection(
    'Y&Dstar+',
    Algorithm=alg_Dst,
    RequiredSelections=[sel_ups, pc.Dstar()])
sel_Lcst = Selection(
    'Y&Lcstar+',
    Algorithm=alg_Lcst,
    RequiredSelections=[sel_ups, pc.LamCstar()])
sel_Sc = Selection(
    'Y&Sigmac',
    Algorithm=alg_Sc,
    RequiredSelections=[sel_ups, pc.SigC()])

# =============================================================================
# 4. make selection sequences and build the stream
# =============================================================================
from PhysSelPython.Wrappers import SelectionSequence, MultiSelectionSequence

ups_charm = MultiSelectionSequence(
    #
    'YC', Sequences=[
        #
        SelectionSequence('SEQ:Y&D0',  TopSelection=sel_D0),
        SelectionSequence('SEQ:Y&D+',  TopSelection=sel_Dp),
        SelectionSequence('SEQ:Y&Ds+',  TopSelection=sel_Ds),
        SelectionSequence('SEQ:Y&Lc+',  TopSelection=sel_Lc),
        SelectionSequence('SEQ:Y&Dstar+',  TopSelection=sel_Dst),
        SelectionSequence('SEQ:Y&Lcstar',  TopSelection=sel_Lcst),
        SelectionSequence('SEQ:Y&Sigmac',  TopSelection=sel_Sc),
        #
    ])

# =============================================================================
# 5. micro-DST writer
# =============================================================================
from DSTWriters.Configuration import (SelDSTWriter,
                                      stripMicroDSTStreamConf,
                                      stripMicroDSTElements)

# Configuration of SelDSTWriter
SelDSTWriterConf = {'default': stripMicroDSTStreamConf(pack=False)}
SelDSTWriterElements = {'default': stripMicroDSTElements(pack=False)}

uDstWriter = SelDSTWriter(
    "MyDSTWriter",
    StreamConf=SelDSTWriterConf,
    MicroDSTElements=SelDSTWriterElements,
    OutputFileSuffix='BandQ',  # output PRE-fix!
    SelectionSequences=[ups_charm]
)

#
# Read only fired events to speed up
#
from PhysConf.Filters import LoKi_Filters
fltrs = LoKi_Filters(
    STRIP_Code="""
    HLT_PASS_RE ( 'Stripping.*%s.*' )
    """ % ups_name
)

# ========================================================================
# The last step: DaVinci & DB
# ========================================================================

##the_year = '2012'
the_year = '2011'

from Configurables import DaVinci
davinci = DaVinci(
    #
    DataType=the_year,  # ATTENTION !!
    #
    EventPreFilters=fltrs.filters('Filters'),
    InputType="DST",
    EvtMax=-1,
    PrintFreq=10000,
    Lumi=True,
)

davinci.appendToMainSequence([uDstWriter.sequence()])

from Configurables import CondDB
CondDB(LatestGlobalTagByDataType=the_year)  # ATTENTION !!


# TEST data sample

# input = [
# '/lhcb/LHCb/Collision12/DIMUON.DST/00020198/0001/00020198_00012742_1.dimuon.dst',
# '/lhcb/LHCb/Collision12/DIMUON.DST/00020198/0001/00020198_00018390_1.dimuon.dst',
# '/lhcb/LHCb/Collision12/DIMUON.DST/00020198/0001/00020198_00016364_1.dimuon.dst',
# '/lhcb/LHCb/Collision12/DIMUON.DST/00020198/0001/00020198_00015767_1.dimuon.dst',
# '/lhcb/LHCb/Collision12/DIMUON.DST/00020198/0002/00020198_00020410_1.dimuon.dst',
# '/lhcb/LHCb/Collision12/DIMUON.DST/00020198/0000/00020198_00007306_1.dimuon.dst',
# '/lhcb/LHCb/Collision12/DIMUON.DST/00020198/0001/00020198_00014550_1.dimuon.dst',
# '/lhcb/LHCb/Collision12/DIMUON.DST/00020198/0001/00020198_00016402_1.dimuon.dst',
# '/lhcb/LHCb/Collision12/DIMUON.DST/00020198/0001/00020198_00014253_1.dimuon.dst',
# '/lhcb/LHCb/Collision12/DIMUON.DST/00020198/0001/00020198_00014185_1.dimuon.dst'
# ]
# Input data
## from GaudiConf import IOHelper
## PFN  = 'PFN:root://eoslhcb.cern.ch//eos'
## ioh = IOHelper()
## ioh.inputFiles ( [ PFN + i for i in input   ] )

## davinci.EvtMax     = 100000
## davinci.PrintFreq  =    100
## davinci.EvtMax     =  50000


# =============================================================================
# The END
# =============================================================================
