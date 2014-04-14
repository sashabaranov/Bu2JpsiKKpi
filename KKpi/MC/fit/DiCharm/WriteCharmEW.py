#!/usr/bin/env gaudirun.py
# =============================================================================
# $Id:$
# =============================================================================
# @file WriteCharmEW.py
#
# Helper script to write muDST with Charm+W/Z0
#
# @author Vanya BELYAEV Ivan.Belyaev@itep.ru
# @date 2013-09-01
#
#                   $Revision:$
# Last modification $Date$
#                by $Author$
# =============================================================================
"""
Helper script to write MicroDST with Charm+W/Z0 selection 
"""
# =============================================================================
__author__ = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__ = "2013-09-01"
__version__ = "$Revision:$"
# =============================================================================
# logging
# =============================================================================
# logging
# =============================================================================
try:
    from AnalysisPython.Logger import getLogger
    _name = __name__
    if 0 == _name.find('__'):
        _name = 'C&W/Z'
    logger = getLogger(_name)
except:
    from Gaudi.Configuration import log as logger
# =============================================================================

from GaudiKernel.SystemOfUnits import MeV, GeV

from PhysSelPython.Wrappers import AutomaticData, Selection, SelectionSequence

from StandardParticles import StdLooseAllPhotons  # needed for chi_b
from GaudiConfUtils.ConfigurableGenerators import FilterDesktop

#
## W&Z-locations  in EW.DST
#
W_Location = '/Event/EW/Phys/WMuLine/Particles'
Z_Location = '/Event/EW/Phys/Z02MuMuLine/Particles'
from PhysSelPython.Wrappers import AutomaticData
W_Strip = AutomaticData(Location=W_Location)
Z_Strip = AutomaticData(Location=Z_Location)

#
# get the prompt charm
#
from StrippingSelections.StrippingPromptCharm import StrippingPromptCharmConf as STRIP
#
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
STRIP.pions = _pions_
STRIP.kaons = _kaons_
STRIP.protons = _protons_
#
logger.warning("Redefine PromptCharm.pions   ")
logger.warning("Redefine PromptCharm.kaons   ")
logger.warning("Redefine PromptCharm.protons ")
#
strip = STRIP('PromptCharm', {
    'pT(D0)':  0.5 * GeV,  # pt-cut for  prompt   D0
    'pT(D0) for D*+':  0.5 * GeV,  # pt-cut for  D0 from  D*+
    'pT(D+)':  0.5 * GeV,  # pt-cut for  prompt   D+
    'pT(Ds+)':  0.5 * GeV,  # pt-cut for  prompt   Ds+
    'pT(Lc+)':  0.5 * GeV,  # pt-cut for  prompt   Lc+
    #
})
#
# the prompt charm selection
#
charm = strip.PromptCharm()
dimu = strip.DiMuon()

# =============================================================================
# Common preambulo
# =============================================================================
EW_preambulo = [
    "pion_cuts  = in_range ( 300 * MeV , PT , 10 * GeV ) & ( CLONEDIST > 5000 ) & ( TRCHI2DOF < 5 ) & ( TRGHOSTPROB < 0.5 ) & ( PERR2/P2 < 0.05**2 ) ",
    "ptCone_    =  SUMCONE (   0.25 , PT , '/Event/Phys/StdAllLoosePions/Particles'               )",
    "ptCone_2   =  SUMCONE (   0.25 , PT , '/Event/Phys/StdAllLoosePions/Particles'   , pion_cuts )",
    "etCone_    =  SUMCONE (   0.25 , PT , '/Event/Phys/StdLooseAllPhotons/Particles'             )",
    "ptCone     =    SINFO (  55001 , ptCone_  , True ) ",
    "ptCone2    =    SINFO (  55003 , ptCone_2 , True ) ",
    "etCone     =    SINFO (  55002 , etCone_  , True ) "
]
# ========================================================================
# good W
# ========================================================================
from GaudiConfUtils.ConfigurableGenerators import FilterDesktop
gW = FilterDesktop(
    Preambulo=EW_preambulo,
    Code="""
    in_range ( 15 * GeV , PT , 100 * GeV ) &
    ( -1e+10 * GeV < ptCone  ) &
    ( -1e+10 * GeV < ptCone2 ) &
    ( -1e+10 * GeV < etCone  ) 
    """
)
W_Data = Selection(
    'W',
    Algorithm=gW,
    RequiredSelections=[W_Strip]
)

# ========================================================================
## charm + W
# ========================================================================
from GaudiConfUtils.ConfigurableGenerators import CombineParticles
cW_ = CombineParticles(
    DecayDescriptors=[
        "[H_10 -> D0        mu+ ]cc",
        "[H_10 -> D0        mu- ]cc",
        "[H_10 -> D+        mu+ ]cc",
        "[H_10 -> D+        mu- ]cc",
        "[H_10 -> D_s+      mu+ ]cc",
        "[H_10 -> D_s+      mu- ]cc",
        "[H_10 -> Lambda_c+ mu+ ]cc",
        "[H_10 -> Lambda_c+ mu- ]cc",
    ],
    #
    CombinationCut=" AALL ",
    MotherCut=" PALL",
    ReFitPVs=True
)
#
cW = Selection(
    'CharmAndW',
    Algorithm=cW_,
    RequiredSelections=[W_Data, charm]
)
# =============================================================================
## dimu + W
# ========================================================================
mW_ = CombineParticles(
    #
    DecayDescriptors=[
        "[H_10 -> J/psi(1S) mu+ ]cc",
    ],
    #
    DaughtersCuts={
        'mu+': 'BPVIPCHI2() < 20 ',
        'J/psi(1S)': 'M < 3.2 * GeV '
    },
    CombinationCut=" AALL ",
    MotherCut=" VFASPF ( VCHI2PDOF ) < 10 ",
    ReFitPVs=True
)
#
mW = Selection(
    'DiMuAndW',
    Algorithm=mW_,
    RequiredSelections=[W_Data, dimu]
)

# =============================================================================
# Good Z
# ========================================================================
gZ = FilterDesktop(
    Preambulo=EW_preambulo,
    Code="""
    CHILDCUT ( ( -1e+10 * GeV < ptCone  ) & ( -1e+10 * GeV < ptCone2 ) & ( -1e+10 * GeV < etCone  ) , 1 ) & 
    CHILDCUT ( ( -1e+10 * GeV < ptCone  ) & ( -1e+10 * GeV < ptCone2 ) & ( -1e+10 * GeV < etCone  ) , 2 ) 
    """
)
Z_Data = Selection(
    'Z',
    Algorithm=gZ,
    RequiredSelections=[Z_Strip]
)

# =============================================================================
## charm + Z
# ========================================================================
cZ_ = CombineParticles(
    #
    DecayDescriptors=[
        "[H_10 -> D0        Z0 ]cc",
        "[H_10 -> D+        Z0 ]cc",
        "[H_10 -> D_s+      Z0 ]cc",
        "[H_10 -> Lambda_c+ Z0 ]cc",
    ],
    CombinationCut=" AALL ",
    MotherCut=" PALL ",
    ReFitPVs=True
)
#
cZ = Selection(
    'CharmAndZ',
    Algorithm=cZ_,
    RequiredSelections=[Z_Data, charm]
)

# =============================================================================
## dimu + Z
# ========================================================================
mZ_ = CombineParticles(
    #
    DecayDescriptor=" H_10 -> J/psi(1S) Z0",
    #
    CombinationCut=" AALL ",
    MotherCut=" PALL ",
    ReFitPVs=True
)
#
mZ = Selection(
    'DiMuAndZ',
    Algorithm=mZ_,
    RequiredSelections=[Z_Data, dimu]
)

s_cW = SelectionSequence("CharmW", TopSelection=cW)
s_cZ = SelectionSequence("CharmZ", TopSelection=cZ)
s_mW = SelectionSequence("MuMuW", TopSelection=mW)
s_mZ = SelectionSequence("MuMuZ", TopSelection=mZ)


from PhysSelPython.Wrappers import MultiSelectionSequence
charmEW = MultiSelectionSequence(
    "CharmEW",
    Sequences=[s_cW,
               s_cZ,
               s_mW,
               s_mZ
               ]
)

from DSTWriters.Configuration import (SelDSTWriter,
                                      stripMicroDSTStreamConf,
                                      stripMicroDSTElements)

# =============================================================================
#
from Configurables import LoKi__CounterAlg as CounterAlg
cnt = CounterAlg(
    'CharmEWCounters',
    Location="Counters/CharmEW",
    Preambulo=[
        "from LoKiPhys.decorators import *",
        "from LoKiCore.functions  import *",
        "pion_cuts  = in_range ( 300 * MeV , PT , 120 * GeV ) & ( CLONEDIST > 5000 ) & ( TRCHI2DOF < 5 ) & ( TRGHOSTPROB < 0.5 ) & ( PERR2/P2 < 0.05**2 ) ",
        "gamma_cuts = in_range ( 300 * MeV , PT ,  10 * GeV )  ",
        "pions      = SOURCE   ( '/Event/Phys/StdAllNoPIDsPions/Particles'  ,  pion_cuts ) ",
        "gammas     = SOURCE   ( '/Event/Phys/StdLooseAllPhotons/Particles' , gamma_cuts ) ",
    ],
    Variables={
        "px_c": " pions  >> sum ( PX ) ",
        "py_c": " pions  >> sum ( PY ) ",
        "px_g": " gammas >> sum ( PX ) ",
        "py_g": " gammas >> sum ( PY ) ",
        "n_c": " pions  >> SIZE       ",
        "g_c": " gammas >> SIZE       ",
    }
)
from Configurables import DataOnDemandSvc
dod = DataOnDemandSvc()
dod.AlgMap['/Event/Counters/CharmEW'] = cnt


def myMicroDSTStreamConf(pack=True):
    conf = stripMicroDSTStreamConf(pack)
    conf.extraItems += ['/Event/Counters#1']
    conf.extraItems += ['/Event/Counters/CharmEW#1']
    return conf

# Configuration of SelDSTWriter
SelDSTWriterConf = {'default':    myMicroDSTStreamConf(pack=False)}
SelDSTWriterElements = {'default': stripMicroDSTElements(pack=False)}

udstWriter = SelDSTWriter(
    "MyDSTWriter",
    StreamConf=SelDSTWriterConf,
    MicroDSTElements=SelDSTWriterElements,
    OutputFileSuffix='EW',
    SelectionSequences=[charmEW]
)

# Read only fired events to speed up
from PhysConf.Filters import LoKi_Filters
fltrs = LoKi_Filters(
    STRIP_Code="HLT_PASS_RE('Stripping.*(WMu|Z02MuMu).*')"
)

## the_year = "2011"
the_year = "2012"

logger.info('Use Year: %s ' % the_year)

from Configurables import DaVinci
dv = DaVinci(
    EventPreFilters=fltrs.filters('Filters'),
    InputType='DST',
    DataType=the_year,
    EvtMax=-1,
    Lumi=True,
    HistogramFile="DVHistos.root",
    #
    PrintFreq=50000
)
#
from Configurables import CondDB
CondDB(LatestGlobalTagByDataType=the_year)
#

dv.appendToMainSequence([udstWriter.sequence()])

# =============================================================================
if '__main__' == __name__:

    print 80 * '*'
    print __doc__
    print ' Author  : ', __author__
    print ' Version : ', __version__
    print ' Date    : ', __date__
    print 80 * '*'

# Input data
## from GaudiConf import IOHelper
# IOHelper().inputFiles([
# 'PFN:root://eoslhcb.cern.ch//eos/lhcb/LHCb/Collision12/EW.DST/00020241/0000/00020241_00000139_1.ew.dst'
# ] )

## dv.EvtMax = 1000

# =============================================================================
# The END
# =============================================================================
