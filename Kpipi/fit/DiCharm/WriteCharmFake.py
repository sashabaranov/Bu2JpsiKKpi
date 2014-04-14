#!/usr/bin/env gaudirun.py
# =============================================================================
# $Id:$
# =============================================================================
# @file WriteCharmFake.py
#
# Helper script to write muDST with Charm+W/Z0
#
# @author Vanya BELYAEV Ivan.Belyaev@itep.ru
# @date 2012-07-19
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
__date__ = "2012-05-02"
__version__ = "$Revision:$"
# =============================================================================

from GaudiKernel.SystemOfUnits import MeV, GeV

from PhysSelPython.Wrappers import AutomaticData, Selection, SelectionSequence

from StandardParticles import StdLooseAllPhotons  # needed for chi_b
from GaudiConfUtils.ConfigurableGenerators import FilterDesktop

#
# All "muons" (needed for fakes)
#
from StandardParticles import StdAllNoPIDsMuons
Fake_Strip = StdAllNoPIDsMuons

#
# get the prompt charm
#
from StrippingSelections.StrippingPromptCharm import StrippingPromptCharmConf as STRIP
#
strip = STRIP('PromptCharm', {
    'pT(D0)':  0.95 * GeV,  # pt-cut for  prompt   D0
    'pT(D0) for D*+':  0.95 * GeV,  # pt-cut for  D0 from  D*+
    'pT(D+)':  0.95 * GeV,  # pt-cut for  prompt   D+
    'pT(Ds+)':  0.95 * GeV,  # pt-cut for  prompt   Ds+
    'pT(Lc+)':  0.95 * GeV,  # pt-cut for  prompt   Lc+
    #
})
#
# the prompt charm selection
#
charm = strip.PromptCharm()
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
    "etCone     =    SINFO (  55002 , etCone_  , True ) ",
]

# ========================================================================
# good W
# ========================================================================
from GaudiConfUtils.ConfigurableGenerators import FilterDesktop
gW = FilterDesktop(
    Preambulo=EW_preambulo,
    Code="""
    in_range ( 15 * GeV , PT , 100 * GeV ) & INMUON & 
    ( -1e+10 * GeV < ptCone  ) &
    ( -1e+10 * GeV < ptCone2 ) &
    ( -1e+10 * GeV < etCone  ) 
    """
)
W_Data = Selection(
    'W',
    Algorithm=gW,
    RequiredSelections=[Fake_Strip]
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

s_cW = SelectionSequence("CharmFake", TopSelection=cW)

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
        "pion_cuts  = in_range ( 300 * MeV , PT , 120 * GeV ) & ( CLONEDIST > 5000 ) & ( TRCHI2DOF < 5 ) ",
        "gamma_cuts = in_range ( 300 * MeV , PT ,  10 * GeV )  ",
        "pions      = SOURCE ( '/Event/Phys/StdAllNoPIDsPions/Particles'  ,  pion_cuts ) ",
        "gammas     = SOURCE ( '/Event/Phys/StdLooseAllPhotons/Particles' , gamma_cuts ) ",
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
    OutputFileSuffix='CFake',
    ## SelectionSequences = [ charmEW ]
    SelectionSequences=[s_cW]
)

# Read only fired events to speed up
from PhysConf.Filters import LoKi_Filters
fltrs = LoKi_Filters(
    STRIP_Code=" HLT_PASS_RE('Stripping.*(D0|Dh|D2|Ds|Lc|Lambdac)*.*')  "
    # , VOID_Code  = " SKIP ( 50 ) "
)

the_year = "2011"
## the_year = "2012"

fltrs = fltrs.filters('Filters')
fltrs.reverse()

from Configurables import DaVinci
dv = DaVinci(
    EventPreFilters=fltrs,
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
# CHARM COMPLETE EVENT
# 'PFN:root://eoslhcb.cern.ch//eos/lhcb/LHCb/Collision11/CHARMCOMPLETEEVENT.DST/00022761/0000/00022761_00005178_1.charmcompleteevent.dst'
# ] )

## dv.EvtMax     = 10000
## dv.PrintFreq  = 500

# =============================================================================
# The END
# =============================================================================
