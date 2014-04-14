#!/usr/bin/env gaudirun.py
# =============================================================================
# $Id:$
# =============================================================================
# @file WriteDsPsi.py
#
# Helper script to write muDST with DsPsi
#
# @author Vanya BELYAEV Ivan.Belyaev@itep.ru
# @date 2012-07-19
#
#                   $Revision:$
# Last modification $Date$
#                by $Author$
# =============================================================================
"""
Helper script to write MicroDST with DsPsiselection 
"""
# =============================================================================
__author__ = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__ = "2012-05-02"
__version__ = "$Revision:$"
# =============================================================================

from GaudiKernel.SystemOfUnits import MeV, GeV

from PhysSelPython.Wrappers import AutomaticData, Selection, SelectionSequence

## needed for chi_b
from StandardParticles import StdLooseAllPhotons
from GaudiConfUtils.ConfigurableGenerators import FilterDesktop

#
# dimuon locations in DIMUON.DST
#
Jpsi_location = '/Event/Dimuon/Phys/FullDSTDiMuonJpsi2MuMuDetachedLine/Particles'
from PhysSelPython.Wrappers import AutomaticData
jpsi = AutomaticData(Location=Jpsi_location)
#
# get the prompt charm
#
from StrippingSelections.StrippingPromptCharm import StrippingPromptCharmConf as PC
#
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

# =============================================================================
# prepare D0
# =============================================================================
D0_alg = CombineParticles(
    # the decays to be reconstructed
    DecayDescriptor="[D0  -> K-  pi+]cc ",
    #
    Preambulo=Preambulo,
    # combination cut : wide mass-cut & PT-cut
    CombinationCut="""
    ( ADAMASS('D0') <  80 * MeV ) &
    ( APT           > 900 * MeV )
    """ ,
    # mother cut
    MotherCut="""
    ( chi2vx       <  9       ) &
    ( PT           >  1 * GeV ) &
    ( ADMASS('D0') < 75 * MeV ) 
    """ ,
    #
    ParticleCombiners={'': 'LoKi::VertexFitter'}
)

# make selection
from PhysSelPython.Wrappers import Selection
D0_sel = Selection(
    'D0ForPsiC',
    Algorithm=D0_alg,
    RequiredSelections=[kaons, pions]
)

# =============================================================================
# prepare D+
# =============================================================================
Dp_alg = CombineParticles(
    # the decays to be reconstructed
    DecayDescriptor="[D+  -> K-  pi+ pi+]cc ",
    #
    Preambulo=Preambulo,
    # combination cut : wide mass-cut & PT-cut
    CombinationCut="""
    ( ADAMASS('D+') <  60 * MeV ) &
    ( APT           > 900 * MeV )
    """ ,
    # mother cut
    MotherCut="""
    ( chi2vx       < 25       ) &
    ( PT           >  1 * GeV ) &
    ( ADMASS('D+') < 50 * MeV ) 
    """ ,
    #
    ParticleCombiners={'': 'LoKi::VertexFitter'}
)

# make selection
Dp_sel = Selection(
    'DpForPsiC',
    Algorithm=Dp_alg,
    RequiredSelections=[kaons, pions]
)

# =============================================================================
# prepare Ds+
# =============================================================================
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


# make selection
Ds_sel = Selection(
    'DsForPsiC',
    Algorithm=Ds_alg,
    RequiredSelections=[kaons, pions]
)

# =============================================================================
BD0_alg = CombineParticles(
    # the decays to be reconstructed
    DecayDescriptor="[ B0 -> J/psi(1S) D0 ]cc ",
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
# ========================================================================
BDp_alg = CombineParticles(
    # the decays to be reconstructed
    DecayDescriptor="[ B+ -> J/psi(1S) D+ ]cc ",
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
BD0_sel = Selection(
    'PsiD0',
    Algorithm=BD0_alg,
    RequiredSelections=[jpsi, D0_sel]
)

# make the selection
BDp_sel = Selection(
    'PsiDp',
    Algorithm=BDp_alg,
    RequiredSelections=[jpsi, Dp_sel]
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

Bc_1PI = SelectionSequence("PSIPi", TopSelection=Bc_1p)
Bc_3PI = SelectionSequence("PSI3Pi", TopSelection=Bc_3p)
B_D0 = SelectionSequence("PSID0", TopSelection=BD0_sel)
B_Dp = SelectionSequence("PSIDp", TopSelection=BDp_sel)
B_Ds = SelectionSequence("PSIDs", TopSelection=BDs_sel)

#
# selection sequence
#
from PhysSelPython.Wrappers import MultiSelectionSequence
B_ALL = MultiSelectionSequence(
    "B2PSIC",
    Sequences=[B_D0,
               B_Dp,
               B_Ds,
               Bc_3PI,
               Bc_1PI]
)

from DSTWriters.Configuration import (SelDSTWriter,
                                      stripMicroDSTStreamConf,
                                      stripMicroDSTElements)

# Configuration of SelDSTWriter
SelDSTWriterConf = {'default': stripMicroDSTStreamConf(pack=False)}
SelDSTWriterElements = {'default': stripMicroDSTElements(pack=False)}

udstWriter = SelDSTWriter(
    "MyMicroDSTWriter",
    StreamConf=SelDSTWriterConf,
    MicroDSTElements=SelDSTWriterElements,
    OutputFileSuffix='PsiC',
    SelectionSequences=[B_ALL]
)

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

## the_year = "2011"
the_year = "2012"

from Configurables import DaVinci
dv = DaVinci(
    EventPreFilters=fltrs_0.filters(
        'Filters0') + fltrs.filters('Filters'),
    InputType='DST',
    DataType=the_year,
    EvtMax=-1,
    Lumi=True,
    HistogramFile="DVHistos.root",
    # dbase
    ## DDDBtag       = "head-20120413"  ,
    ## CondDBtag     = "head-20120724"  ,
    #
    PrintFreq=1000
)
#
from Configurables import CondDB
CondDB(LatestGlobalTagByDataType=the_year)


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
## PFN = 'PFN:castor:/castor/cern.ch/grid/'
# IOHelper().inputFiles([
# stripping 19
##     PFN + '/lhcb/LHCb/Collision12/DIMUON.DST/00018987/0000/00018987_00000094_1.dimuon.dst' ,
##     PFN + '/lhcb/LHCb/Collision12/DIMUON.DST/00018987/0000/00018987_00000099_1.dimuon.dst' ,
##     PFN + '/lhcb/LHCb/Collision12/DIMUON.DST/00018987/0000/00018987_00000108_1.dimuon.dst' ,
##     PFN + '/lhcb/LHCb/Collision12/DIMUON.DST/00018987/0000/00018987_00000109_1.dimuon.dst' ,
##     PFN + '/lhcb/LHCb/Collision12/DIMUON.DST/00018987/0000/00018987_00000110_1.dimuon.dst' ,
##     PFN + '/lhcb/LHCb/Collision12/DIMUON.DST/00018987/0000/00018987_00000777_1.dimuon.dst'
# ] )

## dv.EvtMax    = 100000
## dv.PrintFreq =    100

# =============================================================================
# The END
# =============================================================================
