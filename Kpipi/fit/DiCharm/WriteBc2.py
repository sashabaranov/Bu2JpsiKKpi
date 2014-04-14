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

#
# get the prompt charm
#
from StrippingSelections.StrippingPromptCharm import StrippingPromptCharmConf as PC
#
# ======================================
pc = PC('PromptCharm', {
    'TrackCuts'       : """
    ( TRCHI2DOF   < 5   ) &
    ( TRGHOSTPROB < 1.0 ) &           
    ( PT > 250 * MeV    ) &
    in_range  ( 2 , ETA , 5 ) 
    """ ,
    'BasicCuts': ' ',
    'KaonCuts': ' & in_range ( 3.2 * GeV , P , 150 * GeV ) & ( -5 < PIDK  - PIDpi ) ',
    'PionCuts': ' & in_range ( 3.2 * GeV , P , 150 * GeV ) & ( -5 < PIDpi - PIDK  ) ',
}
)

pions = pc.pions()
kaons = pc.kaons()

Preambulo = [
    # shortcut for chi2 of vertex fit
    'chi2vx = VFASPF(VCHI2) ',
    # shortcut for the c*tau
    "from GaudiKernel.PhysicalConstants import c_light",
    # no embedded cut!
    "ctau   = BPVLTIME () * c_light "  # ATTENTION
]

# ========================================================================
# B+ -> J/psi + K+
# ========================================================================
from GaudiConfUtils.ConfigurableGenerators import CombineParticles
bp_K = CombineParticles(
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
    ( chi2vx  < 25 ) &
    ( ctau    >  0 )
    """ ,
    ParticleCombiners={'': 'LoKi::VertexFitter'},
    ReFitPVs=True
)
#
from PhysSelPython.Wrappers import Selection
Bp_K = Selection(
    'Bp2PsiK',
    Algorithm=bp_K,
    RequiredSelections=[jpsi, kaons]
)
# ===========================================================================
# Bc -> B K pi
# ===========================================================================
bc = CombineParticles(
    #
    DecayDescriptors=[
        '[ B_c+ -> B+ K+ pi-]cc',
        '[ B_c+ -> B+ K- pi+]cc',
    ],
    #
    Preambulo=Preambulo,
    #
    CombinationCut="""
    in_range ( 5.9 * GeV , AM ,  6.6 * GeV ) 
    """ ,
    #
    MotherCut="""
    in_range  ( 5.95 * GeV , M , 6.55 * GeV ) &
    ( chi2vx  < 36 ) &
    ( ctau    >  0 )
    """ ,
    ParticleCombiners={'': 'LoKi::VertexFitter'},
    ReFitPVs=True
)
#
from PhysSelPython.Wrappers import Selection
Bc = Selection(
    'Bc',
    Algorithm=bc,
    RequiredSelections=[Bp_K, kaons, pions]
)

from PhysSelPython.Wrappers import SelectionSequence
Bc_seq = SelectionSequence("BC", TopSelection=Bc)

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
    OutputFileSuffix='Bc2',
    SelectionSequences=[Bc_seq]
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
#
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
    DDDBtag="head-20120413",
    CondDBtag="head-20120724",
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

## dv.EvtMax    = 50000
## dv.PrintFreq =   100

# =============================================================================
# The END
# =============================================================================
