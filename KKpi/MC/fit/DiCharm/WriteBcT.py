#!/usr/bin/env gaudirun.py
# =============================================================================
# $Id:$
# =============================================================================
# @file WriteBcT.py
#
# Helper script to write muDST
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
from GaudiKernel.SystemOfUnits import MeV, GeV, micrometer
from PhysSelPython.Wrappers import AutomaticData, Selection, SelectionSequence
from GaudiConfUtils.ConfigurableGenerators import FilterDesktop

#
# dimuon locations in DIMUON.DST
#
Jpsi_det_location = '/Event/Dimuon/Phys/FullDSTDiMuonJpsi2MuMuDetachedLine/Particles'
Jpsi_unb_location = '/Event/Dimuon/Phys/FullDSTDiMuonJpsi2MuMuTOSLine/Particles'

from PhysSelPython.Wrappers import AutomaticData
jpsi_det = AutomaticData(Location=Jpsi_det_location)
jpsi_unb = AutomaticData(Location=Jpsi_unb_location)

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
    ( PROBNNpi > 0.10      ) 
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
    ( PROBNNk > 0.10       ) 
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
    ( chi2vx < 25  ) & mbc_cut &  ( ctau_25 > %s )
    """ %  ( 40 * micrometer )
)

# B+ :
alg_bu = CombineParticles(
    DecayDescriptor="[B+ -> J/psi(1S) K+ ]cc",
    Preambulo=Preambulo,
    CombinationCut=" mbp_acut ",
    MotherCut="""
    ( chi2vx < 25  ) & mbp_cut &  ( ctau_25 > %s )
    """ %  ( 40 * micrometer )
)


sel_bc_det = Selection(
    "TheBc_det",
    Algorithm=alg_bc,
    RequiredSelections=[jpsi_det, pions]
)

sel_bc_unb = Selection(
    "TheBc_unb",
    Algorithm=alg_bc,
    RequiredSelections=[jpsi_unb, pions]
)


sel_bu_det = Selection(
    "TheB_det",
    Algorithm=alg_bu,
    RequiredSelections=[jpsi_det, kaons]
)


sel_bu_unb = Selection(
    "TheB_unb",
    Algorithm=alg_bu,
    RequiredSelections=[jpsi_unb, kaons]
)

Bc_det_seq = SelectionSequence("Bc_det", TopSelection=sel_bc_det)
Bu_det_seq = SelectionSequence("Bu_det", TopSelection=sel_bu_det)
Bc_unb_seq = SelectionSequence("Bc_unb", TopSelection=sel_bc_unb)
Bu_unb_seq = SelectionSequence("Bu_unb", TopSelection=sel_bu_unb)

#
# selection sequence
#
from PhysSelPython.Wrappers import MultiSelectionSequence
B_SEQ = MultiSelectionSequence(
    "B2PSI",
    Sequences=[Bc_det_seq,
               Bu_det_seq,
               Bc_unb_seq,
               Bu_unb_seq]
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
    OutputFileSuffix='BSEQ',
    SelectionSequences=[B_SEQ]
)

# Read only fired events to speed up
from PhysConf.Filters import LoKi_Filters
fltrs = LoKi_Filters(
    STRIP_Code="HLT_PASS_RE('Stripping.*DiMuonJpsi2MuMu.*')",
    VOID_Code="""
    ( 0.5 < CONTAINS('%s') ) | 
    ( 0.5 < CONTAINS('%s') ) 
    """ %  ( Jpsi_det_location , Jpsi_unb_location )
)
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
    #
    PrintFreq=100000
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
# IOHelper().inputFiles([
# 'PFN:root://eoslhcb.cern.ch//eos/lhcb/LHCb/Collision11/DIMUON.DST/00022727/0000/00022727_0000%04d_1.dimuon.dst' % i for i in range(7,1000)
# ] )
#
## dv.EvtMax    = 10000
## dv.PrintFreq = 100

# =============================================================================
# The END
# =============================================================================
