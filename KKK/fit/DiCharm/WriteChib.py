#!/usr/bin/env gaudirun.py
# =============================================================================
# $Id:$
# =============================================================================
# @file WriteChib.py
#
# Helper script to write MicroDST (with chib selection)
# for Sascha Mazurov
#
# @author Vanya BELYAEV Ivan.Belyaev@itep.ru
# @date 2012-05-02
#
#                   $Revision:$
# Last modification $Date$
#                by $Author$
# =============================================================================
"""

Helper script to write MicroDST (with chib selection)
for Sascha Mazurov

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
# dimuon in stripping DST
#
DiMuLocation = '/Event/Dimuon/Phys/FullDSTDiMuonDiMuonHighMassLine/Particles'
from PhysSelPython.Wrappers import AutomaticData
DiMuData = AutomaticData(Location=DiMuLocation)

# =============================================================================
# Upsilon -> mumu, cuts by  Giulia Manca
# =============================================================================
UpsAlg = FilterDesktop(
    Code="""
    ( M > 7 * GeV ) & 
    DECTREE   ('Meson -> mu+ mu-'  )                      &
    CHILDCUT( 1 , HASMUON & ISMUON )                      &
    CHILDCUT( 2 , HASMUON & ISMUON )                      & 
    ( MINTREE ( 'mu+' == ABSID , PT ) > 1 * GeV         ) &
    ( MAXTREE ( ISBASIC & HASTRACK , TRCHI2DOF ) < 4    ) & 
    ( MINTREE ( ISBASIC & HASTRACK , CLONEDIST ) > 5000 ) & 
    ( VFASPF  ( VPCHI2 ) > 0.5/100 ) 
    & ( abs ( BPV ( VZ    ) ) <  0.5 * meter     ) 
    & (       BPV ( vrho2 )   < ( 10 * mm ) ** 2 ) 
    """ ,
    Preambulo=[
        "vrho2 = VX**2 + VY**2"
    ],
    ReFitPVs=True
)
UpsSel = Selection(
    'Upsilon',
    Algorithm=UpsAlg,
    RequiredSelections=[DiMuData]
)

Ups = SelectionSequence("Upsilons", TopSelection=UpsSel)

# =============================================================================
# chi_b -> Upsilon gamma
# =============================================================================
from GaudiConfUtils.ConfigurableGenerators import CombineParticles
ChibAlg = CombineParticles(
    DecayDescriptor="chi_b1(1P) -> J/psi(1S) gamma",
    DaughtersCuts={
        "gamma": " ( 350 * MeV < PT ) & ( CL > 0.01 )  "
    },
    CombinationCut="""
    ( AM - AM1 ) < 3 * GeV 
    """ ,
    MotherCut=" PALL",
    #
    # we are dealing with photons!
    #
    ParticleCombiners={
        '': 'LoKi::VertexFitter'
    }
)
ChibSel1 = Selection(
    'PreSelChib',
    Algorithm=ChibAlg,
    RequiredSelections=[UpsSel, StdLooseAllPhotons]
)
from GaudiConfUtils.ConfigurableGenerators import Pi0Veto__Tagger
TagAlg = Pi0Veto__Tagger(
    ExtraInfoIndex=25001,  # should be unique!
    MassWindow=20 * MeV,  # cut on delta-mass
    MassChi2=-1,  # no cut for chi2(mass)
)
ChibSel2 = Selection(
    'Chi_b',
    Algorithm=TagAlg,
    RequiredSelections=[ChibSel1]
)

Chib = SelectionSequence("ChiB", TopSelection=ChibSel2)

from PhysSelPython.Wrappers import MultiSelectionSequence
Upsilons = MultiSelectionSequence(
    "Bottomonia",
    Sequences=[Ups, Chib]
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
    OutputFileSuffix='DiMuon',
    SelectionSequences=[Upsilons]
)

#
# Read only fired events to speed up
from PhysConf.Filters import LoKi_Filters
fltrs = LoKi_Filters(
    STRIP_Code="HLT_PASS_RE('Stripping.*DiMuon.*HighMass.*')"
)

the_year = "2011"

from Configurables import DaVinci
dv = DaVinci(
    EventPreFilters=fltrs.filters('Filters'),
    InputType='DST',
    DataType=the_year,
    EvtMax=-1,
    Lumi=True,
    HistogramFile="DVHistos.root",
    # dbase
    DDDBtag="head-20120413",
    CondDBtag="head-20120420",
    #
    PrintFreq=1000
)

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
# 'PFN:castor:/castor/cern.ch/grid/lhcb/LHCb/Collision11/DIMUON.DST/00016768/0000/00016768_00000036_1.dimuon.dst'
# ] )

# =============================================================================
# The END
# =============================================================================
