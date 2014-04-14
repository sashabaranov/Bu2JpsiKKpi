#!/usr/bin/env gaudirun.py
# =============================================================================
# $Id:$
# =============================================================================
# @file WriteX.py
#
# Helper script to write muDST with  B->X(3782)K candidates
#
# @author Vanya BELYAEV Ivan.Belyaev@itep.ru
# @date 2013-06-08
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


#
# J/psi
#
jpsi_name = 'FullDSTDiMuonJpsi2MuMuDetachedLine'
from PhysSelPython.Wrappers import AutomaticData
jpsi = AutomaticData('/Event/Dimuon/Phys/%s/Particles' % jpsi_name)


#from StandardParticles import StdAllNoPIDsPions as _pions
#from StandardParticles import StdAllNoPIDsKaons as _kaons
#from StandardParticles import StdAllLoosePions as _pions
#from StandardParticles import StdAllLooseKaons as _kaons
from StandardParticles import StdLoosePions as _pions
from StandardParticles import StdLooseKaons as _kaons

#
# Pions
#
from GaudiConfUtils.ConfigurableGenerators import FilterDesktop
pion_alg = FilterDesktop(
    Code="""
    ( TRCHI2DOF   < 4   ) & 
    ( TRGHOSTPROB < 0.5 ) &
    ( 4 < MIPCHI2DV ()  ) 
    """
)
from PhysSelPython.Wrappers import Selection
pions = Selection(
    "Pions",
    Algorithm=pion_alg,
    RequiredSelections=[_pions])

#
# Kaons
#
from GaudiConfUtils.ConfigurableGenerators import FilterDesktop
kaon_alg = FilterDesktop(
    Code="""
    ( TRCHI2DOF   < 4   ) & 
    ( TRGHOSTPROB < 0.5 ) &
    ( 4 < MIPCHI2DV ()  ) 
    """
)
kaons = Selection(
    "Kaons",
    Algorithm=kaon_alg,
    RequiredSelections=[_kaons])


#
# make B+ -> X K + candidates
#
from GaudiConfUtils.ConfigurableGenerators import CombineParticles
b_alg = CombineParticles(
    DecayDescriptor="[ B+ -> J/psi(1S) pi+ pi- K+]cc",
    #
    CombinationCut="""
    ( in_range   ( 5.100 * GeV , AM , 5.450 * GeV ) 
    | in_range   ( 6.100 * GeV , AM , 6.450 * GeV ) ) &
    ADOCACHI2CUT ( 16 , '' )
    """ ,
    #
    MotherCut="""
    ( in_range ( 5.150 * GeV ,  M , 5.400 * GeV )
    | in_range ( 6.150 * GeV ,  M , 6.400 * GeV ) ) &   
    ( VFASPF   ( VCHI2PDOF ) < 10 ) &
    ( BPVLTIME ( 400 ) > 0.050 * millimeter / c_light ) 
    """
)
#
# convert algorithm into selection
#
from PhysSelPython.Wrappers import Selection
b_sel = Selection(
    "Bplus",
    Algorithm=b_alg,
    RequiredSelections=[jpsi, kaons, pions]
)


#
# finally build selection sequence
#
from PhysSelPython.Wrappers import SelectionSequence
b_seq = SelectionSequence("B", TopSelection=b_sel)

from DSTWriters.Configuration import (SelDSTWriter,
                                      stripMicroDSTStreamConf,
                                      stripMicroDSTElements)

# Configuration of SelDSTWriter
SelDSTWriterConf = {'default': stripMicroDSTStreamConf(pack=False)}
SelDSTWriterElements = {'default': stripMicroDSTElements(pack=False)}

udstWriter = SelDSTWriter(
    "MyDSTWriter",
    StreamConf=SelDSTWriterConf,
    MicroDSTElements=SelDSTWriterElements,
    OutputFileSuffix='BX',
    SelectionSequences=[b_seq]
)

#
# Read only fired events to speed up
#
from PhysConf.Filters import LoKi_Filters
fltrs1 = LoKi_Filters(
    STRIP_Code=" HLT_PASS_RE('Stripping.*%s.*') " % jpsi_name
)
fltrs2 = LoKi_Filters(
    VOID_Code="""
    ( RECSUMMARY (  0 , -1 ) > 0.5 ) & 
    ( RECSUMMARY ( 10 , -1 ) < 500 ) & 
    ( RECSUMMARY ( 13 , -1 ) < 500 )
    """
)

##the_year = "2011"
the_year = "2012"

from Configurables import DaVinci
dv = DaVinci(
    EventPreFilters=fltrs1.filters(
        'Filters1') + fltrs2.filters('Filters2'),
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
#from GaudiConf import IOHelper
# IOHelper().inputFiles([
#    'PFN:root://eoslhcb.cern.ch//eos/lhcb/LHCb/Collision11/DIMUON.DST/00022727/0000/00022727_0000%04d_1.dimuon.dst' % i for i in range(7,1000)
#    ] )
#
#dv.EvtMax    = 10000
#dv.PrintFreq = 100
