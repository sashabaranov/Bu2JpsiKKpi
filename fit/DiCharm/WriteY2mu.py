#!/usr/bin/env gaudirun.py
# =============================================================================
# $Id:$
# =============================================================================
# @file WriteCharmEW.py
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
# 2x2mu-locations  in LEPTONIC.MDST
#  Stripping 17:
Mu22_Location = '/Event/Leptonic/Phys/DoubleDiMuonForCharmAssociative/Particles'
from PhysSelPython.Wrappers import AutomaticData
Mu22_Strip = AutomaticData(Location=Mu22_Location)

Y2mu_alg = FilterDesktop(
    Code="""
    DECTREE ( ' Meson -> ( Meson ->  mu+ mu- ) ( Meson -> mu+ mu- ) ' ) 
    & (  ( 7 * GeV < M1 )  |  ( 7 * GeV < M2 ) ) 
    """
)
Y2mu = Selection(
    'Ymu2',
    Algorithm=Y2mu_alg,
    RequiredSelections=[Mu22_Strip]
)

seq = SelectionSequence("Y2Mu", TopSelection=Y2mu)

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
    OutputFileSuffix='Y2Muon',
    SelectionSequences=[seq]
)

# Read only fired events to speed up
from PhysConf.Filters import LoKi_Filters
fltrs = LoKi_Filters(
    VOID_Code=" 0.5 < CONTAINS('%s') " % Mu22_Location,
    STRIP_Code=" HLT_PASS_RE('Stripping.*DoubleDiMuon.*')"
)


the_year = "2011"
## the_year = "2012"

from Configurables import EventNodeKiller
killer = EventNodeKiller("KillpRec")
killer.Nodes += ['/Event/pRec']

from Configurables import DaVinci
dv = DaVinci(
    EventPreFilters=fltrs.filters('Filters'),
    InputType='MDST',
    DataType=the_year,
    EvtMax=-1,
    Lumi=True,
    HistogramFile="DVHistos.root",
    # dbase
    # DDDBtag       = "head-20120413"  ,
    # CondDBtag     = "head-20120724"  ,
    #
    PrintFreq=1000
)
#
from Configurables import CondDB
CondDB(LatestGlobalTagByDataType=the_year)

dv.appendToMainSequence([killer] + [udstWriter.sequence()])


# =============================================================================
if '__main__' == __name__:

    print 80 * '*'
    print __doc__
    print ' Author  : ', __author__
    print ' Version : ', __version__
    print ' Date    : ', __date__
    print 80 * '*'

# Input data
# from GaudiConf import IOHelper
# IOHelper().inputFiles(
#     [ 'PFN:castor:/castor/cern.ch/grid/lhcb/LHCb/Collision11/LEPTONIC.MDST/00012565/0000/00012565_00000063_1.leptonic.mdst' ]
#    )

# =============================================================================
# The END
# =============================================================================
