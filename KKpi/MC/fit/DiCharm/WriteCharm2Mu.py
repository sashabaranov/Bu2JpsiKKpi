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
# Charm&2mu-locations  in EW.DST
#
C2Mu_Location = '/Event/Charm/Phys/DiMuonAndCharmForPromptCharm/Particles'
from PhysSelPython.Wrappers import AutomaticData
C2Mu_Strip = AutomaticData(Location=C2Mu_Location)


seq = SelectionSequence("Charm2Mu", TopSelection=C2Mu_Strip)

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
    OutputFileSuffix='Charm2Mu',
    SelectionSequences=[seq]
)

# Read only fired events to speed up
from PhysConf.Filters import LoKi_Filters
fltrs = LoKi_Filters(
    VOID_Code="""
    0.5 < CONTAINS('/Event/Charm/Phys/DiMuonAndCharmForPromptCharm/Particles')
    """ ,
    STRIP_Code="HLT_PASS_RE('Stripping.*(DiMuonAndCharm).*')"
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
## from GaudiConf import IOHelper
# IOHelper().inputFiles([
# 'PFN:castor:/castor/cern.ch/grid/lhcb/LHCb/Collision11/CHARM.MDST/00012561/0000/00012561_00000001_1.charm.mdst'
# /lhcb/LHCb/Collision11/EW.DST/00012564/0000/00012564_00000001_1.ew.dst'
# ] )

# =============================================================================
# The END
# =============================================================================
