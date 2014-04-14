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

the_year = "2012"

from Configurables import DaVinci
dv = DaVinci(
    InputType='DST',
    DataType=the_year,
    EvtMax=-1,
    Lumi=True,
    HistogramFile="DVHisto.root",
    TupleFile="DVTuple.root",
    # dbase
    DDDBtag="head-20120413",
    CondDBtag="head-20120724",
    #
    PrintFreq=1000
)
#
from Configurables import CondDB
CondDB(LatestGlobalTagByDataType=the_year)

# =============================================================================
if '__main__' == __name__:

    print 80 * '*'
    print __doc__
    print ' Author  : ', __author__
    print ' Version : ', __version__
    print ' Date    : ', __date__
    print 80 * '*'

## import shelve
## db     = shelve.open('/afs/cern.ch/user/i/ibelyaev/public/LFNs/lfns.db','r')
## giulia = db ['Giulia/ICHEP']
# db.close()

# Input data
## from GaudiConf import IOHelper
## PFN = 'PFN:castor:/castor/cern.ch/grid/'

## IOHelper().inputFiles ( [ PFN + i for i in giulia ]  )

## dv.EvtMax = 5000
