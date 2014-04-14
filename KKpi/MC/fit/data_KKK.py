import ROOT
import glob
from AnalysisPython.PyRoUts import *
from AnalysisPython.GetLumi import getLumi
import AnalysisPython.ZipShelve as ZipShelve
from AnalysisPython.Logger import getLogger


# =============================================================================
logger = getLogger(__name__)


#
# define the default storage for data set
#
RAD = ROOT.RooAbsData
if RAD.Tree != RAD.getDefaultStorageType():
    print 'DEFINE default storage type to be TTree! '
    RAD.setDefaultStorageType(RAD.Tree)


def _draw_(self, *args, **kwargs):
    t = self.store().tree()
    return t.Draw(*args, **kwargs)
ROOT.RooDataSet.Draw = _draw_


#
tBu2 = ROOT.TChain('MCBu2JpsiKKK/NTuple_for_B-mesons')


outputdir = '/afs/cern.ch/user/a/albarano/cmtuser/Bender_v22r8/Scripts/Bu2JpsiKKK/MC/output/'
files = ['2012-KKK-Pythia8.root']
nFiles = len(files)

for f in files:
    tBu2.Add(outputdir + f)


# =============================================================================


from GaudiKernel.SystemOfUnits import second, nanosecond
from GaudiKernel.PhysicalConstants import c_light

ct_Bu_PDG = VE(0.4911, (0.4911 * 0.011 / 1.638) ** 2)
ct_Bu_MC = 0.492


#
# finally make the class
#
# clname = 'BcTChain'
# logger.info ( 'Finally: prepare&load C++ class %s ' % clname )
# tBc1.MakeClass( clname )
# ROOT.gROOT.LoadMacro( clname + '.C+' )
