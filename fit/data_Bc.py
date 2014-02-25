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
tSelection5 = ROOT.TChain('JpsiKKpi/t')
tSelection6 = ROOT.TChain('JpsiKKpi/t')

#
tLumi = ROOT.TChain('GetIntegratedLuminosity/LumiTuple')

outputdir = '/afs/cern.ch/user/a/albarano/cmtuser/Bender_v22r8/Scripts/testing/output/'

files = ['B2JpsiKKpi-2011.root',  'B2JpsiKKpi-2012.root']
files_sel6 = ['B2JpsiKKpi-2011-sel6.root',  'B2JpsiKKpi-2012-sel6.root']


nFiles = len(files)

for f in files:
    tSelection5.Add(outputdir + f)

for f in files_sel6:
    tSelection6.Add(outputdir + f)    


# mc_outputdir = '/afs/cern.ch/user/a/albarano/cmtuser/Bender_v22r8/Scripts/testing/MC/output/'
# tBu_mc = ROOT.TChain('Bplus/t')
# mc_files = ['2011-Kpipi-Pythia6.root',
#          '2011-Kpipi-Pythia8.root',
#          '2012-Kpipi-Pythia6.root', 
#          '2012-Kpipi-Pythia8.root']

# for f in mc_files:
#     tBu_mc.Add(mc_outputdir + f)






# =============================================================================
# get the luminosity
lumi = getLumi(tLumi)
logger.info(" Luminosity: %s #files %d" % (lumi, nFiles))
logger.info(" Entries  B+ -> J/psi KK pi+: %s" % len(tSelection5))


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
