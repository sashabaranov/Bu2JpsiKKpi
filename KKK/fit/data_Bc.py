import ROOT
import os
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




# =============================================================================
class DataAndLumi(object):
    """
    Contains data and luminosity info inside.
    """

    def __init__(self, 
                branch, 
                lumi_branch='GetIntegratedLuminosity/LumiTuple',
                files=[]):
        self.data = ROOT.TChain(branch)
        self.lumi = ROOT.TChain(lumi_branch)
        self.filelist = []

        if files:
            self.add_files(files)
    
    def add_files(self, files):
        self.filelist += files

        for f in files:
            self.data.Add(f)
            self.lumi.Add(f)

    def get_luminosity(self):
        return getLumi(self.lumi)

    def __str__(self):
        ret = "<"
        ret += "Luminosity: {}; #files: {}; ".format(self.get_luminosity(), len(self.filelist))
        ret += "Entries: {}>".format(len(self.data))

        return ret




# =============================================================================
OUTPUT_DIR = '$KKKdir/output/'

def append_dir(directory, files):
    return [os.path.join(directory, f) for f in files]


data5 = append_dir(OUTPUT_DIR, ['RD-2011.root',  'RD-2012.root'])
data6 = append_dir(OUTPUT_DIR, ['RD-2011-sel6.root',  'RD-2012-sel6.root'])

selection5 = DataAndLumi(branch='JpsiKKK/t', files=data5)
selection6 = DataAndLumi(branch='JpsiKKK/t', files=data6)



logger.info("Selection 5: %s" % selection5)
logger.info("Selection 6: %s" % selection6)

tSelection5 = selection5.data
tSelection6 = selection6.data


# =============================================================================
mc_outputdir = '$KKKdir/MC/output/'

tBu_mc = ROOT.TChain('Bplus/t')
tBu_mc_p6 = ROOT.TChain('Bplus/t')
tBu_mc_p8 = ROOT.TChain('Bplus/t')

mc_pythia6 = ['2011-KKK-Pythia6.root', '2012-KKK-Pythia6.root']
mc_pythia8 = ['2011-KKK-Pythia8.root', '2012-KKK-Pythia8.root']

mc_files = mc_pythia6 + mc_pythia8

for f in mc_pythia6:
    tBu_mc_p6.Add(mc_outputdir + f)

for f in mc_pythia8:
    tBu_mc_p8.Add(mc_outputdir + f)

for f in mc_files:
    tBu_mc.Add(mc_outputdir + f)




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
