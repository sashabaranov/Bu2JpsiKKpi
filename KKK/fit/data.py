import ROOT
import os
import glob
from AnalysisPython.PyRoUts import *
import AnalysisPython.ZipShelve as ZipShelve
from ostap.data import Data, DataAndLumi
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
OUTPUT_DIR = '$KKKdir/output/'

def append_dir(directory, files):
    return [os.path.join(directory, f) for f in files]

def from_ganga(job_id):
    import os
    from glob import glob
    prefix = "$HOME/gangadir/workspace/albarano/LocalXML/"
    pattern = prefix + str(job_id) + "/*/output/*.root"
    return glob(os.path.expandvars(pattern))



# =============================================================================
# files5 = append_dir(OUTPUT_DIR, ['RD-2011.root',  'RD-2012.root'])
# files6 = append_dir(OUTPUT_DIR, ['RD-2011-sel6.root',  'RD-2012-sel6.root'])
files7 = from_ganga(157) + from_ganga(158)


# selection5 = DataAndLumi(branch='JpsiKKK/t', files=files5)
# selection6 = DataAndLumi(branch='JpsiKKK/t', files=files6)
selection7 = DataAndLumi(branch='JpsiKKK/t', files=files7)


# logger.info("Selection 5: %s" % selection5)
# logger.info("Selection 6: %s" % selection6)
logger.info("Selection 7: %s" % selection7)


# tSelection5 = selection5.data
# tSelection6 = selection6.data
tSelection7 = selection7.data


# =============================================================================
mc_files = {
    '2011': {
        'Pythia6': from_ganga(161) + from_ganga(162),
        'Pythia8': from_ganga(163) + from_ganga(164),
    },

    '2012': {
        'Pythia6': from_ganga(159) + from_ganga(160),
        'Pythia8': from_ganga(165) + from_ganga(166),
    }
}


mc_Pythia6 = Data(branch='Bplus/t', files=mc_files['2011']['Pythia6'] + mc_files['2012']['Pythia6'])
mc_Pythia8 = Data(branch='Bplus/t', files=mc_files['2011']['Pythia8'] + mc_files['2012']['Pythia8'])

mc_all_files = []
for year, pythias in mc_files.items():
    for pythia, files in pythias.items():
        mc_all_files += files

mc_total = Data(branch='Bplus/t', files=mc_all_files)


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
