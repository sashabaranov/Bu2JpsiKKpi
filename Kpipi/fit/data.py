import os

import ROOT
from AnalysisPython.PyRoUts import *
from ostap.data import Data, DataAndLumi
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
OUTPUT_DIR = '$Kpipidir/output/'

def append_dir(directory, files):
    return [os.path.join(directory, f) for f in files]


def from_ganga(job_id):
    from glob import glob
    prefix = "/afs/cern.ch/user/a/albarano/gangadir/workspace/albarano/LocalXML/"
    pattern = prefix + str(job_id) + "/*/output/*.root"
    return glob(pattern)



# =============================================================================
files7 = from_ganga(187) + from_ganga(188)

selection7 = DataAndLumi(branch='JpsiKpipi/t', files=files7)

logger.info("Selection 7: %s" % selection7)



# =============================================================================
mc_files = {
    '2011': {
        'Pythia6': from_ganga(182) + from_ganga(183),
        'Pythia8': from_ganga(180) + from_ganga(181),
    },

    '2012': {
        'Pythia6': from_ganga(184) + from_ganga(185),
        'Pythia8': from_ganga(178) + from_ganga(179),
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
# get the luminosity
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
