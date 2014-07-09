import os

import ROOT
from Ostap.PyRoUts import *
from Ostap.Data import Data, DataAndLumi
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
def append_dir(directory, files):
    return [os.path.join(directory, f) for f in files]


def from_ganga(job_id):
    import os
    from glob import glob
    prefix = "$WORKDIR/gangadir/workspace/albarano/LocalXML/"
    pattern = prefix + str(job_id) + "/*/output/*.root"
    return glob(os.path.expandvars(pattern))



# =============================================================================
files7 = from_ganga(483) + from_ganga(484)

selection7 = DataAndLumi('JpsiKpipi/t', files7)

logger.info("Selection 7: %s" % selection7)



# =============================================================================
mc_files = {
    '2011': {
        'Pythia6': from_ganga(508) + from_ganga(509),
        'Pythia8': from_ganga(506) + from_ganga(507),
    },

    '2012': {
        'Pythia6': from_ganga(510) + from_ganga(511),
        'Pythia8': from_ganga(504) + from_ganga(505),
    }
}


mc_Pythia6 = Data('Bplus/t', mc_files['2011']['Pythia6'] + mc_files['2012']['Pythia6'])
mc_Pythia8 = Data('Bplus/t', mc_files['2011']['Pythia8'] + mc_files['2012']['Pythia8'])

mc_all_files = []
for year, pythias in mc_files.items():
    for pythia, files in pythias.items():
        mc_all_files += files

mc_total = Data('Bplus/t', mc_all_files)





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
