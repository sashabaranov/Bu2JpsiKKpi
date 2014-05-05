import ROOT
import glob
from AnalysisPython.PyRoUts import *
from AnalysisPython.GetLumi import getLumi
import AnalysisPython.ZipShelve as ZipShelve
from AnalysisPython.Logger import getLogger
from ostap.data import Data, DataAndLumi

# =============================================================================
def append_dir(directory, files):
    return [os.path.join(directory, f) for f in files]

def from_ganga(job_id):
    from glob import glob
    prefix = "/afs/cern.ch/user/a/albarano/gangadir/workspace/albarano/LocalXML/"
    pattern = prefix + str(job_id) + "/*/output/*.root"
    return glob(pattern)

# =============================================================================
mc_files = {
    '2011': {
        'Pythia6': from_ganga(202) + from_ganga(203),
        'Pythia8': from_ganga(200) + from_ganga(201),
    },

    '2012': {
        'Pythia6': from_ganga(204) + from_ganga(205),
        'Pythia8': from_ganga(198) + from_ganga(199),
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
files7 = from_ganga(83) + from_ganga(84)

selection7 = DataAndLumi(branch='JpsiKKpi/t', files=files7)

logger.info("Selection 7: %s" % selection7)

tSelection7 = selection7.data
