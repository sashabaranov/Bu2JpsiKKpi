import os

import ROOT
from Ostap.Data import Data, DataAndLumi



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
mc_files = {
    '2011': {
        # 'Pythia6': from_ganga(221) + from_ganga(222),
        'Pythia8': from_ganga(581) + from_ganga(582),
    },

    '2012': {
        # 'Pythia6': from_ganga(223) + from_ganga(224),
        'Pythia8': from_ganga(579) + from_ganga(580),
    }
}

# mc_Pythia6 = Data(branch='Bplus/t', files=mc_files['2011']['Pythia6'] + mc_files['2012']['Pythia6'])
mc_Pythia8 = Data('Bplus/t', mc_files['2011']['Pythia8'] + mc_files['2012']['Pythia8'])

# mc_all_files = []
# for year, pythias in mc_files.items():
#     for pythia, files in pythias.items():
#         mc_all_files += files

# mc_total = Data(branch='Bplus/t', files=mc_all_files)

