from tools import *
from data import mc_Pythia6, mc_Pythia8
from fit7 import ds_Bu


db = shelve.open('$KKpidir/fit/histos.shelve')

###############
# Write RD
###############

hists = [
    ("k1", "mass_k1aspi", "SBu_sw"),
    ("k3", "mass_k3aspi", "SBu_sw"),
    ("k1_cuts", "mass_k1aspi", "SBu_sw && ann_kaon_PI_2 > 0.1"),
    ("k3_cuts", "mass_k3aspi", "SBu_sw && ann_kaon_PI_0 > 0.1"),
]

d = {
    'RD': {param[0]: make_hist(ds_Bu, *param) for param in hists}
}



###############
# Write MC
###############


hists = [
    ["k1", "mass_k1aspi", ""],
    ["k3", "mass_k3aspi", ""],
    ["k1_cuts", "mass_k1aspi", "ann_kaon_PI[2] > 0.1"],
    ["k3_cuts", "mass_k3aspi", "ann_kaon_PI[0] > 0.1"],
]


d['MC'] = {}

for param in hists:
    default = param[0]

    param[0] = "p6_" + default
    d['MC'][param[0]] = make_hist_mc(mc_Pythia6.data, *param)

    param[0] = "p8_" + default
    d['MC'][param[0]] = make_hist_mc(mc_Pythia8.data, *param)



db['KKK'] = d

db.sync()
db.close()