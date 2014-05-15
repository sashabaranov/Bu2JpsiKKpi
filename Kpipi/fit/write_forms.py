from tools import *
from data import mc_Pythia6, mc_Pythia8
from fit7 import ds_Bu


db = shelve.open('$KKpidir/fit/histos.shelve')

###############
# Write RD
###############

hists = [
    ("pi1", "mass_pi1ask", "SBu_sw"),
    ("pi1_cuts", "mass_pi1ask", "SBu_sw && ann_pion_K > 0.1"),
]

d = {
    'RD': {param[0]: make_hist(ds_Bu, *param) for param in hists}
}



###############
# Write MC
###############


hists = [
    ["pi1", "mass_pi1ask", ""],
    ["pi1_cuts", "mass_pi1ask", "ann_pion_K[0] > 0.1"],
]

d['MC'] = {}

for param in hists:
    default = param[0]

    param[0] = "p6_" + default
    d['MC']['p6_' + param[0]] = make_hist_mc(mc_Pythia6.data, *param)

    param[0] = "p8_" + default
    d['MC']['p8_' + param[0]] = make_hist_mc(mc_Pythia8.data, *param)



db['Kpipi'] = d

db.sync()
db.close()