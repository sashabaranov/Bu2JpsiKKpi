from tools import *

from variables import *
from cuts import cuts_Bu, prntCuts
from model import model_Bu
from data import tSelection7 as tBu


logger.info('DATA chain name is %s ' % (tBu.GetName()))

for i in prntCuts(cuts_Bu, "  CUTS B+  "):
    logger.info(i)


# Fit histo
# logger.info('Fill control B+  histogram (takes some time)')
# with timing():
#     tBu.Project(h1.GetName(), 'DTFm_b', cuts_Bu)

# with rooSilent():
#     logger.info('Fit Bc+ & B+ histogram (check the model)')
#     r, f = model_Bu.fitHisto(h1)




sel_Bu = SelectorWithVars(
    variables=selector_variables,
    selection=cuts_Bu
)


logger.info('Build RooFit dataset for B+ , it could take as long as 3-5 minutes')

tBu.process(sel_Bu)

ds_Bu = sel_Bu.dataset()
ds_Bu.Print('v')




logger.info('Make unbinned fit for B+')
#
model_Bu.s.setMax(1.2 * len(ds_Bu))
ru, fu = model_Bu.fitTo(ds_Bu, draw=True, nbins=nbin_Bu)

model_Bu.signal.sigma.release()
model_Bu.signal.mean.release()
ru, fu = model_Bu.fitTo(ds_Bu, draw=True, nbins=nbin_Bu)

# model_Bu.signal.aR.release()
# model_Bu.signal.aL.release()
# ru, fu = model_Bu.fitTo(ds_Bu, draw=True, nbins=nbin_Bu)

# model_Bu.signal.nR.release()
# model_Bu.signal.nL.release()

# ru, fu = model_Bu.fitTo(ds_Bu, draw=True, nbins=nbin_Bu)



# print 'FIT#2 results for B+  ', ru(model_Bu.s_name)[0]
# print 'FIT#2 precision:', ru("SBu")[0].prec()



fu.SetXTitle('#Inv.\,mass(J/\psi\,KKK), GeV/c^2')
fu.SetYTitle('Events / (%d \, MeV/c^{2})' % events_binning)

fu.Draw()




logger.info('running sPlot')
model_Bu.sPlot(ds_Bu)



# ===============================================
# Write histos
# ===============================================



hists = [
    ("k1", "mass_k1aspi", "SBu_sw"),
    ("k3", "mass_k3aspi", "SBu_sw"),
    ("k1_cuts", "mass_k1aspi", "SBu_sw && ann_kaon_PI_2 > 0.1"),
    ("k3_cuts", "mass_k3aspi", "SBu_sw && ann_kaon_PI_0 > 0.1"),
    ("k1_cuts_strong", "mass_k1aspi", "SBu_sw && ann_kaon_PI_2 > 0.1 && ann_kaon0 > 0.3"),
    ("k3_cuts_strong", "mass_k3aspi", "SBu_sw && ann_kaon_PI_0 > 0.1 && ann_kaon2 > 0.3"),
]


db = shelve.open('$KKpidir/fit/histos.shelve')

db['KKK'] = {
    'RD': {param[0]: make_hist(ds_Bu, *param) for param in hists},
    'MC': {}
}



# ===============================================
# Drawing
# ===============================================


h1, h2 = db['KKK']['RD']['k1'], db['KKK']['RD']['k3']
#h1, h2 = db['KKK']['RD']['k1_cuts'], db['KKK']['RD']['k3_cuts']


title = '#Inv.\,mass(J/\psi\,K^+\,K^-\,K^+) \,\, with \,\, misid, GeV/c^2'
h1.SetXTitle(title)
h2.SetXTitle(title)

h1.red()
h2.blue()

h1.Draw()
h2.Draw('same')

db.close()
