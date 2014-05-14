from tools import *

from cuts import cuts_Bu, prntCuts, mctrue

from model import model_Bu_mc as model_Bu
from data import mc_Pythia6, mc_Pythia8, mc_total


cuts_Bu += mctrue


# logger.info('DATA chain name is %s ' % (tBu.GetName()))


for i in prntCuts(cuts_Bu, "  CUTS B+  "):
    logger.info(i)


# logger.info('Fill control B+  histogram (takes some time)')
# with timing():
#     tBu.Project(h1.GetName(), 'DTFm_b', cuts_Bu)



# with rooSilent():
#     logger.info('Fit Bc+ & B+ histogram (check the model)')
#     r, f = model_Bu.fitHisto(h1)



# from variables import *

# sel_Bu = SelectorWithVars(
#     variables=selector_variables,
#     selection=cuts_Bu
# )

# logger.info('Build RooFit dataset for B+ , it could take as long as 3-5 minutes')

# tBu.process(sel_Bu)

# ds_Bu = sel_Bu.dataset()
# ds_Bu.Print('v')


# logger.info('Make unbinned fit for B+')
# #
# model_Bu.s.setMax(1.2 * len(ds_Bu))
# ru, fu = model_Bu.fitTo(ds_Bu, draw=True, nbins=nbin_Bu)

# model_Bu.signal.sigma.release()
# model_Bu.signal.mean.release()
# ru, fu = model_Bu.fitTo(ds_Bu, draw=True, nbins=nbin_Bu)

# model_Bu.signal.aR.release()
# model_Bu.signal.aL.release()
# ru, fu = model_Bu.fitTo(ds_Bu, draw=True, nbins=nbin_Bu)

# model_Bu.signal.nR.release()
# model_Bu.signal.nL.release()

# ru, fu = model_Bu.fitTo(ds_Bu, draw=True, nbins=nbin_Bu)


# fu.SetXTitle('#Inv.\,mass(J/\psi\,KKK), GeV/c^2')
# fu.SetYTitle('Events / (%d \, MeV/c^{2})' % events_binning)

# fu.Draw()


# logger.info('running sPlot')
# model_Bu.sPlot(ds_Bu)


# ===============================================
# Write histos
# ===============================================


hists = [
    ["k1", "mass_k1aspi", ""],
    ["k3", "mass_k3aspi", ""],
    ["k1_cuts", "mass_k1aspi", "ann_kaon_PI[2] > 0.1"],
    ["k3_cuts", "mass_k3aspi", "ann_kaon_PI[0] > 0.1"],
]


logger.info('Writing histos')
db = shelve.open('$KKpidir/fit/histos.shelve')

d = db['KKK']

for param in hists:
    default = param[0]

    param[0] = "p6_" + default
    d['MC']['p6_' + param[0]] = make_hist_mc(mc_Pythia6.data, *param)

    param[0] = "p8_" + default
    d['MC']['p8_' + param[0]] = make_hist_mc(mc_Pythia8.data, *param)

db['KKK'] = d




# ===============================================
# Drawing
# ===============================================

histos = sorted([h for h in d['MC'].values() if 'p8' in h.GetName()], key=lambda h: -h.GetMaximum())

title = '#Inv.\,mass(J/\psi\,K^+\,K^-\,K^+) \,\, with \,\, misid, GeV/c^2'


legend = ROOT.TLegend(0.55, 0.65, 0.86, 0.9);
legend.SetFillColor(ROOT.kWhite)



for i, h in enumerate(histos):
    col = i + 1

    h.SetLineColor(col)
    h.SetFillColor(col)
    h.SetMarkerColor(col)

    h.SetXTitle(title)
    h.Draw('same')

    legend.AddEntry(h.GetName(), make_legend(h.GetName()), "P")

legend.Draw()

db.sync()
db.close()

