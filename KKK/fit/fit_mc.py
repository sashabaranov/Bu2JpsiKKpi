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
