from tools import *

from variables import *
from cuts import cuts_Bu, prntCuts
from model import model_Bu
from data import selection7


logger.info('DATA chain name is %s ' % (selection7.chain.GetName()))


for i in prntCuts(cuts_Bu, "  CUTS B+  "):
    logger.info(i)


# logger.info('Fill control B+ histogram (takes some time)')
# with timing():
#     selection7.chain.Project(h1.GetName(), 'DTFm_b', cuts_Bu)


from Ostap.Selectors import SelectorWithVarsCached

sel_Bu = SelectorWithVarsCached(
    variables=selector_variables,
    selection=cuts_Bu,
    files=selection7.files
)

logger.info('Build RooFit dataset for B+ , it could take as long as 3-5 minutes')

if not sel_Bu._loaded_from_cache:
    selection7.chain.process(sel_Bu)


ds_Bu = sel_Bu.dataset()
ds_Bu.Print('v')


logger.info('Make unbinned fit for B+')
#
model_Bu.s.setMax(1.2 * len(ds_Bu))
ru, fu = model_Bu.fitTo(ds_Bu, draw=True, nbins=nbin_Bu)

model_Bu.signal.sigma.release()
model_Bu.signal.mean.release()

ru, fu = model_Bu.fitTo(ds_Bu, draw=True, nbins=nbin_Bu)


fu.SetXTitle('#Inv.\,mass(J/\psi\,K\pi\pi), GeV/c^2')
fu.SetYTitle('Events / (%d MeV/c^2)' % 4)

fu.Draw()


print ru
# print 'FIT results for B+  ', ru(model_Bu.s_name)[0]
# print 'FIT precision:', ru("SBu")[0].prec()

# print 'FIT#2 results for B+  ', ru(model_Bu.s_name)[0]
# print 'FIT#2 precision:', ru("SBu")[0].prec()


logger.info('running sPlot')
model_Bu.sPlot(ds_Bu)
