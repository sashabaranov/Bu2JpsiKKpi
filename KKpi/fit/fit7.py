from tools import *

from histograms import h1
from variables import *
from cuts import *
from model import model_Bu
from data import tSelection7 as tBu


# cuts_Bu += reflection1 + reflection2

for i in prntCuts(cuts_Bu, "  CUTS B+  "):
    logger.info(i)



# Fill control B+  histogram (takes some time)
with timing():
    tBu.Project(h1.GetName(), 'DTFm_b', cuts_Bu)


# Fit Bc+ & B+ histogram (check the model)
r, f = model_Bu.fitHisto(h1)


sel_Bu = SelectorWithVars(
    variables=selector_variables,
    selection=cuts_Bu
)

# Build RooFit dataset for B+ , it could take as long as 3-5 minutes
tBu.process(sel_Bu)

ds_Bu = sel_Bu.dataset()
ds_Bu.Print('v')

# frame = m_Bu.frame(nbin_Bu)
# frame.SetXTitle('#Inv.\,mass(J/\psi \, KK\pi), GeV/c^2')
# frame.SetYTitle('Events / (%d \, MeV/c^{2})' % int(events_binning))

# ds_Bu.plotOn(frame)
# frame.Draw()



logger.info('Make unbinned fit for B+')

model_Bu.s.setMax(1.2 * len(ds_Bu))
ru, fu = model_Bu.fitTo(ds_Bu, draw=True, nbins=nbin_Bu)

model_Bu.signal.mean.release()
ru, fu = model_Bu.fitTo(ds_Bu, draw=True, nbins=nbin_Bu)

# model_Bu.signal.aR.release()
# model_Bu.signal.nR.release()
# ru, fu = model_Bu.fitTo(ds_Bu, draw=True, nbins=nbin_Bu)

# model_Bu.signal.aL.release()
# model_Bu.signal.nL.release()
# ru, fu = model_Bu.fitTo(ds_Bu, draw=True, nbins=nbin_Bu)

model_Bu.signal.sigma.release()
ru, fu = model_Bu.fitTo(ds_Bu, draw=True, nbins=nbin_Bu)

# ru, fu = model_Bu.fitTo(ds_Bu, draw=True, nbins=nbin_Bu)


# print 'FIT#2 results for B+  ', ru(model_Bu.s_name)[0]
# print 'FIT#2 precision:', ru("SBu")[0].prec()


# logger.info('running sPlot')
model_Bu.sPlot(ds_Bu)

# def count_significance():
#     global ds_Bu, nbin_Bu, model_Bu
#     from math import sqrt

#     vals = [
#         model_Bu.s,
#         # model_Bu.s2,
#         # model_Bu.s3,
#         model_Bu.b,
#         model_Bu.background.tau,
#         model_Bu.signal.aL,
#         model_Bu.signal.aR,
#         model_Bu.signal.nL,
#         model_Bu.signal.nR,
#         model_Bu.signal.mean,
#         model_Bu.signal.sigma
#     ]

#     for x in vals:
#         x.fix(x.getVal())

#     ru, fu = model_Bu.fitTo(ds_Bu, draw=False, nbins=nbin_Bu)
#     Lfixed = ru.minNll()

#     model_Bu.s.fix(0)

#     ru, fu = model_Bu.fitTo(ds_Bu, draw=False, nbins=nbin_Bu)
#     Ls0 = ru.minNll()
#     return sqrt(2 * (Ls0 - Lfixed))


# logger.info('=' * 20)
# logger.info("Signficance7 is " + str(count_significance()))
# logger.info('=' * 20)

