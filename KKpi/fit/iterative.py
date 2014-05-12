import ROOT
from AnalysisPython.PyRoot import *
from PyPAW.PyRoUts import *
from PyPAW.Utils import timing
from PyPAW.Utils import RooSilent
from AnalysisPython.Logger import getLogger
from PyPAW.Selectors import SelectorWithVars



from histograms import h1
from variables import m_Bu, nbin_Bu, selector_variables
from cuts import cuts_Bu, prntCuts
from data import tSelection7 as tBu

from math import sqrt

logger = getLogger(__name__)

# sel_Bu = SelectorWithVars(
#     variables=selector_variables,
#     selection=cuts_Bu
# )

# # Build RooFit dataset for B+ , it could take as long as 3-5 minutes
# tBu.process(sel_Bu)

# ds_Bu = sel_Bu.dataset()
# ds_Bu.Print('v')


# logger.info('Make unbinned fit for B+')

from model import db, smear, Models

all_KKK = db['KKK']['RD'].values() + db['KKK']['MC'].values()
all_Kpipi = db['Kpipi']['RD'].values() + db['Kpipi']['MC'].values()

def count_significance(model_Bu, ds_Bu, nbin_Bu):
    with RooSilent():
        try:
            model_Bu.s.setMax(1.2 * len(ds_Bu))

            ru, fu = model_Bu.fitTo(ds_Bu, draw=False, nbins=nbin_Bu)

            model_Bu.signal.sigma.release()
            model_Bu.signal.mean.release()


            ru, fu = model_Bu.fitTo(ds_Bu, draw=False, nbins=nbin_Bu)
            Sfixed = str(ru("SBu")[0])

            vals = [
                model_Bu.background.tau,
                model_Bu.signal.mean,
                model_Bu.signal.sigma
            ] + model_Bu.nums

            for x in vals:
                x.fix(x.getVal())

            ru, fu = model_Bu.fitTo(ds_Bu, draw=False, nbins=nbin_Bu)
            Lfixed = ru.minNll()

            model_Bu.s.fix(0)

            ru, fu = model_Bu.fitTo(ds_Bu, draw=False, nbins=nbin_Bu)
            Ls0 = ru.minNll()
            return {'significance': sqrt(2 * (Ls0 - Lfixed)), 'Nsig': Sfixed}
        except ValueError:
            return -1

print len(all_KKK) + len(all_Kpipi)

# f = open('results.txt', 'w')
# import json

# for kkk_hist in all_KKK:
#     for kpipi_hist in all_Kpipi:
#         s1_Bu = Models.CB2_pdf(
#             'Bu1',
#             m_Bu.getMin(),
#             m_Bu.getMax(),
#             fixMass=5.2797e+00,
#             fixSigma=5.2149e-03,
#             fixAlphaL=2.1,
#             fixAlphaR=1.9,
#             fixNL=1.2,
#             fixNR=2.8,
#             mass=m_Bu
#         )

#         kkk = Models.H1D_pdf(name="KKK", mass=m_Bu, histo=smear(kkk_hist))
#         kpipi = Models.H1D_pdf(name="Kpipi", mass=m_Bu, histo=smear(kpipi_hist))

#         model = Models.Fit1D(
#             signal=s1_Bu,
#             components=[kkk, kpipi],
#             background=Models.Bkg_pdf('BBu', mass=m_Bu), suffix='Bu'
#         )
#         result = "{}\t{}\t{}\n".format(
#             kkk_hist.GetName(),
#             kpipi_hist.GetName(),
#             json.dumps(count_significance(model, ds_Bu, nbin_Bu))
#         )

#         f.write(result)

#         print result

# f.close()
# db.close()
