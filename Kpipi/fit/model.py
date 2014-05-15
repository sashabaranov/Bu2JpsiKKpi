import ROOT
from AnalysisPython.PyRoot import *
from AnalysisPython.PyRoUts import *
from AnalysisPython.Utils import timing
from AnalysisPython.Utils import rooSilent
from AnalysisPython.Logger import getLogger

logger = getLogger(__name__)

import PyPAW.FitModels as Models

from variables import m_Bu


s1_Bu = Models.CB2_pdf(
    'Bu1',
    m_Bu.getMin(),
    m_Bu.getMax(),
    fixMass=5.2792e+0,
    fixSigma=0.008499e+0,
    fixAlphaL=1.89,
    fixAlphaR=1.94,
    fixNL=0.8,
    fixNR=2.4,
    mass=m_Bu
)

# s1_Bu = Models.Gauss_pdf(
#     name='Bu1',
#     mn=m_Bu.getMin(),
#     mx=m_Bu.getMax(),
#     fixMass=5.2792e+0,
#     fixSigma=0.008499e+0,
#     mass=m_Bu,
#     mean=None,
#     sigma=None
# )

alist = ROOT.RooArgList(m_Bu)


model_Bu = Models.Fit1D(
    signal=s1_Bu,
    background=Models.Bkg_pdf('BBu', mass=m_Bu), suffix='Bu'
)

# model_Bu.b.fix(0)
# model_Bu.background.tau.fix(0)