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
    fixAlphaL=1.935,
    fixAlphaR=1.84,
    fixNL=0.695,
    fixNR=3.345,
    mass=m_Bu
)

model_Bu = Models.Fit1D(
    signal=s1_Bu,
    background=Models.Bkg_pdf('BBu', mass=m_Bu, power=1), 
    suffix='Bu'
)

# model_Bu.b.fix(0)
# model_Bu.background.tau.fix(0)