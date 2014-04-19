import ROOT
from AnalysisPython.PyRoot import *
from AnalysisPython.PyRoUts import *
from AnalysisPython.Utils import timing
from AnalysisPython.Utils import rooSilent
from AnalysisPython.Logger import getLogger

logger = getLogger(__name__)

import DiCharm.Analysis.Models as Models

from cuts import m_Bu


# B+ -> J/psi K+ K- K+
s1_Bu = Models.CB2_pdf(
    'Bu1',
    m_Bu.getMin(),
    m_Bu.getMax(),
    fixMass=5.2792e+0,
    fixSigma=0.008499e+0,
    fixAlphaL=0.45575,
    fixAlphaR=0.446336,
    fixNL=4.999,
    fixNR=9.92,
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


model_Bu = Models.Charm3_pdf(
    signal=s1_Bu,
    # signal2=KKK_pdf,
    # signal3=Kpipi_pdf,
    background=Models.Bkg_pdf('BBu', mass=m_Bu), suffix='Bu'
)
