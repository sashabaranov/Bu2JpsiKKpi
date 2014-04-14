import ROOT
from AnalysisPython.PyRoUts import *
from AnalysisPython.Utils import timing
from AnalysisPython.Utils import rooSilent
from AnalysisPython.Logger import getLogger

logger = getLogger(__name__)

import DiCharm.Analysis.Models as Models

from cutsKKK import m_Bu, h1

# B+ -> J/psi KKK
s1_Bu = Models.CB2_pdf(
    'Bu1',
    m_Bu.getMin(),
    m_Bu.getMax(),
    fixMass=5.16e+0,
    fixSigma=0.05e+0,
    fixAlphaL=0.45575,
    fixAlphaR=0.446336,
    fixNL=4.999,
    fixNR=9.92,
    mass=m_Bu)


model_Bu = Models.Charm1_pdf( # Charm1_pdf(
    signal=s1_Bu,
    #hist1=h1,
    suffix='Bu'
    # background = Models.Bkg_pdf ( 'BBu' , mass = m_Bu , power = 2 )
)


model_Bu.signal.    aL     .  release()
model_Bu.signal.    aR     .  release()
model_Bu.signal.  mean     .  release()
model_Bu.signal.    nL     .  release()
model_Bu.signal.    nR     .  release()
model_Bu.signal. sigma     .  release()


model_Bu.signal.    aL     . setVal(1.5653e-07)
model_Bu.signal.    aR     . setVal(1.9342e+00)
model_Bu.signal.  mean     . setVal(5.1531e+00)
model_Bu.signal.    nL     . setVal(2.9020e+00)
model_Bu.signal.    nR     . setVal(1.3168e-06)
model_Bu.signal. sigma     . setVal(4.0348e-02)

model_Bu.b.fix(0.0)
model_Bu.background.tau.fix(0.0)
