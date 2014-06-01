import ROOT
import PyPAW.FitModels as Models
from variables import m_R

s1 = Models.Gauss_pdf(
    'psi2s',
    3.6,
    3.75,
    fixMass=3.686,
    fixSigma=0.005,
    mass=m_R
)

s2 = Models.Gauss_pdf(
    'X',
    3.8,
    3.9,
    fixMass=3.872,
    fixSigma=0.002,
    mass=m_R
)

model_Resonant = Models.Fit1D(
    signal=s1,
    background=Models.Bkg_pdf('BRes', mass=m_R, power=1),
    components=[s2],
    suffix='Res'
)