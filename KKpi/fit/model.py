from tools import *
import PyPAW.FitModels as Models

from variables import m_Bu
from variables import m_Kstar, mm_Kstar, width_Kstar, low_kpi, high_kpi
from variables import m_Phi, mm_Phi, width_Phi, low_kk, high_kk

from MyFitModels import Charm3_pdf


db = shelve.open("$KKpidir/fit/histos.shelve")

# B+ -> J/psi K+ K- pi+
s1_Bu = Models.CB2_pdf(
    'Bu1',
    m_Bu.getMin(),
    m_Bu.getMax(),
    fixMass=5.2792e+0,
    fixSigma=0.008499e+0,
    fixAlphaL=2.1018e+00,
    fixAlphaR=1.9818e+00,
    fixNL=6.1464e-01,
    fixNR=2.1291e+00,
    mass=m_Bu
)


# s1_Bu = Models.Gauss_pdf(
#     'Bu1',
#     m_Bu.getMin(),
#     m_Bu.getMax(),
#     fixMass=5.2792e+0,
#     fixSigma=0.008499e+0,
#     mass=m_Bu,
# )

kkk_hist = db['KKK']['RD']['k3']
kpipi_hist = db['Kpipi']['MC']['p8_pi1']

kkk = Models.H1D_pdf(name=kkk_hist.GetName(), mass=m_Bu, histo=smear_kkk(kkk_hist))
kpipi = Models.H1D_pdf(name=kpipi_hist.GetName(), mass=m_Bu, histo=smear_kpipi(kpipi_hist))

model_Bu = Charm3_pdf(
    signal=s1_Bu,
    signal2=kkk,
    signal3=kpipi,
    background=Models.Bkg_pdf('BBu', mass=m_Bu), suffix='Bu'
)

model_Bu.background.tau.setMax(-2.0)
model_Bu.background.tau.setVal(-1.0)

model_Bu.s.setMin(100)
model_Bu.s.setVal(101)

model_Bu_mc = Charm3_pdf(
    signal=s1_Bu,
    background=Models.Bkg_pdf('BBu', mass=m_Bu), suffix='Bu'
)

model_Bu_mc.b.fix(0)
model_Bu_mc.background.tau.fix(0)

# model_Bu = Models.Fit1D(
#     signal=s1_Bu,
#     background=Models.Bkg_pdf('BBu', mass=m_Bu, power=4), suffix='Bu'
# )

# model_Bu.background.tau.fix(0)
# model_Bu.b.fix(0)

# anti-Kstar0 -> K- pi+
# s1_Kstar = Models.BreitWigner_pdf(
#     'Kstar1',
#     x=m_Kstar,
#     mass=mm_Kstar,
#     width=width_Kstar,
#     m1=0.493,
#     m2=0.139,
#     L=1,
#     rho=1
# )


# model_Kstar = Models.Fit1D(
#     signal=s1_Kstar,
#     # background = Models.Bkg_pdf ( 'BKstar' , mass = m_Kstar , power = 2 ) ,
#     # suffix = 'Kstar'
#     background=Models.PhaseSpace_pdf(
#         'BKstar', x=m_Kstar, low=low_kpi, high=high_kpi, N=2, L=4, power=0),
#     suffix='Kstar'
# )

# model_Kstar.signal.mean.fix(0.892)
# model_Kstar.signal.width.fix(0.051)


# # phi(1020) -> KK
# s1_Phi = Models.BreitWigner_pdf(
#     'Phi1',
#     x=m_Phi,
#     mass=mm_Phi,
#     width=width_Phi,
#     m1=0.493,
#     m2=0.493,
#     L=1,
#     rho=1
# )


# model_Phi = Models.Fit1D(
#     signal=s1_Phi,
#     #background = Models.Bkg_pdf ( 'BPhi' , mass = m_Phi , power = 4 ) ,
#     background=Models.PhaseSpace_pdf(
#         'BPhi', x=m_Phi, low=low_kk, high=high_kk, N=2, L=4, power=0),
#     suffix='Phi'
# )

# model_Phi.signal.mean.fix(1.020)
# model_Phi.signal.width.fix(4.26e-3)

# model_Phi.signal.mean.release()
# model_Phi.signal.width.release()

