import ROOT
from AnalysisPython.PyRoot import *
from AnalysisPython.PyRoUts import *
from AnalysisPython.Utils import timing
from AnalysisPython.Utils import rooSilent
from AnalysisPython.Logger import getLogger

logger = getLogger(__name__)

import PyPAW.FitModels as Models

from variables import m_Bu
from variables import m_Kstar, mm_Kstar, width_Kstar, low_kpi, high_kpi
from variables import m_Phi, mm_Phi, width_Phi, low_kk, high_kk


db = shelve.open("$KKpidir/fit/histos.shelve")

# B+ -> J/psi K+ K- pi+
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
#     'Bu1',
#     m_Bu.getMin(),
#     m_Bu.getMax(),
#     fixMass=5.2792e+0,
#     fixSigma=0.008499e+0,
#     mass=m_Bu,
# )


alist = ROOT.RooArgList(m_Bu)

KKK_h = db['KKK']['RD']['k1_cuts']

KKK_h2 = ROOT.TH1F('KKK_h2', '', 1000, *KKK_h.xminmax())
KKK_h2 += KKK_h

KKK_h3 = ROOT.TH1F('KKK_h3', '', 100, *KKK_h.xminmax())
KKK_h3 += KKK_h2

# KKK_h.smear(0.003)


Kpipi_h = db['Kpipi']['RD']['pi1_cuts']

Kpipi_h2 = ROOT.TH1F('Kpipi_h2', '', 1000, *Kpipi_h.xminmax())
Kpipi_h2 += Kpipi_h

# # Kpipi_h3 = ROOT.TH1F('Kpipi_h3', '', 100, *Kpipi_h.xminmax())
# # Kpipi_h3 += Kpipi_h
# #Kpipi_h.smear(0.003)

KKK_hist   = ROOT.RooDataHist("KKK_hist", "KKK_hist", alist, KKK_h2)
Kpipi_hist = ROOT.RooDataHist("Kpipi_hist", "Kpipi_hist", alist, Kpipi_h2)


KKK_pdf = ROOT.RooHistPdf("KKK_pdf", "Kpipi_pdf", alist, alist, KKK_hist)
Kpipi_pdf = ROOT.RooHistPdf("Kpipi_pdf", "Kpipi_pdf", alist, alist, Kpipi_hist)

kkk = Models.H1D_pdf(name="KKK", mass=m_Bu, histo=db['KKK']['RD']['k1'])
kpipi = Models.H1D_pdf(name="Kpipi", mass=m_Bu, histo=db['Kpipi']['RD']['pi1'])

model_Bu = Models.Fit1D(
    signal=s1_Bu,
    components=[kkk, kpipi],
    background=Models.Bkg_pdf('BBu', mass=m_Bu), suffix='Bu'
)

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


# phi(1020) -> KK
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

# # model_Phi.signal.mean.release()
# # model_Phi.signal.width.release()



# model2_Kpi =  Models.Charm2_pdf(
#     sig_1=s1_Bu,
#     sig_2=s1_Kstar,
#     suffix='2D-Kpi'
# )

# model2_KK =  Models.Charm2_pdf(
#     sig_1=s1_Bu,
#     sig_2=s1_Phi,
#     suffix='KK'
# )
