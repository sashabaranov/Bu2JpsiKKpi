import ROOT
import PyPAW.FitModels as Models


from variables import m_Bu


# B+ -> J/psi K+ K- K+
s1_Bu = Models.CB2_pdf(
    'Bu1',
    m_Bu.getMin(),
    m_Bu.getMax(),
    fixMass=5.2792e+0,
    fixSigma=0.008499e+0,
    fixAlphaL=2.1005,
    fixAlphaR=1.8805,
    fixNL=1.215,
    fixNR=2.765,
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
    # signal2=KKK_pdf,
    # signal3=Kpipi_pdf,
    background=Models.Bkg_pdf('BBu', mass=m_Bu, power=2), suffix='Bu'
)

model_Bu_mc = Models.Fit1D(
    signal=s1_Bu,
    # signal2=KKK_pdf,
    # signal3=Kpipi_pdf,
    background=Models.Bkg_pdf('BBu', mass=m_Bu), suffix='Bu'
)

model_Bu_mc.b.fix(0)
model_Bu_mc.background.tau.fix(0)
