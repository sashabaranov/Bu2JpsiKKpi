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
    fixAlphaL=2.1,
    fixAlphaR=1.915,
    fixNL=1.345,
    fixNR=2.69,
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



model_Bu = Models.Fit1D(
    signal=s1_Bu,
    background=Models.Bkg_pdf('BBu', mass=m_Bu, power=0), suffix='Bu'
)

model_Bu_mc = Models.Fit1D(
    signal=s1_Bu,
    # signal2=KKK_pdf,
    # signal3=Kpipi_pdf,
    background=Models.Bkg_pdf('BBu', mass=m_Bu), suffix='Bu'
)

model_Bu_mc.b.fix(0)
model_Bu_mc.background.tau.fix(0)

# model_Bu.b.fix(0)
model_Bu.background.tau.fix(0)
