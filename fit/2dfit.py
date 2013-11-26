import ROOT

from AnalysisPython.PyRoUts import *
from AnalysisPython.Utils import timing
from AnalysisPython.Utils import rooSilent
from AnalysisPython.Logger import getLogger

from cuts import m_Bu, nbin_Bu, m_Phi
from cuts import cuts_Bu, prntCuts
from cuts import h1

from model import model2_KK, model2_Kpi
from data_Bc import tBu2 as tBu
from Selectors import SBT


logger = getLogger(__name__)
logger.info('DATA chain name is %s ' % (tBu.GetName()))
for i in prntCuts(cuts_Bu, "  CUTS B+  "):
    logger.info(i)


logger.info('Fill control B+  histogram (takes some time)')
with timing():
    tBu.Project(h1.GetName(), 'DTFm_b', cuts_Bu)


sel_Bu = SBT(m_Bu, cuts_Bu)
logger.info('Build RooFit dataset for B+ , it could take as long as 3-5 minutes')
with timing():
    with rooSilent():
        tBu.process(sel_Bu)

ds_Bu = sel_Bu.dataset()
ds_Bu.Print('v')

logger.info('Make unbinned fit for B+')
with rooSilent():
    model2_KK.ss.setMax(1.2 * len(ds_Bu))
    ru = model2_KK.fitTo(ds_Bu)
    # model2_KK.signal.mean .fix(ru('mean_Bu1')[0])
    # model2_KK.signal.sigma.fix(ru('sigma_Bu1')[0])

    ru = model2_KK.fitTo(ds_Bu)



logger.info('running sPlot')
model2_KK.sPlot(ds_Bu)


# xframe = m_Bu.frame(50)
# xframe.SetTitle('J/#psi KK #pi fit projection')
# ds.plotOn(xframe)
# model2_KK.sig1.plotOn(xframe)
# model2_KK.sig1.paramOn(xframe)

# c3 = TCanvas("c3","c3")
# c3.cd()
# xframe.Draw()


ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetOptTitle(1)
ROOT.gStyle.SetOptStat(0)

canvas = ROOT.TCanvas("canvas", "2D fit projections", 1524,768)
canvas.Divide(4, 2)


i = 0
i += 1; canvas.cd(i); fu = model2_KK.ss_pdf.plotOn(m_Bu.frame()); fu.Draw();
i += 1; canvas.cd(i); fu = model2_KK.sb_pdf.plotOn(m_Bu.frame()); fu.Draw();
i += 1; canvas.cd(i); fu = model2_KK.bs_pdf.plotOn(m_Bu.frame()); fu.Draw();
i += 1; canvas.cd(i); fu = model2_KK.bb_pdf.plotOn(m_Bu.frame()); fu.Draw();

i += 1; canvas.cd(i); fu = model2_KK.ss_pdf.plotOn(m_Phi.frame()); fu.Draw();
i += 1; canvas.cd(i); fu = model2_KK.sb_pdf.plotOn(m_Phi.frame()); fu.Draw();
i += 1; canvas.cd(i); fu = model2_KK.bs_pdf.plotOn(m_Phi.frame()); fu.Draw();
i += 1; canvas.cd(i); fu = model2_KK.bb_pdf.plotOn(m_Phi.frame()); fu.Draw();







# pdfs = [model2_KK.ss_pdf, model2_KK.sb_pdf, model2_KK.bs_pdf, model2_KK.bb_pdf]
# frames = [m_Bu.frame(), m_Phi.frame()]

# i = 0
# for pdf in pdfs:
#     i += 1
#     canvas.cd(i)

#     fu = pdf.plotOn(frames[0]); fu.Draw()
#     del fu

# for pdf in pdfs:
#     i += 1
#     canvas.cd(i)

#     fu = pdf.plotOn(frames[1]); fu.Draw()
#     del fu

# canvas.cd(1)
# fu = model2_KK.ss_pdf.plotOn(m_Bu.frame())
# fu.Draw()

# canvas.cd(2)
# fu = model2_KK.ss_pdf.plotOn(m_Phi.frame())
# fu.Draw()

# canvas.cd(3)
# # hhD = ds_Bu.createHistogram(m_Bu, m_Phi, 50, 50)
# # hhD.Draw('lego2')

# h1.SetBins(30, 0.0, 2.0)
# h1.blue()

# ds_Bu.project(h1, "DTFm_kk", "S1S2_sw")
# h1.Draw()


# hh = model2_KK.ss_pdf.createHistogram('DTFm_b:DTFm_KK')
# hh.SetLineColor(63)
# hh.Draw("surf same")

# logger.info('Fit FOM: '+str(ru("SBu")[0].prec()))



# h4 = ROOT.TH1F(hID(), '', 25, m_Bu.getMin(), m_Bu.getMax())
# h4.Sumw2()
# h4.SetXTitle("[pt_K < 0.3]#\,Inv.\,mass(J/\psi K K\pi), GeV/c^2")



# tBu.project(h4, "DTFm_b", cuts_Bu)

# h4.Draw()

  
logger.info('end of the module')