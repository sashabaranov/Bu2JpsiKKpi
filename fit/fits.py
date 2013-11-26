import ROOT

from AnalysisPython.PyRoUts import *
from AnalysisPython.Utils import timing
from AnalysisPython.Utils import rooSilent
from AnalysisPython.Logger import getLogger

from cuts import m_Bu, nbin_Bu
from cuts import cuts_Bu, prntCuts
from cuts import h1, h2, h3

from model import model_Bu
from data_Bc import tBu2 as tBu
from Selectors import SBT


logger = getLogger(__name__)
logger.info('DATA chain name is %s ' % (tBu.GetName()))
for i in prntCuts(cuts_Bu, "  CUTS B+  "):
    logger.info(i)


logger.info('Fill control B+  histogram (takes some time)')
with timing():
    tBu.Project(h1.GetName(), 'DTFm_b', cuts_Bu)


logger.info('Fit Bc+ & B+ histogram (check the model)')
with rooSilent():
    r, f = model_Bu.fitHisto(h1)
    # model_Bu.signal.mean .release()



sel_Bu = SBT(m_Bu, cuts_Bu)
logger.info('Build RooFit dataset for B+ , it could take as long as 3-5 minutes')
with timing():
    with rooSilent():
        tBu.process(sel_Bu)

ds_Bu = sel_Bu.dataset()
ds_Bu.Print('v')


logger.info('Make unbinned fit for B+')
with rooSilent():
    model_Bu.s.setMax(1.2 * len(ds_Bu))
    ru, fu = model_Bu.fitTo(ds_Bu, draw=False, nbins=nbin_Bu)
    
    # model_Bu.signal.mean .fix(ru('mean_Bu1')[0])
    # model_Bu.signal.sigma.fix(ru('sigma_Bu1')[0])
    model_Bu.signal.sigma.release()
    model_Bu.signal.mean.release()

    model_Bu.signal.aR.release()
    model_Bu.signal.aL.release()
    
    model_Bu.signal.nR.release()
    model_Bu.signal.nL.release()

    ru, fu = model_Bu.fitTo(ds_Bu, draw=True, nbins=nbin_Bu)



logger.info('running sPlot')
model_Bu.sPlot(ds_Bu)


# w_s = "SBu_sw"
# w_b = "BBu_sw"
# w_r = "S3Bu_sw"

# h1.SetBins(30, -0.0, 1.0)
# h2.SetBins(30, -0.0, 1.0)
# h3.SetBins(30, -0.0, 1.0)

# h1.blue()
# h2.red()
# h3.SetLineColor(ROOT.kMagenta)

# ds_Bu.project(h1, "ann_kaon*(1-ann_pion)", w_s)
# ds_Bu.project(h2, "ann_kaon*(1-ann_pion)", w_b)
# ds_Bu.project(h3, "ann_kaon*(1-ann_pion)", w_r)


# h1.SetXTitle("PROBNNk*(1-PROBNN\pi)")
# h2.SetXTitle("PROBNNk*(1-PROBNN\pi)")
# h3.SetXTitle("PROBNNk*(1-PROBNN\pi)")

# h2.Draw()
# h3.Draw("SAME")
# h1.Draw("SAME")

# logger.info('Fit FOM: '+str(ru("SBu")[0].prec()))



# h4 = ROOT.TH1F(hID(), '', 25, m_Bu.getMin(), m_Bu.getMax())
# h4.Sumw2()
# h4.SetXTitle("[pt_K < 0.3]#\,Inv.\,mass(J/\psi K K\pi), GeV/c^2")

# tBu.project(h4, "DTFm_b", "S1S2_sw")

# h4.Draw()

  
logger.info('end of the module')
