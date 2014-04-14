import ROOT
from AnalysisPython.PyRoot import *
from AnalysisPython.PyRoUts import *
from AnalysisPython.Utils import timing
from AnalysisPython.Utils import rooSilent
from AnalysisPython.Logger import getLogger
logger = getLogger(__name__)


from cuts import m_Bu, nbin_Bu
from cuts import cuts_Bu, prntCuts
from cuts import h1, h_mc


logger.info('Import models/PDFs')
from model import model_Bu


logger.info('Import data (for the first time it could take some time)')
from data_Bc import tSelection5 as tBu
from data_Bc import tBu_mc, tBu_mc_p6, tBu_mc_p8


logger.info('DATA chain name is %s ' % (tBu.GetName()))


data_cuts_Bu = cuts_Bu

for i in prntCuts(data_cuts_Bu, "  CUTS B+  "):
    logger.info(i)


logger.info('Fill control B+  histogram (takes some time)')
with timing():
    tBu.Project(h1.GetName(), 'DTFm_b', data_cuts_Bu)


with rooSilent():
    logger.info('Fit Bc+ & B+ histogram (check the model)')
    r, f = model_Bu.fitHisto(h1)

from Selectors import SBT

sel_Bu = SBT(m_Bu, cuts_Bu)

# tBu_mc.project(h_mc, "DTFm_b")
# r, f = model_Bu.fitHisto(h_mc, draw=True)
# h_mc.Draw()

logger.info(
    'Build RooFit dataset for B+ , it could take as long as 3-5 minutes')
with timing():
    ## logger.warning('almost skip B+')
    ## with rooSilent () : tBu.process ( sel_Bu , 100000 )
    with rooSilent():
        tBu.process(sel_Bu)

ds_Bu = sel_Bu.dataset()
ds_Bu.Print('v')

with rooSilent():
    logger.info('Make unbinned fit for B+')
    #
    model_Bu.s.setMax(1.2 * len(ds_Bu))
    ru, fu = model_Bu.fitTo(ds_Bu, draw=True, nbins=nbin_Bu)

    model_Bu.signal.sigma.release()
    model_Bu.signal.mean.release()

    ru, fu = model_Bu.fitTo(ds_Bu, draw=True, nbins=nbin_Bu)

    model_Bu.signal.aL.release()
    model_Bu.signal.aR.release()
    ru, fu = model_Bu.fitTo(ds_Bu, draw=True, nbins=nbin_Bu)


    model_Bu.signal.nL.release()
    model_Bu.signal.nR.release()
    ru, fu = model_Bu.fitTo(ds_Bu, draw=True, nbins=nbin_Bu)

print ru
print 'FIT results for B+  ', ru(model_Bu.s_name)[0]
print 'FIT precision:', ru("SBu")[0].prec()


# print 'FIT#2 results for B+  ', ru(model_Bu.s_name)[0]
# print 'FIT#2 precision:', ru("SBu")[0].prec()


logger.info('running sPlot')
model_Bu.sPlot(ds_Bu)

h = ROOT.TH1F("h", '',  30, 5.16, 5.45)
h.Sumw2()

ds_Bu.project(h, "m_b_misid", "SBu_sw")
h.Draw()

d = shelve.open('$KKpidir/fit/histos.shelve')

d['Kpipi_30bin_sel5'] = h


d.close()

# # plotting against MC

# h = ROOT.TH1F("h", '',  50, 5.0, 5.5)
# h.Sumw2()

# ds_Bu.project(h, "m_b_misid", "SBu_sw")
# h.Draw()



# d = shelve.open('new_result.shelve')

# d['Kpipi_hist'] = h

# d.close()


# ROOT.gStyle.SetPalette(1)
# ROOT.gStyle.SetOptTitle(1)
# ROOT.gStyle.SetOptStat(0)

# canvas = ROOT.TCanvas("canvas", "MC vs RD", 700,700)
# canvas.Divide(3,2)


# # variables = [
# #     # (name, binning, range_min, range_max)
# #     # ("pt_kaon", 25, -0.5, 6.0),
# #     # ("eta_kaon", 25, -0.5, 6.0),
# #     # ("pt_pion", 25, -0.5, 6.0),
# #     # ("eta_pion", 25, -0.5, 6.0)

# #     # ("y_jpsi", 25, 0.0, 5.0),
# #     # ("eta_jpsi", 25, -0.5, 10.0),
# #     # ("pt_jpsi", 25, -0.5, 20.0),
# #     # ("y_b", 25, 0.0, 5.0),
# #     # ("eta_b", 25, -0.5, 10.0),
# #     # ("pt_b", 25, -0.5, 40.0),
    
# #     # ('c2ip_b', 25, 0, 20),
# #     # ('m_jpsi', 25, 3.0, 3.2),
# #     # ('minann_K', 25, 0, 1),
# #     # ('minann_mu', 25, 0, 1),
# #     # ('minPt_track', 25, 0, 3),
# #     # ("pt_kaon", 25, -0.5, 6.0)
# #     ("pt_pion", 25, 0.0, 6.0),
# #     ("p_pion", 25, 0.0, 30.0),
# #     ("eta_pion", 25, 1.0, 6.0)

# #     # ('pt_jpsi', 25, 0, 25),
# # ]

# # variables = [
# #     # (name, binning, range_min, range_max)
# #     ('c2ip_b', 25, 0, 10, "\chi^{2}_{IP}(B)"),
# #     ("eta_b", 25, 2.0, 5.5, "\eta_{B}"),
# #     ("y_b", 25, 2.0, 5.0, "y_{B}"),

# #     ("pt_b", 25, -0.5, 20.0, "p_{T}(B)"),
# #     ('minann_K', 25, 0.96, 1, "min(Prob.NN(K))"),
# #     ('minPt_track', 25, 0, 1.5, "min(p_{T}(track))"),
    
# #     ("pt_pion", 25, 0.0, 2.0, "p_{T}(\pi)"),
# #     ("p_pion", 25, 0.0, 30.0, "p(\pi)"),
# #     ("eta_pion", 25, 1.5, 5.0, "\eta_{\pi}")
# # ]


# variables = [
#     # (name, binning, range_min, range_max)
#     ("eta_b", 25, 1.0, 5.5, "\eta_{B}"),
#     ("y_b", 25, 2.0, 5.0, "y_{B}"),
#     ("pt_b", 25, -0.5, 40.0, "p_{T}(B)"),

#     ('minPt_track', 25, 0, 2, "min(p_{T}(track))"),
#     ('y_jpsi', 25, 1.5, 5.0, "y_{J/\psi}"),
#     ('pt_jpsi', 25, 0, 25, "p_{T}(J/\psi)"),

# ]

# hists_rd = []
# hists_mc = []

# mctrue = " && mcTrueB && mcTruePsi && mcTruePi1 && mcTruePi2 && mcTrueK && mcTrueMu1 && mcTrueMu2"

# for v in variables:
#     hh_rd = ROOT.TH1F("rd_h_" + v[0], '', v[1], v[2], v[3])
#     hh_rd.Sumw2()

#     # rd_name = ''
#     # if not '[0]' in v[0]:
#     #     rd_name = v[0] + '[0]'
#     # else:
#     #     rd_name = v[0]

#     # ds_Bu.project(hh_rd, rd_name, "SBu_sw")
#     tBu_mc_p6.project(hh_rd, v[0], cuts_Bu + mctrue)


#     hh_mc = ROOT.TH1F("mc_h_" + v[0], '', v[1], v[2], v[3])
#     hh_mc.Sumw2()
#     # tBu_mc.project(hh_mc, v[0] + "[0]", cuts_Bu + mctrue)
#     tBu_mc_p8.project(hh_mc, v[0], cuts_Bu + mctrue)

#     hh_rd.scale(1.0)
#     hh_mc.scale(1.0)

#     hh_rd.SetXTitle(v[-1])
#     hh_mc.SetXTitle(v[-1])

#     hh_rd.SetMarkerStyle(21)
#     hh_rd.SetMarkerSize(0.5)

#     hh_mc.SetMarkerStyle(21)
#     hh_mc.SetMarkerSize(0.5)

#     hists_rd.append(hh_rd)
#     hists_mc.append(hh_mc)


# for i, name in enumerate(variables):
#     canvas.cd(i+1)
#     hists_rd[i].red()
#     hists_mc[i].blue()


#     hists_mc[i].Draw()
#     hists_rd[i].Draw("same")
#     # else:
#     #     hists_rd[i].Draw()
#     #     hists_mc[i].Draw("same")


# canvas.Update()

  
# logger.info('end of the module')
