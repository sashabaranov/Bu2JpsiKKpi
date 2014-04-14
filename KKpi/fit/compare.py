import ROOT
from fit5 import ds_Bu as ds5
from fit6 import ds_Bu as ds6

from Selectors import filter_name

# plotting against MC

ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetOptTitle(1)
ROOT.gStyle.SetOptStat(0)

canvas = ROOT.TCanvas("canvas", "MC vs RD", 700,700)
canvas.Divide(3,3)

from Variables import comparison_variables

hists_rd = []
hists_mc = []


mctrue = " && mcTrueB && mcTruePsi && mcTrueK1 && mcTrueK2 && mcTrueK3 && mcTrueMu1 && mcTrueMu2"

for v in varlist:
    hh_rd = ROOT.TH1F("rd_h_" + v[0], '', 25, v[2], v[3])
    hh_rd.Sumw2()
    ds5.project(hh_rd, filter_name(v[0]), "SBu_sw")

    hh_mc = ROOT.TH1F("mc_h_" + v[0], '', 25, v[2], v[3])
    hh_mc.Sumw2()

    ds6.project(hh_mc, filter_name(v[0]), "SBu_sw")
    # tBu5.project(hh_rd, v[0], cuts_Bu)
    # tBu6.project(hh_mc, v[0], cuts_Bu)

    hh_rd.scale(1.0)
    hh_mc.scale(1.0)

    hh_rd.SetXTitle(v[-1])
    hh_mc.SetXTitle(v[-1])

    hh_rd.SetMarkerStyle(21)
    hh_rd.SetMarkerSize(0.5)

    hh_mc.SetMarkerStyle(21)
    hh_mc.SetMarkerSize(0.5)

    hists_rd.append(hh_rd)
    hists_mc.append(hh_mc)


for i, name in enumerate(varlist):
    canvas.cd(i+1)
    hists_rd[i].red()
    hists_mc[i].blue()

    hists_mc[i].Draw()
    hists_rd[i].Draw("same")


canvas.Update()

# logger.info('end of the module')
