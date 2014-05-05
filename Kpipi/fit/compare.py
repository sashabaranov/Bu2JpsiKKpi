import ROOT
from AnalysisPython.PyRoot import *
from AnalysisPython.PyRoUts import *


# plotting against MC

db = shelve.open('$KKpidir/fit/histos.shelve')




# ROOT.gStyle.SetPalette(1)
# ROOT.gStyle.SetOptTitle(1)
# ROOT.gStyle.SetOptStat(0)

# canvas = ROOT.TCanvas("canvas", "MC vs RD", 1024, 768)

name = "pi1_cuts"

h_rd = db['Kpipi']['RD'][name]
h_mc = db['Kpipi']['MC'][name]
print h_rd, h_mc

h_rd.red()
h_mc.blue()

h_rd.Scale(1.0/h_rd.Integral())
h_mc.Scale(1.0/h_mc.Integral())

h_rd.SetXTitle('#Inv.\,mass(J/\psi\,K^+\,\pi^-\,\pi^+) \,\, with \,\, misid, GeV/c^2')
h_mc.SetXTitle('#Inv.\,mass(J/\psi\,K^+\,\pi^-\,\pi^+) \,\, with \,\, misid, GeV/c^2')

h_mc.Draw("")
h_rd.Draw("same")


db.close()