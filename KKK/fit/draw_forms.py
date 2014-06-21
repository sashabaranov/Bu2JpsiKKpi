from tools import *

LHCbStyle()
db = shelve.open('$KKpidir/fit/histos.shelve')

histos = sorted([h for h in db['KKK']['MC'].values() if 'p6' in h.GetName()], key=lambda h: -h.GetMaximum())
# histos = sorted([h for h in db['KKK']['RD'].values()], key=lambda h: -h.GetMaximum())
# histos = sorted(db['KKK']['MC'].values(), key=lambda h: -h.GetMaximum())

titleX = '#Inv.\,mass(J/\psi\,KKK) \,\, with \,\, misid, GeV/c^2'
titleY = "Events / (%.1f \, MeV/c^{2})" % (histos[0].GetBinWidth(0) * 1000)

legend = ROOT.TLegend(0.55, 0.65, 0.86, 0.9);
legend.SetFillColor(ROOT.kWhite)



for i, h in enumerate(histos):
    col = i + 1

    h.SetLineColor(col)
    h.SetFillColor(col)
    h.SetMarkerColor(col)

    h.SetXTitle(titleX)
    h.SetYTitle(titleY)
    h.Draw('same')

    legend.AddEntry(h.GetName(), make_legend(h.GetName(), mc=True), "P")

legend.Draw()
