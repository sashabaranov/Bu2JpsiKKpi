from tools import *

db = shelve.open('$KKpidir/fit/histos.shelve')

# a = [(x, True) for x in [h for h in db['KKK']['MC'].values() if 'p6' in h.GetName()][:2]]
# b = [(x, False) for x in [h for h in db['KKK']['RD'].values()][:2]]


# histos = sorted(a+b, key=lambda h: -h[0].GetMaximum());
histos = [(x,False) for x in db['KKK']['RD'].values()]
# histos = sorted(db['KKK']['MC'].values(), key=lambda h: -h.GetMaximum())
# histos = a + b

# histos = [(x, True) for x in [h for h in db['KKK']['MC'].values() if 'p6' in h.GetName()]]


for el in histos:
    h, MC = el
    h.Sumw2()
    h.Scale(1.0/h.Integral())
    h.Smooth(20)
    # for i in xrange(h.GetNbinsX()):
    #     h.SetBinError(i, 0.0)



histos = sorted(histos, key=lambda el: -el[0].GetMaximum())

titleX = '#Inv.\,mass(J/\psi\,KKK) \,\, with \,\, misid, GeV/c^2'
# titleY = "Events / (%.1f \, MeV/c^{2})" % (histos[0][0].GetBinWidth(0) * 1000)
titleY = "Shape(normailzed)"

legend = ROOT.TLegend(0.2, 0.65, 0.4, 0.9);
legend.SetFillColor(ROOT.kWhite)



for i, el in enumerate(histos):
    h, MC = el
    col = i + 1

    h.SetLineWidth(3)
    h.SetLineColor(col)
    # h.SetFillColor(col)
    h.SetMarkerColor(col)

    h.SetXTitle(titleX)
    h.SetYTitle(titleY)

    y_axis = h.GetYaxis()
    y_axis.SetTitleOffset(1.1)

    h.Draw('hist same c')

    legend.AddEntry(h.GetName(), make_legend(h.GetName(), mc=MC), "L")

legend.Draw()
