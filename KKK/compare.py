
ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetOptTitle(1)
ROOT.gStyle.SetOptStat(0)

canvas = ROOT.TCanvas("canvas", "MC vs RD", 700,700)
canvas.Divide(3,2)


[("pt_kaon", 25, -0.5, 6.0),
("eta_kaon", 25, -0.5, 6.0),

("y_jpsi", 25, 0.0, 5.0),
("eta_jpsi", 25, -0.5, 10.0),
("pt_jpsi", 25, -0.5, 20.0),

("y_b", 25, 0.0, 5.0),
("eta_b", 25, -0.5, 10.0),
("pt_b", 25, -0.5, 40.0)]

# ('m_jpsi', 25, 3.0, 3.2),

variables = [
    # (name, binning, range_min, range_max)
    # ('c2ip_b', 25, 0, 20, "\chi^{2}_{IP}(B)"),
    ("eta_b", 25, 1.0, 5.5, "\eta_{B}"),
    ("y_b", 25, 2.0, 5.0, "y_{B}"),

    ("pt_b", 25, -0.5, 40.0, "p_{T}(B)"),
    # ('minann_K', 25, 0, 1, "min(Prob.NN(K))"),
    # ('minann_mu', 25, 0.2, 1, "min(Prob.NN(\mu))"),

    ('minPt_track', 25, 0, 2, "min(p_{T}(track))"),
    ('y_jpsi', 25, 1.5, 5.0, "y_{J/\psi}"),
    ('pt_jpsi', 25, 0, 25, "p_{T}(J/\psi)"),

    # ('ann_kaon1', 25, 0.0, 1.0),
    # ('ann_kaon2', 25, 0.0, 1.0),
    # ('ann_kaon3', 25, 0.0, 1.0),
    # # ("p_b", 25, 0.0, 10.0),
]
hists_rd = []
hists_mc = []


mctrue = " && mcTrueB && mcTruePsi && mcTrueK1 && mcTrueK2 && mcTrueK3 && mcTrueMu1 && mcTrueMu2"

for v in variables:
    hh_rd = ROOT.TH1F("rd_h_" + v[0], '', v[1], v[2], v[3])
    hh_rd.Sumw2()
    # ds_Bu.project(hh_rd, v[0], "SBu_sw")
    tBu_mc_p6.project(hh_rd, v[0], cuts_Bu + mctrue)

    mc_name = v[0]
    # if 'ann_kaon' in v[0]:
    #     mc_name = v[0][:-1] + '[' + str(int(v[0][-1]) - 1) + ']'

    hh_mc = ROOT.TH1F("mc_h_" + v[0], '', v[1], v[2], v[3])
    hh_mc.Sumw2()
    tBu_mc_p8.project(hh_mc, mc_name, cuts_Bu + mctrue)

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


for i, name in enumerate(variables):
    canvas.cd(i+1)
    hists_rd[i].red()
    hists_mc[i].blue()

    hists_mc[i].Draw()
    hists_rd[i].Draw("same")


canvas.Update()
