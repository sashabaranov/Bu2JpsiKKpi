import ROOT

# Bu meson
m_Bu = ROOT.RooRealVar('DTFm_b', 'mass(J/psiKKpi)', 5.2, 5.4)
nbin_Bu = 80

events_binning = (m_Bu.getMax() - m_Bu.getMin()) / nbin_Bu * 1000

m_Bu_a = ROOT.RooRealVar('m_Bu', 'mass(J/psiKKpi)', 5.16, 5.45)


# K*(892)
m_Kstar = ROOT.RooRealVar('DTFm_kpi', 'mass(K+ pi-)', 0.634, 1.2)
nbin_Kstar = 20

mm_Kstar = ROOT.RooRealVar('mm_Kstar', 'mass(Kstar0)', 0.85, 0.9)
mm_Kstar.fix(0.892)

width_Kstar = ROOT.RooRealVar('width_Kstar', 'width(Kstar0)', 0.0, 1.0)

low_kpi = ROOT.RooRealVar('low', 'low_KstarBg', 0.6, 0.7)
low_kpi.fix(0.634)
high_kpi = ROOT.RooRealVar('high', 'high_KstarBg', 1.3, 2.0)



# Phi(1020)
m_Phi = ROOT.RooRealVar('DTFm_kk', 'mass(KK)', 1.0, 1.08)
nbin_Phi = 30

mm_Phi = ROOT.RooRealVar('mm_Phi', 'mass(Phi)', 1.0, 1.1)
mm_Phi.fix(1.02)

width_Phi = ROOT.RooRealVar('width_Phi', 'width(Phi)', 0.0, 1.0)

low_kk = ROOT.RooRealVar('low', 'low_PhiBg', 0.9, 1.0)
low_kk.fix(0.986)
high_kk = ROOT.RooRealVar('high', 'high_PhiBg', 2.0, 2.1)
#high_kk.fix(2.04)




comparison_variables = [
    # (name, binning, range_min, range_max, name)
    # ("DTFchi2ndof", 25, 0.0, 6.0, "\chi^{2}_{DTF}/ndof"),
    # ("DTFctau", 25, 0.2, 1.0, "c tau"),
    # ("pt_b", 25, -0.5, 35.0, "p_{T}(B)"),


    # ('minann_K', 25, 0, 1, "min(Prob.NN(K))"),
    # ('minann_mu', 25, 0.2, 1, "min(Prob.NN(\mu))"),
    # ('minann_pi', 25, 0, 1, "min(Prob.NN(\mu))"),

    # ('minPt_track', 25, 0, 2, "min(p_{T}(track))"),
    # ("eta_b", 25, 1.0, 5.5, "\eta_{B}"),
    # ('pt_jpsi', 25, 0, 15, "p_{T}(J/\psi)"),

    # ("y_b", 25, 2.0, 5.0, "y_{B}"),


    # ("pt_kaon[0]", 25, -0.5, 5.0, "p_{T}(K1)"),
    # ("pt_kaon[1]", 25, -0.5, 5.0, "p_{T}(K2)"),


    # ('y_jpsi', 25, 1.5, 5.0, "y_{J/\psi}"),

    # ('ann_kaon1', 25, 0.0, 1.0),
    # ('ann_kaon2', 25, 0.0, 1.0),
    # ('ann_kaon3', 25, 0.0, 1.0),

    # Kinematics
    # ("pt_pion", "p_{T}(\pi)", -0.5, 10.0, lambda s: s.pt_pion[0]),
    # ("eta_pion[0]", "\eta(\pi)", 2.0, 5.5),
    # ("p_pion[0]", "p(\pi)", 3.0, 20.0),

    # ("pt_kaon[0]", "p_{T}(K1)", -0.5, 10.0),
    # ("eta_kaon[0]", "\eta(K1)", 2.0, 5.5),
    # ("p_kaon[0]", "p(K1)",  3.0, 20.0),

    # ("pt_kaon[1]", "p_{T}(K2)", -0.5, 10.0),
    # ("eta_kaon[1]", "\eta(K2)", 2.0, 5.5),
    # ("p_kaon[1]", "p(K2)",  3.0, 20.0),
    ( 'mass_PI_as_K' , "Inv.mass(misid\,\pi), GeV/c^2", 0.0, 10.0),
    ( 'mass_K_as_PI' , "Inv.mass(misid\,\pi), GeV/c^2", 0.0, 10.0),
    ( 'DTFm_kk', "Inv. \, Mass(KK), GeV/c^2", 0.0, 20.0),
    ( 'DTFm_kpi', "Inv. \, Mass(K\pi), GeV/c^2", 0.0, 20.0),
]

selector_variables = comparison_variables + [(m_Bu, lambda s: s.DTFm_b )]
