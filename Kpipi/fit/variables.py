import ROOT

# Bu meson
m_Bu = ROOT.RooRealVar('m_Bu', '', 5.15, 5.45)
nbin_Bu = 100
binning_b = (m_Bu.getMax() - m_Bu.getMin()) / nbin_Bu * 1000


# Resonanse
m_R = ROOT.RooRealVar('m_R', 'mass(J/\psi \pi \pi)', 3.65, 3.75)
nbin_R = 75

events_binning = int(m_Bu.getBinWidth(nbin_Bu) * 1000)

comparison_variables = [
    # ("DTFchi2ndof", "", -0.5, 20.0),
    # ("DTFctau", "", -0.5, 20.0),
    # ("c2ip_b", "", -0.5, 30.0),
    # ("m_jpsi", "", 2.5, 3.5),

    # ("minann_K", "", -0.5, 1.5),
    # ("ann_kaon", "", -0.5, 1.5, lambda s: s.ann_kaon[0]),
    # ("ann_pion1", "", -0.5, 1.5, lambda s: s.ann_pion[0]),
    # ("ann_pion2", "", -0.5, 1.5, lambda s: s.ann_pion[1]),

    ("ann_pion_K", "", -0.5, 1.5, lambda s: s.ann_pion_K[0]),
    # # ("minann_mu", "", -0.5, 1.5),

    # # # Interesting variables
    # # ("pt_kaon", "", -0.5, 20.0),
    # # ("eta_kaon", "", -0.5, 10.0),

    # # ("y_jpsi", "", 0.0, 7.0),
    # # ("eta_jpsi", "", -0.5, 10.0),
    # # ("pt_jpsi", "", -0.5, 30.0),

    # # ("y_b", "", 0.0, 10.0),
    # # ("eta_b", "", -0.5, 20.0),
    # # ("pt_b", "", -0.5, 40.0),

    # # ("DTFm_KKK", "", -0.5, 40.0),

    # # ("minPt_track", "", -0.5, 20.0),
    ("mass_pi2ask", "", -0.5, 40.0),
    ('m134c2', "", -0.5, 40.0),
    ('m234c2', "", -0.5, 40.0),
    ('m24c2', "", -0.5, 40.0),
    ('m34c2', "", -0.5, 40.0),
]



selector_variables = [(m_Bu, lambda s: s.DTFm_b )] + comparison_variables
