import ROOT

# Bu meson
m_Bu = ROOT.RooRealVar('m_Bu', 'mass(J/\psi K \pi \pi)', 5.2, 5.4)
nbin_Bu = 75

# Resonanse
m_R = ROOT.RooRealVar('m_R', 'mass(J/\psi \pi \pi)', 0.0, 10.0)
nbin_R = 75

events_binning = int(m_Bu.getBinWidth(nbin_Bu) * 1000)

comparison_variables = [
    # ("DTFchi2ndof", "", -0.5, 20.0),
    # ("DTFctau", "", -0.5, 20.0),
    # ("c2ip_b", "", -0.5, 30.0),
    # ("m_jpsi", "", 2.5, 3.5),

    # ("minann_K", "", -0.5, 1.5),
    ("ann_kaon", "", -0.5, 1.5, lambda s: s.ann_kaon[0]),
    ("ann_pion1", "", -0.5, 1.5, lambda s: s.ann_pion[0]),
    ("ann_pion2", "", -0.5, 1.5, lambda s: s.ann_pion[1]),

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
    ("mass_pi1ask", "", -0.5, 10.0),
    ('m1234c2', "", 0.0, 10.0),
    ('m123c2', "", 0.0, 10.0),
    ('m124c2', "", 0.0, 10.0),
    ('m12c2', "", 0.0, 10.0),
    ('m134c2', "", 0.0, 10.0),
    ('m13c2', "", 0.0, 10.0),
    ('m14c2', "", 0.0, 10.0),
    ('m1c2', "", 0.0, 10.0),
    ('m234c2', "", 0.0, 10.0),
    ('m23c2', "", 0.0, 10.0),
    ('m24c2', "", 0.0, 10.0),
    ('m34c2', "", 0.0, 10.0),
]



selector_variables = [(m_Bu, lambda s: s.DTFm_b ), (m_R, lambda s: s.m134c2 )] + comparison_variables
