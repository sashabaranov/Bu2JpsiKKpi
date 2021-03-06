import ROOT

# Bu meson
m_Bu = ROOT.RooRealVar('m_Bu', '', 5.1, 5.5)
nbin_Bu = 100

events_binning = int(m_Bu.getBinWidth(nbin_Bu) * 1000)

comparison_variables = [
    # ("DTFchi2ndof", "", -0.5, 20.0),
    # ("DTFctau", "", -0.5, 20.0),
    # ("c2ip_b", "", -0.5, 30.0),
    # ("m_jpsi", "", 2.5, 3.5),

    # ("minann_K", "", -0.5, 1.5),
    ("ann_kaon0", "", -0.5, 1.5, lambda s: s.ann_kaon[0]),
    ("ann_kaon1", "", -0.5, 1.5, lambda s: s.ann_kaon[1]),
    ("ann_kaon2", "", -0.5, 1.5, lambda s: s.ann_kaon[2]),

    ("ann_kaon_PI_0", "", -0.5, 1.5, lambda s: s.ann_kaon_PI[0]),
    ("ann_kaon_PI_1", "", -0.5, 1.5, lambda s: s.ann_kaon_PI[1]),
    ("ann_kaon_PI_2", "", -0.5, 1.5, lambda s: s.ann_kaon_PI[2]),

    # ("minann_mu", "", -0.5, 1.5),

    # # Interesting variables
    # ("pt_kaon", "", -0.5, 20.0),
    # ("eta_kaon", "", -0.5, 10.0),

    # ("y_jpsi", "", 0.0, 7.0),
    # ("eta_jpsi", "", -0.5, 10.0),
    # ("pt_jpsi", "", -0.5, 30.0),

    # ("y_b", "", 0.0, 10.0),
    # ("eta_b", "", -0.5, 20.0),
    # ("pt_b", "", -0.5, 40.0),

    # ("DTFm_KKK", "", -0.5, 40.0),

    # ("minPt_track", "", -0.5, 20.0),
    ("mass_k1aspi", "", -0.5, 10.0),
    ("mass_k3aspi", "", -0.5, 10.0),
]



selector_variables = [(m_Bu, lambda s: s.DTFm_b )] + comparison_variables
