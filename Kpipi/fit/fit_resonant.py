from fit7 import ds_Bu
from resonant_model import *

masses = [
    ('m123c2', 3.9, 5.5),
    ('m124c2', 3.9, 5.9),
    ('m12c2', 3.0, 5.0),
    ('m134c2', 3.5, 5.0),
    ('m13c2', 3.0, 5.0),
    ('m14c2', 3.0, 5.0),
    ('m234c2', 0.5, 2.5),
    ('m23c2', 0.5, 2.0),
    ('m24c2', 0.5, 2.0),
    ('m34c2', 0.0, 2.0),
]


for m, a, b in masses:
    h = ROOT.TH1F("h_"+m, m, 75, a, b)
    h.Sumw2()
    h.SetXTitle(m)

    ds_Bu.project(h, m, "SBu_sw")

    h.Draw()
    canvas.SaveAs("pics/{}.png".format(m))


# h = ROOT.TH1F("h_R", "", 75, 3.6, 3.9)
# h.Sumw2()

# ds_Bu.project(h, "m_R", "SBu_sw")


# r, f = model_Resonant.fitHisto(h, draw=True)

# # # r("sigma_X")[1].release()
# # # r.signal.sigma.release()

# r, f = model_Resonant.fitHisto(h, draw=True)
# f.SetXTitle("Inv.\, mass(J/\psi \pi \pi), GeV/c^{2}")

# f.Draw()