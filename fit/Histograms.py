import ROOT
from AnalysisPython.PyRoUts import hID

from Variables import *

# Histograms
h1 = ROOT.TH1F(hID(), '', nbin_Bu, m_Bu.getMin(), m_Bu.getMax())
h1.Sumw2()

h1a = ROOT.TH1F(hID(), '', nbin_Bu, m_Bu.getMin(), m_Bu.getMax())
h1a.Sumw2()

h2 = ROOT.TH1F(hID(), '', nbin_Kstar, m_Kstar.getMin(), m_Kstar.getMax())
h2.Sumw2()

h3 = ROOT.TH1F(hID(), '', nbin_Phi, m_Phi.getMin(), m_Phi.getMax())
h3.Sumw2()

h4 = ROOT.TH1F(hID(), '', nbin_Bu, 1.0, 3.5)
h4.Sumw2()
