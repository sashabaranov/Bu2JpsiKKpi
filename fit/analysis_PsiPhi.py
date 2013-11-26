from ROOT import *
import ROOT
from AnalysisPython import Models

RAD = ROOT.RooAbsData

#------------------- no statistic ---------------
gStyle.SetOptStat(kFALSE)
gStyle.SetOptFit(kTRUE)

#------------------- define storage type --------
if RAD.Tree != RAD.getDefaultStorageType():
    print 'DEFINE default storage type to be TTree!'
    RAD.setDefaultStorageType (RAD.Tree)

#------------------- production dataset ---------
from data_Bc import tBu2 as chain

print "Variables"
#Psi Range
min_Kpi = 0.892 - 0.05
max_Kpi = 0.892 + 0.05

#Phi Range
min_KK = 0.990
max_KK = 1.060

mb = RooRealVar ( "mb", "Gev/c^{2}", 5.16, 5.45) 
m_Kstar = RooRealVar("m_Kstar","Mev/c^{2}",min_Kpi, max_Kpi)
m_KK = RooRealVar("m_KK", "Mev/c^{2}", min_KK, max_KK)

print "Make VarSet"
varset  = RooArgSet  (mb, m_Kstar, m_KK)
dataset = RooDataSet("dataset","Dataset",varset)

print "Make DataSet"
nEvents = chain.GetEntries()
print "Num events = ", nEvents

for evt in range(nEvents):
	if evt%1000 == 0: print "Event",evt,'of',nEvents

	if chain.GetEntry(evt) <= 0:
        	print "EVENT ", evt
        	continue

	mBu = chain.DTFm_b
	cuts = [chain.DTFchi2ndof > 0,
			chain.DTFchi2ndof < 5,
			chain.DTFctau > 0.2,
			chain.pt_kaon[0] > 0.6 ,
			chain.pt_kaon[1] > 0.6,
			chain.pt_pion > 0.3,
			chain.m_jpsi    > 3.020 ,
			chain.m_jpsi    < 3.135,
			chain.minann_K  > 0.2,
			(chain.psi3_l0tos & 2) == 2]

	if not all(cuts):
		continue

	mb.setVal(mBu) 
	m_Kstar.setVal(chain.m_Kstar)
	m_KK.setVal(chain.m_KK)
	dataset.add(varset)

#------------------- Cuts ------ ------------------
cut_Kpi = 'DTFm_kpi> {0} && DTFm_kpi< {1}'.format(min_Kpi, max_Kpi)
cut_KK = 'DTFm_kk> {0} && m_kk< {1}'.format(min_KK, max_KK)
cuts = cut_Kpi+' && '+cut_KK

#------------------- Parameters
mean_Bu = RooRealVar("mean_Bu", "mean value Bs", 5.279, 5.25, 5.45)
sigma_B = RooRealVar("sigma_B", "sigma B", 0.005, 0.1 * 0.005, 5 * 0.005)
coef = RooRealVar( "coef","coef of exp", -1, -10., 20.)

mean_Kpi = RooRealVar("mean_Kpi","mean of Gauss",3686.1, min_Kpi, max_Kpi)
sigma_Kpi =  RooRealVar("sigma_Kpi","sigma of gaussian", 2.93, 1.0, 10.) 

mean_KK = RooRealVar("mean_KK","mean of BW",1019.4, 1000., 1050.)
width_KK =  RooRealVar("width_KK","width of gaussian", 4.26, 1., 10.) 
L_KK = 1

m1 = 493.7
m2 = 493.7
m3 = 3686.1
mB = 5366.77
L = 0

p1 = RooRealVar( "p1","coef of pol", 1., -1000., 1000.)

#------------------- Yields
sBs   = RooRealVar ( "SBs"  , "Signal Bs", 3500, 0, 15000)
bgd_B     = RooRealVar ( "Bgd"    , "Background" , 3000, 0, 14000)
sB0   = RooRealVar ( "SB0"  , "Signal B0", 500, 0, 14000)

ss = RooRealVar ( "Psi_Phi"  , "Signal Psi Phi", 600, 0, 15000)
sb = RooRealVar ("Psi_ps", "Background" , 500, 0, 15000)
bs = RooRealVar ("pol_Phi", "Background" , 500, 0, 15000)
bb = RooRealVar ("pol_ps", "Background" , 500, 0, 15000)

#------------------- Different functions ------
gauss_Bs = RooGaussian ("gauss_Bs", "Signal gauss", mb, mean_Bu, sigma_B)
exp_bg = RooExponential ( "exp_bg","Exp background", mb, coef)
gauss_B0 = RooGaussian ("gauss_B0", "Signal gauss", mb, mean_B0, sigma_B)

gauss_Psi = RooGaussian ("gauss_Psi", "Signal gauss", m_Kstar, mean_Kpi, sigma_Kpi)
pol_Psi = Analysis.Models.PolyPositive("pol_Psi","pol", m_Kstar, p1, min_Kpi, max_Kpi)

BW_KK = Analysis.Models.BreitWigner('BW_KK','BreitWigner pdf', m_KK, mean_KK, width_KK, m1, m2, L_KK)
ps_KK = Analysis.Models.PhaseSpace23L('ps_KK','pss pdf', m_KK, m1, m2, m3, mB, L)
#------------------- Models -------------------
print "Create model"
pdfB = RooAddPdf("pdfB", "signal+bg+signal",RooArgList(gauss_Bs, exp_bg, gauss_B0),RooArgList (sBs, bgd_B, sB0))

gauss_BW_prod = RooProdPdf("gauss_BW_prod","BW's prod",gauss_Psi,BW_KK)
gauss_ps_prod = RooProdPdf("gauss_ps_prod","prod",gauss_Psi,ps_KK)
pol_BW_prod = RooProdPdf("pol_BW_prod","prod",pol_Psi,BW_KK)
pol_ps_prod = RooProdPdf("pol_ps_prod","prod",pol_Psi,ps_KK)


pdfSignal = RooAddPdf("pdfSignal", "signal+bg",RooArgList(gauss_BW_prod, pol_ps_prod, gauss_ps_prod, pol_BW_prod),RooArgList (ss, bb, sb, bs))
#------------------- Fit ----------------------
print "Fit B signal"

fit_res = pdfB.fitTo(dataset, RooFit.Save())

#------------------- plot B0 signal and fit result

mframe = mb.frame(40)
mframe.SetTitle('J/#psiK#piK#pi')
dataset.plotOn(mframe)
pdfB.plotOn(mframe)
pdfB.paramOn(mframe)


c1 = TCanvas("c1","signal")
c1.cd()
mframe.Draw()

#------------------- splot --------------------
print "Jpsipipi sPlot"
splot = RooStats.SPlot("splot", "An sPlot", dataset, pdfB,RooArgList(sBs, bgd_B, sB0))

ds = RooDataSet("ds","Dataset", dataset, dataset.get(), cuts,'SBs_sw')

histD = ds.createHistogram(m_Kstar,m_KK, 50, 50)
histD.SetTitle("#Psi(2S) / #Phi")
c2 = TCanvas("c2","c2")
c2.cd()
histD.Draw('lego2')

#------------------- fit splot result ----------
print "fit sPlot result"

fit_res = pdfSignal.fitTo(ds, RooFit.SumW2Error(kTRUE), RooFit.Minos(kFALSE), RooFit.Save())

xframe = m_Kstar.frame(50)
xframe.SetTitle('J/#psi#pi#pi fit projection')
ds.plotOn(xframe)
pdfSignal.plotOn(xframe)
pdfSignal.paramOn(xframe)

c3 = TCanvas("c3","c3")
c3.cd()
xframe.Draw()

yframe = m_KK.frame(50)
yframe.SetTitle('KK fit projection')
ds.plotOn(yframe)
pdfSignal.plotOn(yframe)
pdfSignal.paramOn(yframe)

c4 = TCanvas("c4","c4")
c4.cd()
yframe.Draw()

c2.cd()
hh = pdfSignal.createHistogram('m_Kstar:m_KK', 50, 50)
hh.SetLineColor(63)
hh.Draw("surf same")

#------------------- print plots -----------------
print "print plots"

#c1.Print('./pictures/signal.gif')
#c2.Print('./pictures/histD.gif')
#c3.Print('./pictures/histPsi.gif')
#c4.Print('./pictures/histPhi.gif')
