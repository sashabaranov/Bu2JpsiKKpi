#!/usr/bin/python
# =============================================================================
# @file DiCharm/Analysis/Models.py
#
#  Set of PDFs for 2xCharm fits
#
#  @author Vanya BELYAEV Ivan.Belyaeve@Cern.ch
#  @date 2011-07-25
#
#                    $Revision$
#  Last modification $Date$
#                 by $Author$
# =============================================================================
"""
Set of PDFs for 2xCharm fits
"""
# =============================================================================
__version__ = "$Revision:"
__author__ = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__ = "2011-07-25"
# =============================================================================
__all__ = (
    #
    'Bukin_pdf',  # generic Bukin PDF
    'CB2_pdf',  # double-sided cristal ball
    'Gauss_pdf',  # simple  Gauss
    'Bkg_pdf',  # Background: exponential
    #
    'D0_pdf',  # PDF for D0        : Bukin
    'Dp_pdf',  # PDF for D+        : Bukin
    'Ds_pdf',  # PDF for Ds+       : Bukin
    'Lc_pdf',  # PDF for Lambda_c+ : Gauss
    'UpsFit_pdf',  # PDF for Upsilon -> mu+mu- spectrum
    #
    'Charm1_pdf',  # PDF
    'Charm2_pdf',  # PDF for 2xCharm fit
    'Charm3_pdf',
)
# =============================================================================
import ROOT
import math
from AnalysisPython.PyRoUts import VE, hID, allInts, cpp, rootID
from AnalysisPython.Logger import getLogger
# =============================================================================
logger = getLogger(__name__)
# =============================================================================
_nemax = 80000  # number of events per CPU-core
_ncmax = 4  # maximal number of CPUs: there are some problems with >= 7
                # @see https://sft.its.cern.ch/jira/browse/ROOT-4897
#


def ncpu(events):
    #
    n = events // _nemax
    if n <= 1:
        return ROOT.RooFit.NumCPU(1)
    #
    import multiprocessing
    n_cores = multiprocessing.cpu_count()
    #
    return ROOT.RooFit.NumCPU(min(n, n_cores, _ncmax))


class ModelProxy(object):
    """In case we have only pdf, not model itself."""
    def __init__(self, pdf):
        self.pdf = pdf


# =============================================================================
# @class Bukin_pdf
#  simple wrapper over Bukin-pdf
#  @see RooBukinPdf
#  @author Vanya BELYAEV Ivan.Belyaeve@Cern.ch
#  @date 2011-07-25
class Bukin_pdf(object):

    """
    Define PDFs for D0,D+,Ds+
    """

    def __init__(self,
                 name,
                 mn,
                 mx,
                 fixMass=None,
                 fixSigma=None,
                 fixXi=None,
                 fixRhoL=None,
                 fixRhoR=None,
                 mass=None,
                 mean=None,
                 sigma=None,
                 xi=None,
                 rhol=None,
                 rhor=None):

        mean_0 = 1.8627
        sigma_0 = 0.00729757

        xi_0 = -2.8893e-02
        rhol_0 = -1.00529e-02
        rhor_0 = -2.55883e-01

        if 0 <= name.find("Dp") or 0 <= name.find("D+"):
            mean_0 = 1.86760e+00
            sigma_0 = 6.93069e-03
        elif 0 <= name.find("Ds") or 0 <= name.find("D_s+"):
            mean_0 = 1.96714e+00
            sigma_0 = 5.75354e-03
        elif 0 <= name.find("Lc") or 0 <= name.find("Lc+"):
            mean_0 = 2.28596e+00
            sigma_0 = 5.11791e-03
        elif 0 <= name.find("Bc"):
            mean_0 = 6.2765832270928321
            sigma_0 = 1.33172e-02

        if fixMass and isinstance(fixMass, float):
            mean_0 = fixMass
        if fixSigma and isinstance(fixSigma, float):
            sigma_0 = fixSigma
        if fixXi and isinstance(fixXi, float):
            xi_0 = fixXi
        if fixRhoL and isinstance(fixRhoL, float):
            rhol_0 = fixRhoL
        if fixRhoR and isinstance(fixRhoR, float):
            rhor_0 = fixRhoR

        if not mass:
            self.mass = ROOT.RooRealVar(
                "m_%s" % name, "mass(%s)" % name,            mn, mx)
        else:
            self.mass = mass

        if not mean:
            self.mean = ROOT.RooRealVar(
                "mean_%s" % name, "mean(%s)" % name, mean_0, mean_0 - sigma_0, mean_0 + sigma_0)
            if fixMass:
                self. mean  . setConstant(True)
        else:
            self.mean = mean

        if not sigma:
            self.sigma = ROOT.RooRealVar("sigma_%s" % name, "sigma(%s)" %
                                         name, sigma_0,    0.5 * sigma_0,    2.5 * sigma_0)
            if fixSigma:
                self. sigma . setConstant(True)
        else:
            self.sigma = sigma

        if not xi:
            self.xi = ROOT.RooRealVar(
                "xi_%s" % name, "xi(%s)" % name, xi_0, -1, 1)
            if fixXi:
                self. xi    . setConstant(True)
        else:
            self.xi = xi

        if not rhol:
            #
            if isinstance(fixRhoL, (int, long, float)):
                self.rhol = ROOT.RooRealVar(
                    "rhol_%s" % name, "rhol(%s)" % name, fixRhoL, -2.0, 0.0)
            else:
                self.rhol = ROOT.RooRealVar(
                    "rhol_%s" % name, "rhol(%s)" % name, 0.0, -2.0, 0.0)
            #
            if fixRhoL or isinstance(fixRhoL, (int, long, float)):
                self.rhol.setConstant(True)
        else:
            self.rhol = rhol

        if not rhor:

            if isinstance(fixRhoR, (int, long, float)):
                self.rhor = ROOT.RooRealVar(
                    "rhor_%s" % name, "rhor(%s)" % name, fixRhoR, -2.0, 0.0)
            else:
                self.rhor = ROOT.RooRealVar(
                    "rhor_%s" % name, "rhor(%s)" % name, 0.0, -2.0, 0.0)
             #
            if fixRhoR or isinstance(fixRhoR, (int, long, float)):
                self.rhor.setConstant(True)
            #
        else:
            self.rhor = rhor

        #self.rhol . setConstant ( True )
        #self.rhol . setConstant ( True )

        self.pdf = ROOT.RooBukinPdf("bkn_" + name,
                                    "Bukin(%s)" % name,
                                    self.mass,
                                    self.mean,
                                    self.sigma,
                                    self.xi,
                                    self.rhol,
                                    self.rhor)

        # fit
        def fitTo(self, dataset, *args, **kwargs):
            """
            Perform the fit
            """
            return self.pdf.fitTo(dataset, *args, **kwargs)


# =============================================================================
# @class CB2_pdf
#  simple wrapper over double-sided Cristal Ball function
#  @see Analysis::Models::CrystalBallDS
#  @author Vanya BELYAEV Ivan.Belyaev@Cern.ch
#  @date 2011-07-25
class CB2_pdf(object):

    """
    Define double sided Crystall Ball
    """

    def __init__(self,
                 name,
                 mn,
                 mx,
                 fixMass=None,
                 fixSigma=None,
                 fixAlphaL=None,
                 fixAlphaR=None,
                 fixNL=None,
                 fixNR=None,
                 mass=None,
                 mean=None,
                 sigma=None,
                 alphaL=None,
                 alphaR=None,
                 nL=None,
                 nR=None):

        mean_0 = 1.8627
        sigma_0 = 0.00729757

        aL_0 = 1
        nL_0 = 0

        aR_0 = 1
        nR_0 = 1

        if 0 <= name.find("Dp") or 0 <= name.find("D+"):
            mean_0 = 1.86760e+00
            sigma_0 = 6.93069e-03
        elif 0 <= name.find("Ds") or 0 <= name.find("D_s+"):
            mean_0 = 1.96714e+00
            sigma_0 = 5.75354e-03
        elif 0 <= name.find("Lc") or 0 <= name.find("Lc+"):
            mean_0 = 2.28596e+00
            sigma_0 = 5.11791e-03
        elif 0 <= name.find("Bc"):
            mean_0 = 6.2767
            sigma_0 = 0.013
        elif 0 <= name.find("Bs"):
            mean_0 = 5.366
            sigma_0 = 0.013
        elif 0 <= name.find("B"):
            mean_0 = 5.278
            sigma_0 = 0.013

        if fixMass and isinstance(fixMass, float):
            mean_0 = fixMass
        if fixSigma and isinstance(fixSigma, float):
            sigma_0 = fixSigma
        if fixAlphaL and isinstance(fixAlphaL, float):
            aL_0 = fixAlphaL
        if fixAlphaR and isinstance(fixAlphaR, float):
            aR_0 = fixAlphaR
        if fixNL and isinstance(fixNL, float):
            nL_0 = fixNL
        if fixNR and isinstance(fixNR, float):
            nR_0 = fixNR

        if not mass:
            self.mass = ROOT.RooRealVar(
                "m_%s" % name, "mass(%s)" % name,            mn, mx)
        else:
            self.mass = mass

        if not mean:
            self.mean = ROOT.RooRealVar(
                "mean_%s" % name, "mean(%s)" % name, mean_0, mean_0 - sigma_0, mean_0 + sigma_0)
            if fixMass:
                self. mean  . setConstant(True)
        else:
            self.mean = mean

        if not sigma:
            self.sigma = ROOT.RooRealVar("sigma_%s" % name, "sigma(%s)" %
                                         name, sigma_0,    0.5 * sigma_0,    2.5 * sigma_0)
            if fixSigma:
                self. sigma . setConstant(True)
        else:
            self.sigma = sigma

        if not alphaL:
            self.aL = ROOT.RooRealVar(
                "aL_%s" % name, "#alpha_{L}(%s)" % name, aL_0, 0, 10)
            if fixAlphaL:
                self. aL    . setConstant(True)
        else:
            self.aL = alphaL

        if not alphaR:
            self.aR = ROOT.RooRealVar(
                "aR_%s" % name, "#alpha_{R}(%s)" % name, aR_0, 0, 10)
            if fixAlphaR:
                self. aR    . setConstant(True)
        else:
            self.aR = alphaR

        if not nL:
            self.nL = ROOT.RooRealVar(
                "nL_%s" % name, "n_{L}(%s)" % name, nL_0, 0, 10)
            if fixNL:
                self. nL    . setConstant(True)
        else:
            self.nL = nL

        if not nR:
            self.nR = ROOT.RooRealVar(
                "nR_%s" % name, "n_{R}(%s)" % name, nR_0, 0, 50)
            if fixNR:
                self. nR    . setConstant(True)
        else:
            self.nR = nR

        self.pdf = cpp.Analysis.Models.CrystalBallDS(
            "cb2_" + name,
            "CB_{2}(%s)" % name,
            self.mass,
            self.mean,
            self.sigma,
            self.aL,
            self.nL,
            self.aR,
            self.nR)

        # fit
        def fitTo(self, dataset, *args, **kwargs):
            """
            Perform the fit
            """
            return self.pdf.fitTo(dataset, *args, **kwargs)


# =============================================================================
# @class Gauss_pdf
#  simple wrapper over Gaussian-pdf
#  @see RooGaussian
#  @author Vanya BELYAEV Ivan.Belyaeve@Cern.ch
#  @date 2011-07-25
class Gauss_pdf(object):

    """
    Simple Gaussian
    """

    def __init__(self,
                 name,
                 mn,
                 mx,
                 fixMass,
                 fixSigma,
                 mass,
                 mean=None,
                 sigma=None):

        mean_0 = 1.864
        sigma_0 = 0.0075
        if 0 <= name.find("Dp") or 0 <= name.find("D+"):
            mean_0 = 1.8676
            sigma_0 = 0.0068
        elif 0 <= name.find("Ds") or 0 <= name.find("Ds+"):
            mean_0 = 1.9672
            sigma_0 = 0.0068
        elif 0 <= name.find("Lc") or 0 <= name.find("Lc+"):
            mean_0 = 2.28510e+00
            sigma_0 = 5.93874e-03
        elif 0 <= name.find("Bc"):
            mean_0 = 5.278 + 0.9946
            sigma_0 = 0.009
        elif 0 <= name.find("B"):
            mean_0 = 5.278
            sigma_0 = 0.009

        if isinstance(fixMass, float):
            mean_0 = fixMass
        if isinstance(fixSigma, float):
            sigma_0 = fixSigma

        if fixMass:
            mean_0 = fixMass
        if fixSigma:
            sigma_0 = fixSigma

        if not mass:
            self.mass = ROOT.RooRealVar(
                "m_%s" % name, "mass(%s)" % name,            mn, mx)
        else:
            self.mass = mass

        if not mean:
            self.mean = ROOT.RooRealVar(
                "mean_%s" % name, "mean(%s)" % name, mean_0, mean_0 - sigma_0, mean_0 + sigma_0)
            if fixMass:
                self.mean .setConstant(True)
        else:
            self.mean = mean

        if not sigma:
            self.sigma = ROOT.RooRealVar("sigma_%s" % name, "sigma(%s)" %
                                         name, sigma_0,    0.25 * sigma_0,    4.0 * sigma_0)
            if fixSigma:
                self.sigma.setConstant(True)
        else:
            self.sigma = sigma

        self.pdf = ROOT.RooGaussian("gau_" + name,
                                    "Gauss(%s)" % name,
                                    self.mass,
                                    self.mean,
                                    self.sigma)

        # fit
        def fitTo(self, dataset, *args, **kwargs):
            """
            Perform the fit
            """
            return self.pdf.fitTo(dataset, *args, **kwargs)

# =============================================================================
# @class Gauss_pdf
#  simple wrapper over exponential
#  @see RooGaussian
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date 2011-07-25


class Bkg_pdf(object):

    """
    Define pdf for background : Exponential
    """
    # constructor

    def __init__(self,
                 name,
                 mass,
                 power=0):

        self.power = power

        import math
        change = 1.e+6

        taumin = math.log(change) / (mass.getMax() - mass.getMin())
        taumin = -1 * abs(taumin)
        taumax = abs(taumin)

        self.tau = ROOT.RooRealVar(
            "tau_%s" % name, "tau(%s)" % name, 0, taumin, taumax)
        self.phis = []

        if 1 == power:
            self.phi = ROOT.RooRealVar(
                'phi_%s' % name, 'phi(%s)' % name,  0, -3.5, 3.5)
            self.phis.append(self.phi)
            self.pdf = cpp.Analysis.Models.ExpoPositive('expopos_%s' % name,
                                                        'expopos(%s)' % name,
                                                        mass,
                                                        self.tau,
                                                        self.phi,
                                                        mass.getMin(),
                                                        mass.getMax())
        elif 2 == power:
            self.phi1 = ROOT.RooRealVar(
                'phi1_%s' % name, 'phi1(%s)' % name,  0, -3.5, 3.5)
            self.phi2 = ROOT.RooRealVar(
                'phi2_%s' % name, 'phi2(%s)' % name,  0, -3.5, 3.5)
            self.phis.append(self.phi1)
            self.phis.append(self.phi2)
            self.pdf = cpp.Analysis.Models.ExpoPositive('expopos_%s' % name,
                                                        'expopos(%s)' % name,
                                                        mass,
                                                        self.tau,
                                                        self.phi1,
                                                        self.phi2,
                                                        mass.getMin(),
                                                        mass.getMax())
        elif 3 == power:
            self.phi1 = ROOT.RooRealVar(
                'phi1_%s' % name, 'phi1(%s)' % name,  0, -3.5, 3.5)
            self.phi2 = ROOT.RooRealVar(
                'phi2_%s' % name, 'phi2(%s)' % name,  0, -3.5, 3.5)
            self.phi3 = ROOT.RooRealVar(
                'phi3_%s' % name, 'phi3(%s)' % name,  0, -3.5, 3.5)
            self.phis.append(self.phi1)
            self.phis.append(self.phi2)
            self.phis.append(self.phi3)
            self.pdf = cpp.Analysis.Models.ExpoPositive('expopos_%s' % name,
                                                        'expopos(%s)' % name,
                                                        mass,
                                                        self.tau,
                                                        self.phi1,
                                                        self.phi2,
                                                        self.phi3,
                                                        mass.getMin(),
                                                        mass.getMax())
        elif 4 == power:
            self.phi1 = ROOT.RooRealVar(
                'phi1_%s' % name, 'phi1(%s)' % name,  0, -3.5, 3.5)
            self.phi2 = ROOT.RooRealVar(
                'phi2_%s' % name, 'phi2(%s)' % name,  0, -3.5, 3.5)
            self.phi3 = ROOT.RooRealVar(
                'phi3_%s' % name, 'phi3(%s)' % name,  0, -3.5, 3.5)
            self.phi4 = ROOT.RooRealVar(
                'phi4_%s' % name, 'phi4(%s)' % name,  0, -3.5, 3.5)
            self.phis.append(self.phi1)
            self.phis.append(self.phi2)
            self.phis.append(self.phi3)
            self.phis.append(self.phi4)
            self.pdf = cpp.Analysis.Models.ExpoPositive('expopos_%s' % name,
                                                        'expopos(%s)' % name,
                                                        mass,
                                                        self.tau,
                                                        self.phi1,
                                                        self.phi2,
                                                        self.phi3,
                                                        self.phi4,
                                                        mass.getMin(),
                                                        mass.getMax())
        else:
            self.pdf = ROOT.RooExponential(
                'exp_%s' % name, 'exp(%s)' % name, mass, self.tau)

    # fit
    def fitTo(self, dataset, *args, **kwargs):
            """
            Perform the fit
            """
            return self.pdf.fitTo(dataset, *args, **kwargs)

# =============================================================================
# @class D0_pdf
#  simple wrapper over Bukin-pdf
#  @see RooBukinPdf
#  @author Vanya BELYAEV Ivan.Belyaeve@Cern.ch
#  @date 2011-07-25


class D0_pdf(Bukin_pdf):

    """
    Define PDF for D0
    """

    def __init__(self,
                 name,
                 fixMass=None,
                 fixSigma=None,
                 fixXi=-0.00044,
                 fixRhoL=-0.170619,
                 fixRhoR=-0.289993,
                 mass=None,
                 mean=None,
                 sigma=None,
                 xi=None,
                 rhol=None,
                 rhor=None):

        Bukin_pdf.__init__(self,
                           name,
                           1.82,
                           1.92,
                           fixMass,
                           fixSigma,
                           fixXi,
                           fixRhoL,
                           fixRhoR,
                           mass,
                           mean,
                           sigma,
                           xi,
                           rhol,
                           rhor)

# =============================================================================
# @class Dp_pdf
#  simple wrapper over Bukin-pdf
#  @see RooBukinPdf
#  @author Vanya BELYAEV Ivan.Belyaeve@Cern.ch
#  @date 2011-07-25


class Dp_pdf(Bukin_pdf):

    """
    Define PDF for D+
    """

    def __init__(self,
                 name,
                 fixMass=None,
                 fixSigma=None,
                 fixXi=-2.44474e-04,
                 fixRhoL=-1.07796e-01,
                 fixRhoR=-2.84132e-01,
                 mass=None,
                 mean=None,
                 sigma=None,
                 xi=None,
                 rhol=None,
                 rhor=None):

        Bukin_pdf.__init__(self,
                           name,
                           1.82,
                           1.91,
                           fixMass,
                           fixSigma,
                           fixXi,
                           fixRhoL,
                           fixRhoR,
                           mass,
                           mean,
                           sigma,
                           xi,
                           rhol,
                           rhor)

# =============================================================================
# @class Ds_pdf
#  simple wrapper over Bukin-pdf
#  @see RooBukinPdf
#  @author Vanya BELYAEV Ivan.Belyaeve@Cern.ch
#  @date 2011-07-25


class Ds_pdf(Bukin_pdf):

    """
    Define PDF for Ds+
    """

    def __init__(self,
                 name,
                 fixMass=1.9672,
                 fixSigma=0.0068,
                 fixXi=-6.45755e-04,
                 fixRhoL=-9.25349e-02,
                 fixRhoR=-1.86051e-01,
                 mass=None,
                 mean=None,
                 sigma=None,
                 xi=None,
                 rhol=None,
                 rhor=None):

        Bukin_pdf.__init__(self,
                           name,
                           1.90,
                           2.04,
                           fixMass,
                           fixSigma,
                           fixXi,
                           fixRhoL,
                           fixRhoR,
                           mass,
                           mean,
                           sigma,
                           xi,
                           rhol,
                           rhor)


# =============================================================================
# @class Lc_pdf
#  simple wrapper over Bukin-pdf
#  @see RooBukinPdf
#  @author Vanya BELYAEV Ivan.Belyaeve@Cern.ch
#  @date 2011-07-25
class Lc_pdf(Bukin_pdf):

    """
    Define PDF for Lc+
    """

    def __init__(self,
                 name,
                 fixMass=2.28590e+00,
                 fixSigma=5.11874e-03,
                 fixXi=1.82493e-03,
                 fixRhoL=-1.83140e-01,
                 fixRhoR=-2.40318e-01,
                 mass=None,
                 mean=None,
                 sigma=None,
                 xi=None,
                 rhol=None,
                 rhor=None):

        Bukin_pdf.__init__(self,
                           name,
                           2.24,
                           2.33,
                           fixMass,
                           fixSigma,
                           fixXi,
                           fixRhoL,
                           fixRhoR,
                           mass,
                           mean,
                           sigma,
                           xi,
                           rhol,
                           rhor)

# =============================================================================
# @class Jpsi_pdf
#  Wrapper over CrystalBall function, using
#   alpha/sigma parameterization by Matt Needham
#  @see RooBukinPdf
#  @author Vanya BELYAEV Ivan.Belyaeve@Cern.ch
#  @date 2011-07-25


class Jpsi_pdf(object):

    """
    Define PDF for Jpsi
    """

    def __init__(self,
                 name,
                 fixMass,
                 fixSigma,
                 fixN=1.0,
                 mass=None):

        mean_0 = 3.093
        sigma_0 = 0.013

        mn = 3.0
        mx = 3.2

        if fixMass:
            mean_0 = fixMass
        if fixSigma:
            sigma_0 = fixSigma

        n_0 = 1
        if fixN:
            n_0 = fixN

        n_0 = max(n_0, 1)

        if mass:
            self.mass = mass
        else:
            self.mass = ROOT.RooRealVar(
                "m_%s" % name, "mass(%s)" % name,            mn, mx)

        self.mean = ROOT.RooRealVar("means", "mean(%s)" %
                                    name, mean_0, mean_0 - sigma_0, mean_0 + sigma_0)
        self.sigma = ROOT.RooRealVar("sigma", "sigma(%s)" %
                                     name, sigma_0,    0.5 * sigma_0,    2.5 * sigma_0)

        # alpha/sigma parameterization by Matt Needham
        _unit = 1000
        if self.mass.getMin() <= 3096 <= self.mass.getMax():
            _unit = 1
        #
        self.a0 = ROOT.RooRealVar("a0", "a0(%s)" % name,  1.975, -10, 10)
        self.a1 = ROOT.RooRealVar(
            "a1", "a1(%s)" % name, -0.0011 * _unit, -10 * _unit, 10 * _unit)
        self.a2 = ROOT.RooRealVar(
            "a2", "a2(%s)" % name, -0.00018 * _unit ** 2, -10 * _unit ** 2, 10 * _unit ** 2)

        self.a0.setConstant(True)
        self.a1.setConstant(True)
        self.a2.setConstant(True)

        self.n = ROOT.RooRealVar(
            "n_%s" % name, "n(%s)" % name, n_0, 1, 5)

        if fixN:
            self.n.setConstant(True)

        self.aset = ROOT.RooArgList(self.a0,
                                    self.a1,
                                    self.a2,
                                    self.sigma)

        self.alpha = ROOT.RooFormulaVar("alpha_%s" % name,
                                        "alpha(%s)" % name,
                                        "a0 + a1*sqrt(sigma*sigma) + a2*sigma*sigma",
                                        self.aset)

        self.pdf = ROOT.RooCBShape("pdf_%s" % name, "CrystalBall(%s)" % name,
                                   self.mass,
                                   self.mean,
                                   self.sigma,
                                   self.alpha,
                                   self.n)

        if fixMass:
            self . mean  . setConstant(True)
        if fixSigma:
            self . sigma . setConstant(True)

    # fit
    def fitTo(self, dataset, *args, **kwargs):
            """
            Perform the fit
            """
            return self.pdf.fitTo(dataset, *args, **kwargs)


# =============================================================================
# @class Z0_pdf
#  Wrapper over CrystalBall function
#  @author Vanya BELYAEV Ivan.Belyaeve@Cern.ch
#  @date 2011-07-25
class Z0_pdf(object):

    """
    Define PDF for Z0
    """

    def __init__(self,
                 name,
                 mass=None,
                 fixMass=True,
                 fixSigma=True):

        if mass:
            self.m = mass
        else:
            self.m = ROOT.RooRealVar(
                "m_%s" % name, "mass(%s)" % name,        60, 120)

        self.mass = self.m

        self.mean = ROOT.RooRealVar(
            "mean_%s" % name, "mean(%s)" % name, 90.7, 85,  95)
        self.sigma = ROOT.RooRealVar(
            "sigma_%s" % name, "sigma(%s)" % name,  3.0, 2.0, 5.0)
        self.n = ROOT.RooRealVar(
            "n_%s" % name, "n(%s)" % name, 1, 1, 5)
        self.alpha = ROOT.RooRealVar(
            "a_%s" % name, "alpha(%s)" % name, 1.8, 1.5, 3.0)

        self.alpha .setConstant(True)
        self.n     .setConstant(True)

        if fixMass:
            self.mean .setConstant(True)
        if fixSigma:
            self.sigma.setConstant(True)

        self.pdf = ROOT.RooCBShape("pdf_%s" % name, "CrystalBall(%s)" % name,
                                   self.m,
                                   self.mean,
                                   self.sigma,
                                   self.alpha,
                                   self.n)

    # fit
    def fitTo(self, dataset, *args, **kwargs):
            """
            Perform the fit
            """
            return self.pdf.fitTo(dataset, *args, **kwargs)


# =============================================================================
# @class Charm_2
#  the final full PDF for 2xCharm fit
#  @author Vanya BELYAEV Ivan.Belyaev@Cern.ch
#  @date 2011-07-25
class Charm2_pdf (object):

    """
    """

    def __init__(self,
                 sig_1,
                 sig_2,
                 suffix=''):

        self.sig1 = sig_1
        self.sig2 = sig_2

        self.m1 = sig_1.mass
        self.m2 = sig_2.mass

        self.bkg1 = Bkg_pdf('Bkg(1)' + suffix, self.m1, power=4)
        self.bkg2 = Bkg_pdf('Bkg(2)' + suffix, self.m2, power=4)

        self.bkgA = Bkg_pdf('Bkg(A)' + suffix, self.m1)
        self.bkgB = Bkg_pdf('Bkg(B)' + suffix, self.m2)

        self.ss_pdf = ROOT.RooProdPdf(
            "S1S2pdf", "Sig(1) x Sig(2)", self.sig1.pdf, self.sig2.pdf)
        self.sb_pdf = ROOT.RooProdPdf(
            "S1B2pdf", "Sig(1) x Bkg(2)", self.sig1.pdf, self.bkg2.pdf)
        self.bs_pdf = ROOT.RooProdPdf(
            "B1S2pdf", "Bkg(1) x Sig(2)", self.bkg1.pdf, self.sig2.pdf)
        self.bb_pdf = ROOT.RooProdPdf(
            "B1B2pdf", "Bkg(1) x Bkg(2)", self.bkg1.pdf, self.bkg2.pdf)

        self.ss = ROOT.RooRealVar(
            "S1S2", "Sig(1) & Sig(2)", 1000,  0,  1.e+5)
        self.sb = ROOT.RooRealVar(
            "S1B2", "Sig(1) & Bkg(2)",   50,  0,  1.e+5)
        self.bs = ROOT.RooRealVar(
            "B1S2", "Bkg(1) & Sig(2)",   50,  0,  1.e+5)
        self.bb = ROOT.RooRealVar(
            "B1B2", "Bkg(A) & Bkg(B)",    5,  0,  1.e+5)

        self.alist1 = ROOT.RooArgList(
            self.ss_pdf,
            self.sb_pdf,
            self.bs_pdf,
            self.bb_pdf)
        self.alist2 = ROOT.RooArgList(
            self.ss,
            self.sb,
            self.bs,
            self.bb)

        self.pdf = ROOT.RooAddPdf("model" + suffix,
                                  "model()%s" % suffix,
                                  self.alist1,
                                  self.alist2)

        self._splots = []

    # fit
    def fitTo(self, dataset, *args):
            """
            Perform the fit
            """
            return self.pdf.fitTo(dataset,
                                  ROOT.RooFit.Save(),
                                  ncpu(len(dataset)),
                                  *args)

    def draw(self, drawvar, dataset=None, nbins=100):

        frame = drawvar.frame(nbins)

        if dataset:
            dataset  .plotOn(frame)

        self.pdf .plotOn(frame,
                         ROOT.RooFit.Components(self.sb_pdf.GetName()),
                         ROOT.RooFit.LineStyle(ROOT.kDashed),
                         ROOT.RooFit.LineColor(ROOT.kGreen))

        self.pdf .plotOn(frame,
                         ROOT.RooFit.Components(self.bs_pdf.GetName()),
                         ROOT.RooFit.LineStyle(ROOT.kDotted),
                         ROOT.RooFit.LineColor(ROOT.kMagenta))

        self.pdf .plotOn(frame,
                         ROOT.RooFit.Components(self.bb_pdf.GetName()),
                         ROOT.RooFit.LineWidth(1),
                         ROOT.RooFit.LineColor(ROOT.kBlack))

        self.pdf .plotOn(frame,
                         ROOT.RooFit.Components(self.ss_pdf.GetName()),
                         ROOT.RooFit.LineWidth(1),
                         ROOT.RooFit.LineColor(ROOT.kRed))

        self.pdf .plotOn(frame,
                         ROOT.RooFit.LineColor(ROOT.kRed))

        frame.SetXTitle('')
        frame.SetYTitle('')
        frame.SetZTitle('')

        frame.Draw()

        return frame

    #
    def sPlot(self,
              dataset,
              *args):

        splot = ROOT.RooStats.SPlot(rootID("sPlot_"),
                                    "sPlot",
                                    dataset,
                                    self.pdf,
                                    self.alist2)

        self._splots += [splot]

        return splot



# =============================================================================
# @class Charm_1
#  the final full PDF for 1xCharm fit
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date 2011-08-02
class Charm1_pdf (object):

    """
    """

    def __init__(self,
                 signal,
                 background=None,
                 signal2=None,
                 suffix=''):

        self.signal = signal
        self.mass = self.signal.mass
        #
        if background:
            self.background = background
        else:
            self.background = Bkg_pdf('Background' + suffix, self.mass)
        #
        self.s = ROOT.RooRealVar(
            "S" + suffix, "Signal", 1000,  0,  1.e+8)
        self.b = ROOT.RooRealVar(
            "B" + suffix, "Background",   10,  0,  1.e+8)
        #
        self.s_name = self.s.GetName()
        self.b_name = self.b.GetName()
        #
        if not signal2:

            self.alist1 = ROOT.RooArgList(
                self.signal     .pdf,
                self.background .pdf
            )
            self.alist2 = ROOT.RooArgList(
                self.s,
                self.b
            )
        else:

            self.signal2 = signal2
            self.s2 = ROOT.RooRealVar(
                "S2" + suffix, "Signal(2)", 1000,  0,  1.e+9)

            self.alist1 = ROOT.RooArgList(
                self.signal     .pdf,
                self.signal2    .pdf,
                self.background .pdf
            )
            self.alist2 = ROOT.RooArgList(
                self.s,
                self.s2,
                self.b
            )

        self.pdf = ROOT.RooAddPdf("model" + suffix,
                                  "model(%s)" % suffix,
                                  self.alist1,
                                  self.alist2)

        self._splots = []

    # fit
    def fitTo(self, dataset, draw=False, nbins=100, refit=True, *args):
        """
        Perform the fit
        """
        size = len(dataset)
        if 1 > size:
            logger.error('Empty data set!')

        nmax = max(1.2 * (size + 10 * math.sqrt(size)), 100)

        for s in self.alist2:
            s.setMin(0.0)
            s.setMax(nmax)

        result = self.pdf.fitTo(dataset,
                                ROOT.RooFit.Save(),
                                ## ncpu ( len ( dataset ) ) ,
                                *args)

        if 0 != result.status():
            result = self.pdf.fitTo(dataset,
                                    ROOT.RooFit.Save(),
                                    ## ncpu ( len ( dataset ) ) ,
                                    *args)

        if draw:
            frame = self.mass.frame(nbins)
            dataset  .plotOn(frame, ROOT.RooFit.Name("data"))
            self.pdf .plotOn(frame,
                             ROOT.RooFit.Components(
                                 self.background.pdf.GetName()),
                             ROOT.RooFit.Name("background"),
                             ROOT.RooFit.LineStyle(ROOT.kDashed),
                             ROOT.RooFit.LineColor(ROOT.kBlue))

            self.pdf .plotOn(frame,
                             ROOT.RooFit.Components(
                                 self.signal.pdf.GetName()),
                             ROOT.RooFit.Name("signal"),
                             ROOT.RooFit.LineStyle(ROOT.kDashed),
                             ROOT.RooFit.LineColor(ROOT.kGreen))

            if hasattr(self, 'signal2'):
                self.pdf .plotOn(frame,
                                 ROOT.RooFit.Components(
                                     self.signal2.pdf.GetName()),
                                 ROOT.RooFit.Name("signal2"),
                                 ROOT.RooFit.LineStyle(ROOT.kDashed),
                                 ROOT.RooFit.LineColor(ROOT.kGreen))

            self.pdf .plotOn(frame,
                             ROOT.RooFit.Name("total"),
                             ROOT.RooFit.LineColor(ROOT.kRed))

            self.legend = ROOT.TLegend(0.65, 0.73, 0.86, 0.87)
            self.legend.SetFillColor(ROOT.kWhite)
            
            self.legend.AddEntry("data", "Data", "P")
            self.legend.AddEntry("total", "Signal + Background", "L")

            self.legend.AddEntry("signal", "Signal", "L")
            self.legend.AddEntry("background", "Background", "L")

            frame.SetXTitle('#Inv.\,mass(J/\psi\,K^+\,K^-\,\pi^+), GeV/c^2')

            frame.SetYTitle('')
            frame.SetZTitle('')

            frame.Draw()
            self.legend.Draw()
            return result, frame

        return result, None

    #
    def sPlot(self,
              dataset,
              *args):

        splot = ROOT.RooStats.SPlot(rootID("sPlot_"),
                                    "sPlot",
                                    dataset,
                                    self.pdf,
                                    self.alist2)

        self._splots += [splot]

        return splot

    def fitHisto(self, histo, draw=False, refit=True, *args):
        """
        Fit the histogram
        """
        self.lst = ROOT.RooArgList(self.mass)
        self.imp = ROOT.RooFit.Import(histo)

        self.hset = ROOT.RooDataHist(
            rootID('hds_'),
            "Data set for histogram '%s'" % histo.GetTitle(),
            self.lst,
            self.imp
        )

        result = self.pdf.fitTo(self.hset,
                                ROOT.RooFit.Save(),
                                *args)

        frame = None

        if draw:

            frame = self.mass.frame()
            self.hset.plotOn(frame, ROOT.RooFit.Name("data"),)

            self.pdf .plotOn(frame,
                             ROOT.RooFit.Components(
                                 self.background.pdf.GetName()),
                             ROOT.RooFit.Name("background"),
                             ROOT.RooFit.LineStyle(ROOT.kDashed),
                             ROOT.RooFit.LineColor(ROOT.kBlue))

            self.pdf .plotOn(frame,
                             ROOT.RooFit.Components(
                                 self.signal.pdf.GetName()),
                             ROOT.RooFit.Name("signal"),
                             ROOT.RooFit.LineStyle(ROOT.kDashed),
                             ROOT.RooFit.LineColor(ROOT.kGreen))

            if hasattr(self, 'signal2'):
                self.pdf .plotOn(frame,
                                 ROOT.RooFit.Components(
                                     self.signal2.pdf.GetName()),
                                 ROOT.RooFit.Name("signal2"),
                                 ROOT.RooFit.LineStyle(ROOT.kDashed),
                                 ROOT.RooFit.LineColor(ROOT.kGreen))

            self.pdf .plotOn(frame, ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name("total"),)

            self.legend = ROOT.TLegend(0.65, 0.73, 0.86, 0.87)
            self.legend.SetFillColor(ROOT.kWhite)
            
            self.legend.AddEntry("data", "Data", "P")
            self.legend.AddEntry("total", "Signal + Background", "L")

            self.legend.AddEntry("signal", "Signal", "L")
            self.legend.AddEntry("background", "Background", "L")

            frame.SetXTitle('Inv.\,Mass(KK),\,GeV/c^2')
            frame.SetYTitle('')
            frame.SetZTitle('')

            frame.Draw()

        return result, frame



# =============================================================================
# @class Charm_3
#  PDF for fitting signal+background + 2 mirrors
#  @author Alexander BARANOV a.baranov@cern.ch
#  @date 2013-09-14
class Charm3_pdf (object):

    """
    """

    def __init__(self,
                 signal,
                 background=None,
                 signal2=None,
                 signal3=None,
                 suffix=''):

        self.signal = signal
        self.mass = self.signal.mass

        if background:
            self.background = background
        else:
            self.background = Bkg_pdf('Background' + suffix, self.mass)

        self.s = ROOT.RooRealVar(
            "S" + suffix, "Signal", 1000,  0,  1.e+8)
        self.b = ROOT.RooRealVar(
            "B" + suffix, "Background",   10,  0,  1.e+8)

        self.s_name = self.s.GetName()
        self.b_name = self.b.GetName()

        if not signal2:
            self.alist1 = ROOT.RooArgList(self.signal     .pdf, self.background .pdf)
            self.alist2 = ROOT.RooArgList(self.s, self.b)
        else:
            self.signal2 = signal2
            self.s2 = ROOT.RooRealVar(
                "S2" + suffix, "Signal(2)", 1000,  0,  1.e+9)

            if not signal3:
                self.alist1 = ROOT.RooArgList(self.signal.pdf, self.signal2.pdf, self.background.pdf)
                self.alist2 = ROOT.RooArgList(self.s, self.s2, self.b)
            else:
                self.signal3 = signal3
                self.s3 = ROOT.RooRealVar("S3" + suffix, "Signal(3)", 1000, 0, 1.e+9)

                self.alist1 = ROOT.RooArgList(
                    self.signal.pdf, 
                    self.signal2,
                    self.signal3,
                    self.background.pdf
                )
                self.alist2 = ROOT.RooArgList(
                    self.s,
                    self.s2,
                    self.s3,
                    self.b
                )



        self.pdf = ROOT.RooAddPdf("model" + suffix,
                                  "model(%s)" % suffix,
                                  self.alist1,
                                  self.alist2)

        self._splots = []

    # fit
    def fitTo(self, dataset, draw=False, nbins=100, refit=True, *args):
        """
        Perform the fit
        """
        size = len(dataset)
        if 1 > size:
            logger.error('Empty data set!')

        nmax = max(1.2 * (size + 10 * math.sqrt(size)), 100)

        for s in self.alist2:
            s.setMin(0.0)
            s.setMax(nmax)

        result = self.pdf.fitTo(dataset,
                                ROOT.RooFit.Save(),
                                ## ncpu ( len ( dataset ) ) ,
                                *args)

        if 0 != result.status():
            result = self.pdf.fitTo(dataset,
                                    ROOT.RooFit.Save(),
                                    ## ncpu ( len ( dataset ) ) ,
                                    *args)

        if draw:
            self.frame = self.mass.frame(nbins)
            dataset  .plotOn(self.frame, ROOT.RooFit.Name("data"))
            self.pdf .plotOn(self.frame,
                             ROOT.RooFit.Components(
                                 self.background.pdf.GetName()),
                             ROOT.RooFit.Name("background"),
                             ROOT.RooFit.LineStyle(ROOT.kDashed),
                             ROOT.RooFit.LineColor(ROOT.kBlue))

            self.pdf .plotOn(self.frame,
                             ROOT.RooFit.Components(
                                 self.signal.pdf.GetName()),
                             ROOT.RooFit.Name("signal"),
                             ROOT.RooFit.LineStyle(ROOT.kDashed),
                             ROOT.RooFit.LineColor(ROOT.kGreen))

            if hasattr(self, 'signal2'):
                self.pdf .plotOn(self.frame,
                                 ROOT.RooFit.Components(
                                     self.signal2.GetName()),
                                 ROOT.RooFit.Name("signal2"),
                                 ROOT.RooFit.LineStyle(ROOT.kDashed),
                                 ROOT.RooFit.LineColor(ROOT.kCyan))
                if hasattr(self, 'signal3'):
                    self.pdf .plotOn(self.frame,
                                 ROOT.RooFit.Components(
                                     self.signal3.GetName()),
                                 ROOT.RooFit.Name("signal3"),
                                 ROOT.RooFit.LineStyle(ROOT.kDashed),
                                 ROOT.RooFit.LineColor(ROOT.kMagenta))


            self.pdf.plotOn(self.frame, ROOT.RooFit.Name("total"), ROOT.RooFit.LineColor(ROOT.kRed))

            self.legend = ROOT.TLegend(0.65, 0.73, 0.9, 0.87);
            self.legend.SetFillColor(ROOT.kWhite);
            
            self.legend.AddEntry("data", "Data", "P")
            self.legend.AddEntry("total", "Signal + Background + Reflections", "L")

            self.legend.AddEntry("signal", "Signal", "L")
            self.legend.AddEntry("background", "Exponential background", "L")

            self.legend.AddEntry("signal2", "B^+ \\to J/\psi K^+ K^- K^+ reflection", "L")
            self.legend.AddEntry("signal3", "B^+ \\to J/\psi K^+ \pi^+ \pi^- reflection", "L")

            self.frame.SetXTitle('#Inv.\,mass(J/\psi\,K^+\,K^-\,\pi^+), GeV/c^2')
            self.frame.SetYTitle('')
            self.frame.SetZTitle('')

            self.frame.Draw()
            self.legend.Draw()

            return result, self.frame

        return result, None

    #
    def sPlot(self,
              dataset,
              *args):

        splot = ROOT.RooStats.SPlot(rootID("sPlot_"),
                                    "sPlot",
                                    dataset,
                                    self.pdf,
                                    self.alist2)

        self._splots += [splot]

        return splot

    def fitHisto(self, histo, draw=False, refit=True, *args):
        """
        Fit the histogram
        """
        self.lst = ROOT.RooArgList(self.mass)
        self.imp = ROOT.RooFit.Import(histo)

        self.hset = ROOT.RooDataHist(
            rootID('hds_'),
            "Data set for histogram '%s'" % histo.GetTitle(),
            self.lst,
            self.imp
        )

        result = self.pdf.fitTo(self.hset,
                                ROOT.RooFit.Save(),
                                *args)

        self.frame = None

        if draw:
            self.frame = self.mass.frame()
            self.hset.plotOn(self.frame)

            self.pdf .plotOn(self.frame,
                             ROOT.RooFit.Components(
                                 self.background.pdf.GetName()),
                             ROOT.RooFit.Name("background"),
                             ROOT.RooFit.LineStyle(ROOT.kDashed),
                             ROOT.RooFit.LineColor(ROOT.kBlue))

            self.pdf .plotOn(self.frame,
                             ROOT.RooFit.Components(
                                 self.signal.pdf.GetName()),
                             ROOT.RooFit.Name("signal"),
                             ROOT.RooFit.LineStyle(ROOT.kDashed),
                             ROOT.RooFit.LineColor(ROOT.kGreen))

            if hasattr(self, 'signal2'):
                self.pdf .plotOn(self.frame,
                                 ROOT.RooFit.Components(
                                     self.signal2.GetName()),
                                 ROOT.RooFit.LineStyle(ROOT.kDashed),
                                 ROOT.RooFit.LineColor(ROOT.kPink))
                if hasattr(self, 'signal3'):
                    self.pdf .plotOn(self.frame,
                                 ROOT.RooFit.Components(
                                     self.signal3.GetName()),
                                 ROOT.RooFit.LineStyle(ROOT.kDashed),
                                 ROOT.RooFit.LineColor(ROOT.kMagenta))



            self.pdf .plotOn(self.frame, ROOT.RooFit.LineColor(ROOT.kRed))

            self.frame.SetXTitle('')
            self.frame.SetYTitle('')
            self.frame.SetZTitle('')

            self.frame.Draw()

        return result, self.frame


# =============================================================================
# @class One
#  the final full PDF for 1xCharm fit
#  @author Vanya BELYAEV Ivan.Belyaeve@cern.ch
#  @date 2011-08-02
class One_pdf (object):

    """
    """

    def __init__(self,
                 signal,
                 background=None,
                 components=[],  # components
                 histos=[],  # histos as components
                 suffix=''):

        self.signal = signal
        self.mass = self.signal.mass
        #
        if background:
            self.background = background
        else:
            self.background = Bkg_pdf('Background', self.mass)
        #
        self.s = ROOT.RooRealVar(
            "S" + suffix, "Signal", 1000,  0,  1.e+6)
        self.b = ROOT.RooRealVar(
            "B" + suffix, "Background",   10,  0,  1.e+6)

        self.alist1 = ROOT.RooArgList(
            self.signal     .pdf,
            self.background .pdf
        )

        self.alist2 = ROOT.RooArgList(
            self.s,
            self.b
        )

        self.nums = [self.s, self.b]

        i = 1
        for c in components:
            si = ROOT.RooRealVar(
                "S_%d" % i + suffix, "Signal", 10,  0,  1.e+6)

            self.alist2.add(si)
            self.alist1.add(c)

            setattr(self, si.GetName(), si)

            self.nums.append(si)
            i += 1

        self._histos = histos

        self._h_lsts = []
        self._h_imps = []
        self._h_hset = []
        self._h_sets = []
        self._h_pdfs = []
        self._h_vars = []

        for histo in self._histos:

            h_lst = ROOT.RooArgList(self.mass)
            h_imp = ROOT.RooFit.Import(histo)
            h_hset = ROOT.RooDataHist(
                rootD('hds_'),
                "Data set for histogram '%s'" % histo.GetTitle(),
                h_lst,
                h_imp)
            h_set = ROOT.RooArgSet(self.mass)
            h_pdf = ROOT.RooHistPdf(
                'hpdf%d' % i + suffix,
                'h-pdf',  # + histo.GetName() + '/' + histo.GetTitle() ,
                h_set,
                h_hset)

            self._h_lsts.append(h_lst)
            self._h_imps.append(h_imp)
            self._h_hset.append(h_hset)
            self._h_sets.append(h_set)
            self._h_pdfs.append(h_pdf)

            hi = ROOT.RooRealVar(
                "H_%d" % i + suffix, "Signal(%d)" % i, 10,  0,  1e+6)
            self._h_vars.append(hi)

            setattr(self, hi.GetName(), hi)

            self.nums.append(hi)

            self.alist2.add(hi)
            self.alist1.add(h_pdf)

            i += 1

        self.pdf = ROOT.RooAddPdf("model",
                                  "model",
                                  self.alist1,
                                  self.alist2)

        self._splots = []

    # fit
    def fitTo(self, dataset, draw=False, nbins=100, *args):
        """
        Perform the fit
        """
        result = self.pdf.fitTo(dataset,
                                ROOT.RooFit.Save(),
                                ncpu(len(dataset)),
                                *args)
        if draw:
            frame = self.mass.frame(nbins)
            dataset  .plotOn(frame)
            self.pdf .plotOn(frame,
                             ROOT.RooFit.Components(
                                 self.background.pdf.GetName()),
                             ROOT.RooFit.LineStyle(ROOT.kDashed),
                             ROOT.RooFit.LineColor(ROOT.kBlue))
            self.pdf .plotOn(frame,
                             ROOT.RooFit.LineColor(ROOT.kRed))

            frame.SetXTitle('')
            frame.SetYTitle('')
            frame.SetZTitle('')

            frame.Draw()
            return result, frame

        return result, None

    #
    def sPlot(self,
              dataset,
              *args):

        splot = ROOT.RooStats.SPlot(rootID("sPlot_"),
                                    "sPlot",
                                    dataset,
                                    self.pdf,
                                    self.alist2)

        self._splots += [splot]

        return splot

    def fitHisto(self, histo, draw=False, *args):
        """
        Fit the histogram
        """
        self.lst = ROOT.RooArgList(self.mass)
        self.imp = ROOT.RooFit.Import(histo)
        self.hset = ROOT.RooDataHist(
            rootID('hds_'),
            "Data set for histogram '%s'" % histo.GetTitle(),
            self.lst,
            self.imp
        )

        result = self.pdf.fitTo(self.hset,
                                ROOT.RooFit.Save(),
                                *args)

        if draw:

            frame = self.mass.frame()
            self.hset.plotOn(frame)

            self.pdf .plotOn(frame,
                             ROOT.RooFit.Components(
                                 self.background.pdf.GetName()),
                             ROOT.RooFit.LineStyle(ROOT.kDashed),
                             ROOT.RooFit.LineColor(ROOT.kBlue))

            self.pdf .plotOn(frame,
                             ROOT.RooFit.Components(
                                 self.signal.pdf.GetName()),
                             ROOT.RooFit.LineStyle(ROOT.kDashed),
                             ROOT.RooFit.LineColor(ROOT.kGreen))

            self.pdf .plotOn(frame,
                             ROOT.RooFit.LineColor(ROOT.kRed))

            frame.SetXTitle('')
            frame.SetYTitle('')
            frame.SetZTitle('')

            frame.Draw()
            return result, frame

        return result, None


# =============================================================================
# @class
#  the final full PDF for 1xCharm
#  @author Vanya BELYAEV Ivan.Belyaeve@cern.ch
#  @date 2011-08-02
class Bc_pdf (object):

    """
    """

    def __init__(self,
                 signal,
                 background=None,
                 suffix='',
                 hist1=None,
                 hist2=None,
                 frac=False,
                 hist3=None,
                 hist4=None):

        self.signal = signal
        self.mass = self.signal.mass

        if background:
            self.background = background
        else:
            self.background = Bkg_pdf('Background', self.mass)

        self.s = ROOT.RooRealVar("S" + suffix, "Signal",  20,  0,  1.e+6)
        self.b = ROOT.RooRealVar(
            "B" + suffix, "Background",  50,  0,  1.e+6)

        self.alist1 = ROOT.RooArgList(
            self.signal     .pdf,
            self.background .pdf
        )
        self.alist2 = ROOT.RooArgList(
            self.s,
            self.b
        )

        if hist1:

            self.h1_lst = ROOT.RooArgList(self.mass)
            self.h1_imp = ROOT.RooFit.Import(hist1)
            self.h1_hset = ROOT.RooDataHist(
                rootID('hds_'),
                "Data set for histogram '%s'" % hist1.GetTitle(),
                self.h1_lst,
                self.h1_imp
            )
            self.h1_set = ROOT.RooArgSet(self.mass)
            self.h1_pdf = ROOT.RooHistPdf('h1_pdf' + suffix,
                                          'h1_pdf',
                                          self.h1_set,
                                          self.h1_hset)

        if hist2:

            self.h2_lst = ROOT.RooArgList(self.mass)
            self.h2_imp = ROOT.RooFit.Import(hist2)
            self.h2_hset = ROOT.RooDataHist(
                rootID('hds_'),
                "Data set for histogram '%s'" % hist2.GetTitle(),
                self.h2_lst,
                self.h2_imp
            )
            self.h2_set = ROOT.RooArgSet(self.mass)
            self.h2_pdf = ROOT.RooHistPdf('h2_pdf' + suffix,
                                          'h2_pdf',
                                          self.h2_set,
                                          self.h2_hset)
        if hist3:

            self.h3_lst = ROOT.RooArgList(self.mass)
            self.h3_imp = ROOT.RooFit.Import(hist3)
            self.h3_hset = ROOT.RooDataHist(
                rootID('hds_'),
                "Data set for histogram '%s'" % hist3.GetTitle(),
                self.h3_lst,
                self.h3_imp
            )
            self.h3_set = ROOT.RooArgSet(self.mass)
            self.h3_pdf = ROOT.RooHistPdf('h3_pdf' + suffix,
                                          'h3_pdf',
                                          self.h3_set,
                                          self.h3_hset)
        if hist4:

            self.h4_lst = ROOT.RooArgList(self.mass)
            self.h4_imp = ROOT.RooFit.Import(hist4)
            self.h4_hset = ROOT.RooDataHist(
                rootID('hds_'),
                "Data set for histogram '%s'" % hist4.GetTitle(),
                self.h4_lst,
                self.h4_imp
            )
            self.h4_set = ROOT.RooArgSet(self.mass)
            self.h4_pdf = ROOT.RooHistPdf('h4_pdf' + suffix,
                                          'h4_pdf',
                                          self.h4_set,
                                          self.h4_hset)

        if frac and hist1 and hist2:

            self.sb = ROOT.RooRealVar(
                "SB" + suffix, "Signal/SB",  20,  0,  1.e+6)
            self.fr = ROOT.RooRealVar(
                "FR" + suffix, "Signal/FR",  0.9,  0,  1)
            self.frlst = ROOT.RooArgList(self.sb, self.fr)

            self.s1 = ROOT.RooFormulaVar("S1_%s" % suffix,
                                         "S1(%s)" % suffix,
                                         "%s*%s" % (
                                             self.sb.GetName(
                                             ), self.fr.GetName()),
                                         self.frlst)

            self.s2 = ROOT.RooFormulaVar("S2_%s" % suffix,
                                         "S2(%s)" % suffix,
                                         "%s*(1-%s)" % (
                                             self.sb.GetName(
                                             ), self.fr.GetName()),
                                         self.frlst)

            self.alist1.add(self.h1_pdf)
            self.alist2.add(self.s1)
            self.alist1.add(self.h2_pdf)
            self.alist2.add(self.s2)
        else:

            if hist1:
                self.s1 = ROOT.RooRealVar(
                    "S1" + suffix, "Signal/H1",  1,  0,  1.e+6)
                self.alist1.add(self.h1_pdf)
                self.alist2.add(self.s1)

            if hist2:
                self.s2 = ROOT.RooRealVar(
                    "S2" + suffix, "Signal/H2",  20,  0,  1.e+6)
                self.alist1.add(self.h2_pdf)
                self.alist2.add(self.s2)

        if frac and hist3 and hist4:

            self.sb_2 = ROOT.RooRealVar(
                "SB2" + suffix, "Signal/SB",  20,  0,  1.e+6)
            self.fr_2 = ROOT.RooRealVar(
                "FR2" + suffix, "Signal/FR",  0.9,  0,  1)
            self.frlst_2 = ROOT.RooArgList(self.sb_2, self.fr_2)

            self.s3 = ROOT.RooFormulaVar("S3_%s" % suffix,
                                         "S3(%s)" % suffix,
                                         "%s*%s" % (
                                             self.sb_2.GetName(
                                             ), self.fr_2.GetName()),
                                         self.frlst_2)

            self.s4 = ROOT.RooFormulaVar("S4_%s" % suffix,
                                         "S4(%s)" % suffix,
                                         "%s*(1-%s)" % (
                                             self.sb_2.GetName(
                                             ), self.fr_2.GetName()),
                                         self.frlst_2)

            self.alist1.add(self.h3_pdf)
            self.alist2.add(self.s3)
            self.alist1.add(self.h4_pdf)
            self.alist2.add(self.s4)
        else:

            if hist3:
                self.s3 = ROOT.RooRealVar(
                    "S3" + suffix, "Signal/H3",  20,  0,  1.e+6)
                self.alist1.add(self.h3_pdf)
                self.alist2.add(self.s3)

            if hist4:
                self.s4 = ROOT.RooRealVar(
                    "S4" + suffix, "Signal/H4",  20,  0,  1.e+6)
                self.alist1.add(self.h4_pdf)
                self.alist2.add(self.s4)

        self.pdf = ROOT.RooAddPdf("model",
                                  "model",
                                  self.alist1,
                                  self.alist2)

        self._splots = []

    # fit
    def fitTo(self, dataset, draw=False, nbins=100, *args):
        """
        Perform the fit
        """
        result = self.pdf.fitTo(dataset,
                                ROOT.RooFit.Save(),
                                ncpu(len(dataset)),
                                *args)

        if draw:
            frame = self.mass.frame(nbins)
            dataset  .plotOn(frame)

            self.pdf .plotOn(frame,
                             ROOT.RooFit.Components(
                                 self.background.pdf.GetName()),
                             ## ROOT.RooFit.LineWidth  ( 1 ) ,
                             ROOT.RooFit.LineStyle(ROOT.kDashed),
                             ROOT.RooFit.LineColor(ROOT.kBlue))

            if hasattr(self, 'h1_pdf') and self.s1.getVal() > 0.1:
                self.pdf .plotOn(frame,
                                 ROOT.RooFit.Components(
                                     self.h1_pdf.GetName()),
                                 ## ROOT.RooFit.LineWidth  ( 1         ) ,
                                 ROOT.RooFit.LineStyle(ROOT.kDotted),
                                 ROOT.RooFit.LineColor(ROOT.kGreen + 1))

            if hasattr(self, 'h2_pdf') and self.s2.getVal() > 0.1:
                self.pdf .plotOn(frame,
                                 ROOT.RooFit.Components(
                                     self.h2_pdf.GetName()),
                                 ## ROOT.RooFit.LineWidth  ( 1         ) ,
                                 ROOT.RooFit.LineStyle(ROOT.kDotted),
                                 ROOT.RooFit.LineColor(ROOT.kOrange))

            if hasattr(self, 'h3_pdf') and self.s3.getVal() > 0.1:
                self.pdf .plotOn(frame,
                                 ROOT.RooFit.Components(
                                     self.h3_pdf.GetName()),
                                 ROOT.RooFit.LineStyle(ROOT.kDashDotted),
                                 ## ROOT.RooFit.LineWidth  ( 3 ) ,
                                 ROOT.RooFit.LineColor(ROOT.kMagenta))

            if hasattr(self, 'h4_pdf') and self.s4.getVal() > 0.1:
                self.pdf .plotOn(frame,
                                 ROOT.RooFit.Components(
                                     self.h4_pdf.GetName()),
                                 ROOT.RooFit.LineStyle(ROOT.kDashDotted),
                                 ## ROOT.RooFit.LineWidth  ( 3 ) ,
                                 ROOT.RooFit.LineColor(ROOT.kCyan))

            self.pdf .plotOn(frame,
                             ## ROOT.RooFit.LineWidth  ( 1         ) ,
                             ROOT.RooFit.LineColor(ROOT.kRed))

            frame.SetXTitle('')
            frame.SetYTitle('')
            frame.SetZTitle('')

            frame.Draw()
            return result, frame

        return result, None

    #
    def sPlot(self,
              dataset,
              *args):

        splot = ROOT.RooStats.SPlot(rootID("sPlot_"),
                                    "sPlot",
                                    dataset,
                                    self.pdf,
                                    self.alist2)

        self._splots += [splot]

        return splot

    def fitHisto(self, histo, draw=False, *args):
        """
        Fit the histogram
        """
        self.lst = ROOT.RooArgList(self.mass)
        self.imp = ROOT.RooFit.Import(histo)
        self.hset = ROOT.RooDataHist(
            rootID('hds_'),
            "Data set for histogram '%s'" % histo.GetTitle(),
            self.lst,
            self.imp
        )

        result = self.pdf.fitTo(self.hset,
                                ROOT.RooFit.Save(),
                                *args)

        if draw:
            frame = self.mass.frame()
            self.hset.plotOn(frame)

            if 1.e-3 < self.b.getVal():
                self.pdf .plotOn(frame,
                                 ROOT.RooFit.Components(
                                     self.background.pdf.GetName()),
                                 ROOT.RooFit.LineStyle(ROOT.kDashed),
                                 ROOT.RooFit.LineColor(ROOT.kBlue))

            if hasattr(self, 'h1_pdf'):
                self.pdf .plotOn(frame,
                                 ROOT.RooFit.Components(
                                     self.h1_pdf.GetName()),
                                 ROOT.RooFit.LineStyle(ROOT.kDotted),
                                 ROOT.RooFit.LineColor(ROOT.kGreen))

            if hasattr(self, 'h2_pdf'):
                self.pdf .plotOn(frame,
                                 ROOT.RooFit.Components(
                                     self.h2_pdf.GetName()),
                                 ROOT.RooFit.LineStyle(ROOT.kDotted),
                                 ROOT.RooFit.LineColor(ROOT.kYellow))

            self.pdf .plotOn(frame,
                             ROOT.RooFit.LineColor(ROOT.kRed))

            frame.SetXTitle('')
            frame.SetYTitle('')
            frame.SetZTitle('')

            frame.Draw()

            return result, frame

        return result, None


# =============================================================================
# @class UpsFit
#  the final full PDF for Upsilons fit
#  @author Vanya BELYAEV Ivan.Belyaeve@cern.ch
#  @date 2011-08-02
class UpsFit_pdf (object):

    """
    """

    def __init__(self,
                 name='Y',
                 mass=None,
                 gev=True,
                 power=0):

        mn = 8.5
        mx = 12.5
        m_y1s = 9.4647e+00
        s_y1s = 4.3679e-02
        dm_y2s = 10.023 - m_y1s
        dm_y3s = 10.355 - m_y1s

        if not gev:
            mn *= 1000
            mx *= 1000
            m_y1s *= 1000
            s_y1s *= 1000
            dm_y2s *= 1000
            dm_y3s *= 1000

        if mass:
            self.mass = mass
        else:
            self.mass = ROOT.RooRealVar(
                "m_%s" % name, "mass(%s)" % name,       mn, mx)

        #
        self.m1s = ROOT.RooRealVar("m1s", "mean (%s)" %
                                   name, m_y1s, m_y1s - 0.2 * s_y1s, m_y1s + 0.2 * s_y1s)
        self.sigma = ROOT.RooRealVar("sigma", "sigma(%s)" %
                                     name, s_y1s,         0.5 * s_y1s,         2.0 * s_y1s)

        self.n = ROOT.RooRealVar(
            "n_%s" % name, "n(%s)" % name, 1, 1, 10)
        self.alpha = ROOT.RooRealVar(
            "a_%s" % name, "n(%s)" % name, 1.5, 0.5, 5.0)

        self.m = self.mass

        self.n     .setConstant(True)
        self.alpha .setConstant(True)

        self.y1s = ROOT.RooCBShape("y1s_%s" % name,
                                   "CristalBall(%s)" % name,
                                   self.mass,
                                   self.m1s,
                                   self.sigma,
                                   self.alpha,
                                   self.n)

        # self.y1s   = ROOT.RooGaussian   ( "y1s_%s" % name ,
        ##                                           "CristalBall(%s)" % name ,
        ##                                           self.mass  ,
        ##                                           self.m1s   ,
        # self.sigma )

        #
        # Y(3S)
        #
        self.dm2s = ROOT.RooRealVar("dm2s",
                                    "dm2s(%s)" % name,
                                    dm_y2s,
                                    dm_y2s - 0.20 * s_y1s,
                                    dm_y2s + 0.20 * s_y1s)

        self.aset11 = ROOT.RooArgList(self.m1s, self.dm2s)
        self.m2s = ROOT.RooFormulaVar("m2s",
                                      "m2s(%s)" % name,
                                      "m1s+dm2s",
                                      self.aset11)

        self.aset12 = ROOT.RooArgList(self.sigma, self.m1s, self.m2s)
        self.s2s = ROOT.RooFormulaVar("s2s" + name,
                                      "s2s(%s)" % name,
                                      "sigma*m2s/m1s",
                                      self.aset12)

        self.y2s = ROOT.RooCBShape("y2s_%s" % name,
                                   "CristalBall(%s)" % name,
                                   self.mass,
                                   self.m2s,
                                   self.s2s,
                                   ## self.sigma         ,
                                   self.alpha,
                                   self.n)

        #
        # Y(3S)
        #
        self.dm3s = ROOT.RooRealVar("dm3s",
                                    "dm3s(%s)" % name,
                                    dm_y3s,
                                    dm_y3s - 0.20 * s_y1s,
                                    dm_y3s + 0.20 * s_y1s)

        self.aset21 = ROOT.RooArgList(self.m1s, self.dm3s)
        self.m3s = ROOT.RooFormulaVar("m3s",
                                      "m3s(%s)" % name,
                                      "m1s+dm3s",
                                      self.aset21)

        self.aset22 = ROOT.RooArgList(self.sigma, self.m1s, self.m3s)
        self.s3s = ROOT.RooFormulaVar("s3s" + name,
                                      "s3s(%s)" % name,
                                      "sigma*m3s/m1s",
                                      self.aset22)

        self.y3s = ROOT.RooCBShape("y3s_%s" % name,
                                   "CristalBall(%s)" % name,
                                   self.mass,
                                   self.m3s,
                                   self.s3s,
                                   ## self.sigma ,
                                   self.alpha,
                                   self.n)

        self.background = Bkg_pdf('Bkg%s' % name, self.mass, power=power)

        self.n1s = ROOT.RooRealVar(
            "N1S" + name, "Signal(Y1S)",  100,  0,  1.e+7)
        self.n2s = ROOT.RooRealVar(
            "N2S" + name, "Signal(Y2S)",  100,  0,  1.e+7)
        self.n3s = ROOT.RooRealVar(
            "N3S" + name, "Signal(Y3S)",  100,  0,  1.e+7)
        self.b = ROOT.RooRealVar(
            "B" + name, "Background",   10,  0,  1.e+8)

        self.alist1 = ROOT.RooArgList(self.y1s, self.y2s, self.y3s)
        self.alist2 = ROOT.RooArgList(self.n1s, self.n2s, self.n3s)

        self.alist1 . add(self.background.pdf)
        self.alist2 . add(self.b)

        self.pdf = ROOT.RooAddPdf("model%s" % name,
                                  "model(%s)" % name,
                                  self.alist1,
                                  self.alist2)

        self.bset = ROOT.RooArgSet(self.background.pdf)
        self.yset = ROOT.RooArgSet(self.y1s, self.y2s, self.y3s)

        self._splots = []

    # fit
    def fitTo(self, dataset, draw=False, nbins=100, *args):
        """
        Perform the fit
        """
        result = self.pdf.fitTo(dataset,
                                ROOT.RooFit.Save(),
                                ncpu(len(dataset)),
                                *args)
        ntot = 0
        for i in self.alist2:
            print ' i: %s (%s+-%s) ' % (i.GetName(), i.getVal(), i.getError())
            _ntot = ntot
            try:
                ntot = ntot + result(i.GetName())[0]
            except:
                ntot = _ntot

        result. ntot = ntot

        print ' #TOTAL ', ntot, len(dataset)

        frame = None
        if draw:
            frame = self.mass.frame(nbins)
            dataset  .plotOn(frame)

            self.pdf .plotOn(frame,
                             ROOT.RooFit.Components(self.bset),
                             ROOT.RooFit.LineStyle(ROOT.kDashed),
                             ROOT.RooFit.LineColor(ROOT.kBlue))

            #
            self.pdf .plotOn(frame,
                             ROOT.RooFit.Components(self.yset),
                             ROOT.RooFit.LineWidth(2),
                             ROOT.RooFit.LineStyle(ROOT.kDotted),
                             ROOT.RooFit.LineColor(ROOT.kMagenta))

            self.pdf .plotOn(frame,
                             ROOT.RooFit.LineColor(ROOT.kRed))

            frame.SetXTitle('')
            frame.SetYTitle('')
            frame.SetZTitle('')

            frame.Draw()

        return result, frame
       #

    def sPlot(self,
              dataset,
              *args):

        splot = ROOT.RooStats.SPlot(rootID("sPlot_"),
                                    "sPlot",
                                    dataset,
                                    self.pdf,
                                    self.alist2)

        self._splots += [splot]

        return splot

    def fitHisto(self, histo, draw=False, nbins=100, *args):
        """
        Fit the histogram
        """
        self.lst = ROOT.RooArgList(self.mass)
        self.imp = ROOT.RooFit.Import(histo)
        self.hset = ROOT.RooDataHist(
            rootID('hds_'),
            "Data set for histogram '%s'" % histo.GetTitle(),
            self.lst,
            self.imp
        )

        result = self.pdf.fitTo(self.hset,
                                ROOT.RooFit.Save(),
                                *args)

        if not draw:
            return result, None

        frame = self.mass.frame(nbins)
        self.hset.plotOn(frame)
        self.pdf .plotOn(frame,
                         ROOT.RooFit.Components(self.bset),
                         ROOT.RooFit.LineStyle(ROOT.kDashed),
                         ROOT.RooFit.LineColor(ROOT.kBlue))
        self.pdf .plotOn(frame,
                         ROOT.RooFit.Components(self.yset),
                         ROOT.RooFit.LineWidth(2),
                         ROOT.RooFit.LineStyle(ROOT.kDotted),
                         ROOT.RooFit.LineColor(ROOT.kMagenta))
        self.pdf .plotOn(frame,
                         ROOT.RooFit.LineColor(ROOT.kRed))

        frame.SetXTitle('')
        frame.SetYTitle('')
        frame.SetZTitle('')

        frame.Draw()

        return result, frame


# =============================================================================
# @class Ups2
#  the final full PDF for 2xCharm fit
#  @author Vanya BELYAEV Ivan.Belyaev@Cern.ch
#  @date 2011-07-25
class Ups2_pdf (object):

    """
    """

    def __init__(self,
                 sig_1,
                 ups_2=UpsFit_pdf('Y'),
                 power=0,
                 suffix=''):

        self.sig1 = sig_1
        self.sig2 = ups_2

        self.m1 = sig_1.mass
        self.m2 = ups_2.mass

        ## backgrund + signal
        self.bkg1 = Bkg_pdf(
            'Bkg(B1)' + suffix, self.m1, power=0)
        ## signal    + background
        self.bkg2 = Bkg_pdf(
            'Bkg(B2)' + suffix, self.m2, power)

        # self.bkg1s   = Bkg_pdf        ( 'Bkg(1s)' + suffix , self.m1 )
        # self.bkg2s   = Bkg_pdf        ( 'Bkg(2s)' + suffix , self.m1 )
        # self.bkg3s   = Bkg_pdf        ( 'Bkg(3s)' + suffix , self.m1 )

        self.bkgA = Bkg_pdf(
            'Bkg(1)' + suffix, self.m1, power=0)
        self.bkgB = Bkg_pdf('Bkg(Y)' + suffix, self.m2, power)

        self.sy1_pdf = ROOT.RooProdPdf(
            "S1Y1pdf" + suffix, "Sig(1) x Y(1S)", self.sig1.pdf, self.sig2.y1s)
        self.sy2_pdf = ROOT.RooProdPdf(
            "S1Y2pdf" + suffix, "Sig(1) x Y(2S)", self.sig1.pdf, self.sig2.y2s)
        self.sy3_pdf = ROOT.RooProdPdf(
            "S1Y3pdf" + suffix, "Sig(1) x Y(3S)", self.sig1.pdf, self.sig2.y3s)

        # self.by1_pdf = ROOT.RooProdPdf  ( "B1Y1pdf" + suffix , "Bkg(1) x Y(1s)"   , self.bkg1s.pdf , self.sig2.y1s  )
        # self.by2_pdf = ROOT.RooProdPdf  ( "B1Y2pdf" + suffix , "Bkg(1) x Y(2s)"   , self.bkg2s.pdf , self.sig2.y2s  )
        # self.by3_pdf = ROOT.RooProdPdf  ( "B1Y3pdf" + suffix , "Bkg(1) x Y(3s)"   , self.bkg3s.pdf , self.sig2.y3s  )

        self.by1_pdf = ROOT.RooProdPdf(
            "B1Y1pdf" + suffix, "Bkg(1) x Y(1S)", self.bkg1.pdf, self.sig2.y1s)
        self.by2_pdf = ROOT.RooProdPdf(
            "B1Y2pdf" + suffix, "Bkg(1) x Y(2S)", self.bkg1.pdf, self.sig2.y2s)
        self.by3_pdf = ROOT.RooProdPdf(
            "B1Y3pdf" + suffix, "Bkg(1) x Y(3S)", self.bkg1.pdf, self.sig2.y3s)

        self.sb_pdf = ROOT.RooProdPdf(
            "S1B2pdf" + suffix, "Sig(1) x Bkg(Y)", self.sig1.pdf, self.bkg2.pdf)

        self.bb_pdf = ROOT.RooProdPdf(
            "B1B2pdf" + suffix, "Bkg(1) x Bkg(Y)", self.bkgA.pdf, self.bkgB.pdf)

        self.sy1 = ROOT.RooRealVar(
            "S1Y1" + suffix, "Sig(1) & Y(1S)", 200,  0,  5.e+5)
        self.sy2 = ROOT.RooRealVar(
            "S1Y2" + suffix, "Sig(1) & Y(2S)",  50,  0,  5.e+5)
        self.sy3 = ROOT.RooRealVar(
            "S1Y3" + suffix, "Sig(1) & Y(3S)",  10,  0,  5.e+5)

        self.by1 = ROOT.RooRealVar(
            "B1Y1" + suffix, "Bkg(1) & Y(1S)", 200,  0,  5.e+5)
        self.by2 = ROOT.RooRealVar(
            "B1Y2" + suffix, "Bkg(1) & Y(2S)",  50,  0,  5.e+5)
        self.by3 = ROOT.RooRealVar(
            "B1Y3" + suffix, "Bkg(1) & Y(3S)",  10,  0,  5.e+5)

        self.sb = ROOT.RooRealVar(
            "S1B2" + suffix, "Sig(1) & Bkg(Y)",  50,  0,  5.e+5)
        self.bb = ROOT.RooRealVar(
            "B1B2" + suffix, "Bkg(1) & Bkg(Y)",   5,  0,  5.e+5)

        self.alist1 = ROOT.RooArgList(
            #
            self.sy1_pdf,
            self.sy2_pdf,
            self.sy3_pdf,
            #
            self.by1_pdf,
            self.by2_pdf,
            self.by3_pdf,
            #
            self.sb_pdf,
            self.bb_pdf)
        self.alist2 = ROOT.RooArgList(
            #
            self.sy1,
            self.sy2,
            self.sy3,
            #
            self.by1,
            self.by2,
            self.by3,
            #
            self.sb,
            self.bb)

        self.pdf = ROOT.RooAddPdf("model" + suffix,
                                  "model(%s)" % suffix,
                                  self.alist1,
                                  self.alist2)

        self._splots = []

    # fit
    def fitTo(self, dataset, draw=False, *args):
            """
            Perform the fit
            """
            result = self.pdf.fitTo(dataset,
                                    ROOT.RooFit.Save(),
                                    ncpu(len(dataset)),
                                    *args)

            ntot = 0
            for i in self.alist2:
                print ' i: %s (%s+-%s) %s ' % (i.GetName(), i.getVal(), i.getError(), result(i.GetName())[0])
                ntot = ntot + result(i.GetName())[0]

            result. ntot = ntot

            print ' #TOTAL ', ntot, len(dataset)

            return result

    def draw(self, drawvar, dataset=None, nbins=100, *args):

        frame = drawvar.frame(nbins)

        if dataset:
            dataset  .plotOn(frame, *args)

        if self.m1 == drawvar:

            self.set1 = ROOT.RooArgSet(self.sy1_pdf,
                                       self.sy2_pdf,
                                       self.sy3_pdf,
                                       self.by1_pdf,
                                       self.by2_pdf,
                                       self.by3_pdf)
            self.set11 = ROOT.RooArgSet(self.by1_pdf,
                                        self.by2_pdf,
                                        self.by3_pdf)
            self.set2 = ROOT.RooArgSet(self.sb_pdf,
                                       self.bb_pdf)
            self.set21 = ROOT.RooArgSet(self.bb_pdf)

            self.bset = ROOT.RooArgSet(self.by1_pdf,
                                       self.by2_pdf,
                                       self.by3_pdf,
                                       self.bb_pdf)

        elif self.m2 == drawvar:

            self.set1 = ROOT.RooArgSet(self.sy1_pdf,
                                       self.sy2_pdf,
                                       self.sy3_pdf,
                                       self.sb_pdf)
            self.set11 = ROOT.RooArgSet(self.sb_pdf)

            self.set2 = ROOT.RooArgSet(self.by1_pdf,
                                       self.by2_pdf,
                                       self.by3_pdf,
                                       self.bb_pdf)
            self.set21 = ROOT.RooArgSet(self.bb_pdf)

            self.bset = ROOT.RooArgSet(self.sb_pdf,
                                       self.bb_pdf)

        self.pdf .plotOn(frame,
                         ROOT.RooFit.Components(self.set11),
                         ROOT.RooFit.LineWidth(2),
                         ROOT.RooFit.LineStyle(ROOT.kDashed),
                         ROOT.RooFit.LineColor(ROOT.kMagenta), *args)
        self.pdf .plotOn(frame,
                         ROOT.RooFit.Components(self.set1),
                         ROOT.RooFit.LineWidth(2),
                         ROOT.RooFit.LineColor(ROOT.kMagenta), *args)

        self.pdf .plotOn(frame,
                         ROOT.RooFit.Components(self.set21),
                         ROOT.RooFit.LineWidth(2),
                         ROOT.RooFit.LineStyle(ROOT.kDashed),
                         ROOT.RooFit.LineColor(ROOT.kBlue), *args)

        self.pdf .plotOn(frame,
                         ROOT.RooFit.Components(self.set2),
                         ROOT.RooFit.LineWidth(2),
                         ROOT.RooFit.LineColor(ROOT.kBlue), *args)

        self.pdf .plotOn(frame,
                         ROOT.RooFit.Components(self.bset),
                         ROOT.RooFit.LineStyle(ROOT.kDashed),
                         ROOT.RooFit.LineColor(ROOT.kRed), *args)

        self.pdf .plotOn(frame,
                         ROOT.RooFit.LineColor(ROOT.kRed), *args)

        frame.SetXTitle('')
        frame.SetYTitle('')
        frame.SetZTitle('')

        frame.Draw()

        return frame

    #
    def sPlot(self,
              dataset,
              *args):

        self._splots = []

        splot = ROOT.RooStats.SPlot(rootID("sPlot_"),
                                    "sPlot",
                                    dataset,
                                    self.pdf,
                                    self.alist2)

        self._splots += [splot]

        return splot

    def fitHisto(self, histo, draw=False, *args):
        """
        Fit the histogram
        """
        self.lst = ROOT.RooArgList(self.m1, self.m2)
        self.imp = ROOT.RooFit.Import(histo)
        self.hset = ROOT.RooDataHist(
            rootID('hds_'),
            "Data set for histogram '%s'" % histo.GetTitle(),
            self.lst,
            self.imp
        )

        result = self.pdf.fitTo(self.hset,
                                ROOT.RooFit.Save(),
                                *args)

        if not draw:
            return result, None

        frame = self.m1.frame()
        self.hset.plotOn(frame)
        self.pdf .plotOn(frame,
                         ROOT.RooFit.Components(
                             self.background.pdf.GetName()),
                         ROOT.RooFit.LineStyle(ROOT.kDashed),
                         ROOT.RooFit.LineColor(ROOT.kBlue))
        self.pdf .plotOn(frame,
                         ROOT.RooFit.LineColor(ROOT.kRed))
        frame.Draw()

        return result, frame

# =============================================================================
# @class BW_pdf
#  simple wrapper over Breit-Wigner PDF
#  @see Analysis::Models::BreitWigner
#  @author Alexander BARANOV a.baranov@cern.ch
#  @date 2013-09-12


class BW_pdf(object):

    """
    Define Breit-Wigner
    """

    def __init__(self,
                 name,
                 x=None,
                 mass=None,
                 width=None,
                 m1=None,
                 m2=None,
                 L=0.0,
                 rho=0):
        self.mass = x
        self.mean = mass
        self.width = width
        self.m1 = m1
        self.m2 = m2
        self.L = L
        self.rho = rho

        self.pdf = cpp.Analysis.Models.BreitWigner(
            "bw_" + name,
            "BW(%s)" % name,
            self.mass,
            self.mean,
            self.width,
            self.m1,
            self.m2,
            self.L,
            self.rho)

        # fit
        def fitTo(self, dataset, *args, **kwargs):
            """
            Perform the fit
            """
            return self.pdf.fitTo(dataset, *args, **kwargs)


# =============================================================================
# @class PhaseSpace_pdf
#  simple wrapper over phase-space pdf
#  @see Analysis::Models::PhaseSpaceNL
#  @author Alexander BARANOV a.baranov@cern.ch
#  @date 2013-09-12
class PhaseSpace_pdf(object):

    """
    Define pdf for background : Exponential
    """
    # constructor

    def __init__(self, name, x, low, high, N, L, power=0):

        self.name = name
        self.mass = x
        self.low = low
        self.high = high
        self.N = N
        self.L = L
        self.power = power

        self.phase_space_pdf = cpp.Analysis.Models.PhaseSpaceNL(
            "PhaseSpace_" + name,
            "PhaseSpace(%s)" % name,
            self.mass,
            self.low,
            self.high,
            self.N,
            self.L)

        if power == 0:
            self.pdf = self.phase_space_pdf
        else:
            self.poly_arglist = ROOT.RooArgList()
            self.phis = [ROOT.RooRealVar('phi%d_%s' % (i, name), 'phi%d(%s)' % (i, name), 0.0, -3.5, 3.5)
                         for i in xrange(power)]

            for i in xrange(power):
                self.poly_arglist.add(self.phis[i])

            self.poly_pdf = cpp.Analysis.Models.PolyPositive(
                "poly_" + name,
                "poly%d(%s)" % (power, name),
                self.mass,
                self.poly_arglist,
                0.986,
                2.04)

            self.pdf = cpp.Analysis.Models.Product(
                "polybw_" + name,
                "polyBW(%s)" % name,
                self.poly_pdf,
                self.phase_space_pdf
            )

    # fit
    def fitTo(self, dataset, *args, **kwargs):
            """
            Perform the fit
            """
            return self.pdf.fitTo(dataset, *args, **kwargs)


# =============================================================================
def _get_Signal_(tree,
                 expr2='S1S2_sw',
                 expr1='S1S2_sw*weight'):
    """
    Extract the signal components from weighted tree 
    """
    htmp1 = ROOT.TH1D(hID(), '', 4, -1, 3)
    htmp1 . Sumw2()
    htmp2 = ROOT.TH1D(hID(), '', 4, -1, 3)
    htmp2 . Sumw2()

    tree.Draw(' 1 >> %s' % htmp1.GetName(), expr1)
    tree.Draw(' 1 >> %s' % htmp2.GetName(), expr2)

    return htmp1.accumulate(), htmp2.accumulate()

ROOT.TTree.getSignal = _get_Signal_

# =============================================================================
if '__main__' == __name__:

    print 80 * '*'
    print __doc__
    print ' Author  : ', __author__
    print ' Version : ', __version__
    print ' Date    : ', __date__
    print ' Symbols : ', __all__
    print 80 * '*'

# =============================================================================
# The END
# =============================================================================
