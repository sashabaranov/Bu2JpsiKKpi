#!/usr/bin/env python
# =============================================================================
# @file DiCharm/GoodParticles.py
#
# Helper module to collect the definitions of ``good particles''
#  for 2xCharm(&Onium) analyses
#
# @author Vanya BELYAEV Ivan.Belyaev@cern.ch
# @date   2011-06-21
#
#                   $Revision$
# Last modification $Date$
#                by $Author$
# =============================================================================
"""
Helper module to collect the definitions of ``good particles'' for
2xCharm(&Onium) analyses 
"""
# =============================================================================
__version__ = "$Revision: 124897 $"
__author__ = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__ = "2011-06-07"
# =============================================================================
__all__ = (
    #
    'prompt',  # ``promptness''
    'pidMuon',  # MuonID
    'pidKaon',  # KaonID
    'pidPion',  # PionID
    'pidProton',  # PionID
    #
    'basicDiMu',  # basic dimu     , no PID(mu) , no promptness
    'basicJpsi',  # basic J/psi    , no PID(mu) , no promptness
    'basicPsiPrime',  # basic   psi'   , no PID(mu) , no promptness
    'basicY',  # basic Upsilons , no PID(mu) , no promptness
    #
    'basicW',  # basic W
    'basicZ',  # basic Z0
    #
    'basicD0',  # basis D0        -> K- pi+           , no PID cuts
    'basicDstar',  # basis D*+ -> ( D0 -> K- pi+ ) pi+   , no PID cuts
    'basicDplus',  # basis D+        -> K- pi+ pi+       , no PID cuts
    'basicDp',  # basis D+        -> K- pi+ pi+       , no PID cuts
    'basicDs',  # basis D+/Ds+    -> (K+K-)pi+        , no PID cuts
    'basicLamC',  # basis Lambda_c+ -> pKpi             , no PID cuts
    'basicLc',  # basis Lambda_c+ -> pKpi             , no PID cuts
    #
)
# =============================================================================
from LoKiPhys.decorators import *
from LoKiCore.functions import *
from GaudiKernel.SystemOfUnits import GeV, MeV, mm
from GaudiKernel.PhysicalConstants import c_light
# =============================================================================

# =============================================================================
# Define the standard muon identification


def pidMuon():
    """
    Define the standard muon identification 
    """
    return 0 < MINTREE('mu+' == ABSID, PIDmu - PIDpi)

# =============================================================================
# Define the standard kaon identification


def pidKaon():
    """
    Define the standard kaon identification 
    """
    return 2 < MINTREE('K+' == ABSID, PIDK - PIDpi)

# =============================================================================
# Define the standard pion identification


def pidPion():
    """
    Define the standard kaon identification 
    """
    return 2 < MINTREE('pi+' == ABSID, PIDpi - PIDK)

# =============================================================================
# Define the standard proton identification


def pidProton():
    """
    Define the standard ``tight'' proton identification 
    """
    #
    fun = MINTREE('p+' == ABSID, PIDp - PIDK) > 10
    return fun & (MINTREE('p+' == ABSID, PIDp - PIDpi) > 10)

# =============================================================================
# define the proper ``prompt'' particle
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2011-06-21


def prompt():
    """
    Define proper ``prompt'' particle 
    """
    return in_range(0, DTF_CHI2NDOF(True), 5)

# ============================================================================
# define the basicDiMu
#
#   - \f$ J/\psi \rightarrow \mu^+ \mu^-\f$
#   - \f$  3 \le \mathrm{M}(\mu_1\mu_2) \le 3.2~\mathrm{GeV}/c^2 \f$
#   - \f$  \min ( p^{\mathrm{T}}(\mu_1) , p^{\mathrm{T}}(\mu_2) ) \ge 650~\mathrm{MeV}/c \f$
#   - \f$  \min ( \Delta_{KL}(\mu_1) , \Delta^{KL}(\mu_2}) \ge 5000 \f$
#   - \f$  \max ( \chi^2_{TR}(\mu_1) , \chi^2_{TR}(\mu_2)) \ge 5000 \f$
#   - \f$  \mu_1,\mu_2: \mathrm{ISMUON} \f$
#
#  @attention there is no ``prompt'' cut here!!!
#
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2011-06-21


def basicDiMu():
    """
    Define the basic  DiMuon :
    - J/psi -> mu+ mu-
    - 3 < M                < 3.2 GeV/c2
    - chi2_vx              < 20         # actually no-effect here
    - muons : INMUON & ISMUON 
    - pt(mu)               > 650 MeV
    - Kullback-Liebler(mu) > 5000
    - chi2_tr(mu)/nDof     < 5 
    """
    #
    dimu = ID == strings(['J/psi(1S)',
                          'psi(2S)',
                          'Upsilon(1S)',
                          'Upsilon(2S)',
                          'Upsilon(3S)'])
    #
    dimu = dimu & (in_range(3.0 * GeV, M,  3.2 * GeV) |
                   in_range(3.5 * GeV, M,  3.9 * GeV) |
                   in_range(8.5 * GeV, M, 12.5 * GeV))
    #
    dimu = dimu & DECTREE('Meson -> mu+ mu-')
    #
    dimu = dimu & CHILDCUT(
        1, HASMUON & ISMUON) & CHILDCUT(2, HASMUON & ISMUON)
    #
    dimu = dimu & (MINTREE('mu+' == ABSID, PT) > 650 * MeV)
    #
    dimu = dimu & (MAXTREE(ISBASIC & HASTRACK, TRCHI2DOF) < 5)
    #
    dimu = dimu & (MINTREE(ISBASIC & HASTRACK, CLONEDIST) > 5000)
    #
    dimu = dimu & (VFASPF(VCHI2) < 20)  # no effect for ``prompt'' cut
    #
    return dimu

# ============================================================================
# define the basic J/psi
#
#   - \f$ J/\psi \rightarrow \mu^+ \mu^-\f$
#   - \f$  3 \le \mathrm{M}(\mu_1\mu_2) \le 3.2~\mathrm{GeV}/c^2 \f$
#   - \f$  \min ( p^{\mathrm{T}}(\mu_1) , p^{\mathrm{T}}(\mu_2) ) \ge 650~\mathrm{MeV}/c \f$
#   - \f$  \min ( \Delta_{KL}(\mu_1) , \Delta^{KL}(\mu_2}) \ge 5000 \f$
#   - \f$  \max ( \chi^2_{TR}(\mu_1) , \chi^2_{TR}(\mu_2)) \ge 5000 \f$
#   - \f$  \mu_1,\mu_2: \mathrm{ISMUON} \f$
#
#  @attention there is no ``prompt'' cut here!!!
#
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2011-06-21


def basicJpsi():
    """
    Define the basic  J/psi :
    - J/psi -> mu+ mu-
    - 3 < M                < 3.2 GeV/c2
    - chi2_vx              < 20         # actually no-effect here
    - muons : INMUON & ISMUON 
    - pt(mu)               > 650 MeV
    - Kullback-Liebler(mu) > 5000
    - chi2_tr(mu)/nDof     < 5 
    """
    #
    psi = 'J/psi(1S)' == ID
    #
    psi = psi & in_range(3.0 * GeV, M, 3.2 * GeV)
    #
    return psi & basicDiMu()

# ============================================================================
# define the basic psi'
#
#   - \f$ J/\psi \rightarrow \mu^+ \mu^-\f$
#   - \f$  3 \le \mathrm{M}(\mu_1\mu_2) \le 3.2~\mathrm{GeV}/c^2 \f$
#   - \f$  \min ( p^{\mathrm{T}}(\mu_1) , p^{\mathrm{T}}(\mu_2) ) \ge 650~\mathrm{MeV}/c \f$
#   - \f$  \min ( \Delta_{KL}(\mu_1) , \Delta^{KL}(\mu_2}) \ge 5000 \f$
#   - \f$  \max ( \chi^2_{TR}(\mu_1) , \chi^2_{TR}(\mu_2)) \ge 5000 \f$
#   - \f$  \mu_1,\mu_2: \mathrm{ISMUON} \f$
#
#  @attention there is no ``prompt'' cut here!!!
#
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2011-06-21


def basicPsiPrime():
    """
    Define the basic psi' :
    - J/psi -> mu+ mu-
    - 3.5 < M              < 3.9 GeV/c2
    - chi2_vx              < 20         # actually no-effect here
    - muons : INMUON & ISMUON 
    - pt(mu)               > 650 MeV
    - Kullback-Liebler(mu) > 5000
    - chi2_tr(mu)/nDof     < 5 
    """
    #
    psi_p = in_range(3.5 * GeV, M, 3.9 * GeV)
    psi_p = psi_p & (ID == strings('J/psi(1S)', 'psi(2S)'))
    #
    return psi_p & basicDiMu()

# ============================================================================
# define the basic Upsilons
#
#   - \f$ J/\psi \rightarrow \mu^+ \mu^-\f$
#   - \f$  3 \le \mathrm{M}(\mu_1\mu_2) \le 3.2~\mathrm{GeV}/c^2 \f$
#   - \f$  \min ( p^{\mathrm{T}}(\mu_1) , p^{\mathrm{T}}(\mu_2) ) \ge 650~\mathrm{MeV}/c \f$
#   - \f$  \min ( \Delta_{KL}(\mu_1) , \Delta^{KL}(\mu_2}) \ge 5000 \f$
#   - \f$  \max ( \chi^2_{TR}(\mu_1) , \chi^2_{TR}(\mu_2)) \ge 5000 \f$
#   - \f$  \mu_1,\mu_2: \mathrm{ISMUON} \f$
#
#  @attention there is no ``prompt'' cut here!!!
#
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2011-06-21


def basicY():
    """
    Define the basic Upsilons :
    - J/psi -> mu+ mu-
    - 8.5 < M              < 12.5 GeV/c2
    - chi2_vx              < 20         # actually no-effect here
    - muons : INMUON & ISMUON 
    - pt(mu)               > 650 MeV
    - Kullback-Liebler(mu) > 5000
    - chi2_tr(mu)/nDof     < 5 
    """
    #
    upsilons = in_range(8.0 * GeV, M, 16 * GeV)
    upsilons = upsilons & (ID == strings(['J/psi(1S)',
                                          'psi(2S)',
                                          'Upsilon(1S)',
                                          'Upsilon(2S)',
                                          'Upsilon(3S)',
                                          'Z0']))

    upsilons = upsilons & basicDiMu()
    #
    # stripping 17
    #upsilons = upsilons & ( MINTREE ( 'mu+' == ABSID , PT ) > 1 * GeV )
    #upsilons = upsilons & ( MINTREE ( 'mu+' == ABSID , P  ) > 8 * GeV )
    #
    return upsilons

# ========================================================================
# define ``basic-charm''


def basicCharm():
    """
    """
    charm = MAXTREE(ISBASIC & HASTRACK, TRCHI2DOF) < 5
    charm = charm & (MAXTREE(ISBASIC & HASTRACK, CLONEDIST) > 5000)
    #
    charm = charm & (BPVLTIME(9) > 100 * micrometer / c_light)
    charm = charm & (MINTREE(ISBASIC & HASTRACK, MIPCHI2DV()) > 9)
    #
    return charm
# ========================================================================
# define the basic \f$ D^0 \rightarrow K^- \pi^+\f$
#  - \f$ D^0 \rightarrow K^- \pi^+ \f$
#  - tracks:
#  -- \f$ \chi^2_{TR}/nDoF <    5 \f$
#  -- \f$ \Delta^{KL}      > 5000 \f$
#  -- \f$ \chi^2_{IP}      >    9 \f$
#  -- \f$\pi,K\f$ : HASRICH
#  -  \f$ \left| \Delta \mathrm{M} \right| < 75 \mathrm{MeV}/c^2 \f$
#  =  \f$ \left| \cos \theta^* \right| < 0.9 \f$
#  -  \f$ \chi^2_{VX} < 9 \f$
#  -  \f$ \chi^2_{IP}(D^0) < 9 \f$
#  -  \f$ c\times\tau(D^0) > 100 \mu{\mathrm{m}} \f$
#
#  @attention there is no PID cuts here
#
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2011-06-21


def basicD0():
    """
    Define the basic cuts for D0 -> K- pi+ 
    """
    #
    d0 = 'D0' == ABSID
    #
    ## d0  =  d0 & ( PT > 2 * GeV )
    d0 = d0 & (ADMASS('D0') < 75 * MeV)
    d0 = d0 & (VFASPF(VCHI2) < 9)
    d0 = d0 & DECTREE('[D0 -> K- pi+]CC')
    #
    d0 = d0 & (abs(LV01) < 0.9)
    #
    # added to avoid pion kaon PID
    d0 = d0 & (MINTREE(ISBASIC & HASTRACK, P) > 3.2 * GeV)
    d0 = d0 & (MAXTREE(ISBASIC & HASTRACK, P) < 100 * GeV)
    d0 = d0 & (MINTREE(ISBASIC & HASTRACK, ETA) > 2)
    d0 = d0 & (MAXTREE(ISBASIC & HASTRACK, ETA) < 5.0)
    #
    d0 = d0 & CHILDCUT(1, HASRICH) & CHILDCUT(2, HASRICH)
    #
    return d0 & basicCharm()


# ========================================================================
# define the basic \f$ D^{*+} \rightarrow D^0 \pi^+ \f$
#  - \f$ D_s^+ \rightarrow ( K^- K+ ) \pi^+\f$
#  - tracks:
#  -- \f$ \chi^2_{TR}/nDoF <    5 \f$
#  -- \f$ \Delta^{KL}      > 5000 \f$
#  -- \f$ \chi^2_{IP}      >    9 \f$
#  -- \f$\pi,K\f$ : HASRICH
#  -  \f$ \left| \Delta \mathrm{M} \right| < 75 \mathrm{MeV}/c^2 \f$
#  -  \f$ \chi^2_{VX} < 9 \f$
#  -  \f$ \chi^2_{IP}(D^+) < 9 \f$
#  -  \f$ c\times\tau(D^+) > 100 \mu{\mathrm{m}} \f$
#
#  @attention there is no PID cuts here
#
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2011-06-21
def basicDstar():
    """
    Define the basic cuts for D0 -> K- pi+ pi+ 
    
    """
    dstar = 'D*(2010)+' == ABSID
    dstar = dstar & DECTREE('[D*(2010)+ -> ( D0 -> K- pi+)  pi+]CC')
    dstar = dstar & (M - M1 < 155 * MeV)
    dstar = dstar & (VFASPF(VCHI2) < 36)
    dstar = dstar & CHILDCUT(1, basicD0())
    #
    return dstar

# ========================================================================
# define the basic \f$ D^+\rightarrow K^- \pi^+ \pi^+ \f$
#  - \f$ D^+ \rightarrow K^- \pi^+ \pi^+\f$
#  - tracks:
#  -- \f$ \chi^2_{TR}/nDoF <    5 \f$
#  -- \f$ \Delta^{KL}      > 5000 \f$
#  -- \f$ \chi^2_{IP}      >    9 \f$
#  -- \f$\pi,K\f$ : HASRICH
#  -  \f$ \left| \Delta \mathrm{M} \right| < 75 \mathrm{MeV}/c^2 \f$
#  -  \f$ \chi^2_{VX} < 25 \f$
#  -  \f$ \chi^2_{IP}(D^+) < 9 \f$
#  -  \f$ c\times\tau(D^+) > 100 \mu{\mathrm{m}} \f$
#
#  @attention there is no PID cuts here
#
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2011-06-21


def basicDplus():
    """
    Define the basic cuts for D0 -> K- pi+ pi+ 
    
    """
    dp = 'D+' == ABSID
    ## dp  =  dp & ( PT > 2 * GeV )
    dp = dp & (ADMASS('D+') < 50 * MeV)
    dp = dp & (VFASPF(VCHI2) < 25)
    dp = dp & DECTREE('[D+ -> K- pi+ pi+]CC')
    #
    dp = dp & (MINTREE(ISBASIC & HASTRACK, P) > 3.2 * GeV)
    dp = dp & (MAXTREE(ISBASIC & HASTRACK, P) < 100 * GeV)
    dp = dp & (MINTREE(ISBASIC & HASTRACK, ETA) > 2)
    dp = dp & (MAXTREE(ISBASIC & HASTRACK, ETA) < 5.0)
    #
    dp = dp & CHILDCUT(1, HASRICH) & CHILDCUT(
        2, HASRICH) & CHILDCUT(3, HASRICH)
    #
    return dp & basicCharm()

basicDp = basicDplus
# ========================================================================
# define the basic \f$ D_s^+\rightarrow ( K^- K^+)  \pi^+ \f$
#  - \f$ D_s^+ \rightarrow ( K^- K+ ) \pi^+\f$
#  - tracks:
#  -- \f$ \chi^2_{TR}/nDoF <    5 \f$
#  -- \f$ \Delta^{KL}      > 5000 \f$
#  -- \f$ \chi^2_{IP}      >    9 \f$
#  -- \f$\pi,K\f$ : HASRICH
#  -  \f$ \left| \Delta \mathrm{M} \right| < 75 \mathrm{MeV}/c^2 \f$
#  -  \f$ \chi^2_{VX} < 25 \f$
#  -  \f$ \chi^2_{IP}(D^+) < 9 \f$
#  -  \f$ c\times\tau(D^+) > 100 \mu{\mathrm{m}} \f$
#
#  @attention there is no PID cuts here
#
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2011-06-21


def basicDs():
    """
    Define the basic cuts for D0 -> K- pi+ pi+ 
    
    """
    ds = 'D_s+' == ABSID

    ## ds  =  ds & ( PT > 2 * GeV )
    ds = ds & ((ADMASS('D_s+') < 75 * MeV) |
               (ADMASS('D+') < 75 * MeV))
    ds = ds & (VFASPF(VCHI2) < 25)
    ds = ds & DECTREE('[D_s+ -> K- K+ pi+]CC')
    #
    ds = ds & (M12 < 1040 * MeV)
    #
    ds = ds & (MINTREE(ISBASIC & HASTRACK, P) > 3.2 * GeV)
    ds = ds & (MAXTREE(ISBASIC & HASTRACK, P) < 100 * GeV)
    ds = ds & (MINTREE(ISBASIC & HASTRACK, ETA) > 2)
    ds = ds & (MAXTREE(ISBASIC & HASTRACK, ETA) < 5.0)
    #
    ds = ds & CHILDCUT(1, HASRICH) & CHILDCUT(
        2, HASRICH) & CHILDCUT(3, HASRICH)
    #
    return ds & basicCharm()

# ========================================================================
# define the basic \f$ \Lambda_c^+\rightarrow p K^-\pi^+ \f$
#  - \f$ \Lambda_c^+ \rightarrow p K^-  \pi^+ \f$
#  - tracks:
#  -- \f$ \chi^2_{TR}/nDoF <    5 \f$
#  -- \f$ \Delta^{KL}      > 5000 \f$
#  -- \f$ \chi^2_{IP}      >    9 \f$
#  -- \f$p,\pi,K\f$ : HASRICH
#  -  \f$ \left| \Delta \mathrm{M} \right| < 75 \mathrm{MeV}/c^2 \f$
#  -  \f$ \chi^2_{VX} < 25 \f$
#  -  \f$ \chi^2_{IP}(\Lambda_c^+) < 9 \f$
#  -  \f$ c\times\tau(\Lambda_c^+) > 100 \mu{\mathrm{m}} \f$
#
#  @attention there is no PID cuts here
#
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2011-06-21


def basicLamC():
    """
    Define the basic cuts for Lamda_c+ -> p K- pi+ 
    
    """
    lc = 'Lambda_c+' == ABSID

    ## lc  =  lc & ( PT > 2 * GeV )
    lc = lc & (ADMASS('Lambda_c+') < 75 * MeV)
    lc = lc & (VFASPF(VCHI2) < 25)
    lc = lc & DECTREE('[Lambda_c+ -> p+ K- pi+]CC')
    #
    #
    #
    # specific for Lambda_c+
    #
    lc = lc & (MINTREE(ISBASIC & HASTRACK, PT) > 0.5 * GeV)
    lc = lc & (BPVLTIME(9) < (500 * micrometer / c_light))
    lc = lc & (CHILD(P, 'p+' == ABSID) > 10 * GeV)
    #
    lc = lc & (MINTREE(ISBASIC & HASTRACK, P) > 3.2 * GeV)
    lc = lc & (MAXTREE(ISBASIC & HASTRACK, P) < 100 * GeV)
    lc = lc & (MINTREE(ISBASIC & HASTRACK, ETA) > 2)
    lc = lc & (MAXTREE(ISBASIC & HASTRACK, ETA) < 5.0)
    #
    lc = lc & CHILDCUT(1, HASRICH) & CHILDCUT(
        2, HASRICH) & CHILDCUT(3, HASRICH)
    #
    return lc & basicCharm()


basicLc = basicLamC


# ========================================================================
# define the basic \f$ W^{\pm}\f$
#
#  @attention there is no PID cuts here
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-21
def basicW():
    """
    Define the basic cuts for W+ 
    
    """
    #
    w = 'mu+' == ABSID
    #
    w = w & (PT > 10 * GeV)
    w = w & in_range(2, ETA, 4.5)
    w = w & (PERR2 < (0.1 ** 2 * P ** 2))
    #
    return w

# ========================================================================
# define the basic \f$ Z^{0}\f$
#
#  @attention there is no PID cuts here
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-21


def basicZ():
    """
    Define the basic cuts for Z0 
    
    """
    #
    z0 = in_range(20 * GeV, M, 160 * GeV)
    z0 = z0 & (ID == strings(['Z0',
                              'J/psi(1S)',
                              'psi(2S)',
                              'Upsilon(1S)',
                              'Upsilon(2S)',
                              'Upsilon(3S)']))
    z0 = z0 & DECTREE(' X0 -> mu+ mu-')
    #
    z0 = z0 & (MINTREE('mu+' == ABSID, PT) > 5 * GeV)
    #
    return z0


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
