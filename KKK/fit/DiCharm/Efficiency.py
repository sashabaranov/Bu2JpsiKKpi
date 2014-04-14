#!/usr/bin/env python
# =============================================================================
# @file DiCharm/Efficiency.py
#
# Helper module to collect various efficiencies for 2xCharm(&Onium) studies
#
# @author Vanya BELYAEV Ivan.Belyaev@cern.ch
# @date   2011-06-19
#
#                   $Revision$
# Last modification $Date$
#                by $Author$
# =============================================================================
"""
Helper module to collect various efficiencies for 2xCharm(&onium) studies
"""
# =============================================================================
__version__ = "$Revision: 124897 $"
__author__ = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__ = "2011-06-07"
# =============================================================================
__all__ = (
    # =========================================================================
    # Efficiency of Generator Cuts
    # =========================================================================
    'genEff_D0',  # generator level cuts for D0   evttype 27163003
    'genEff_Dstar',  # generator level cuts for D*+  evttype 27163003
    'genEff_Dplus',  # generator level cuts for D+   evttype 21263010
    'genEff_Dp',  # generator level cuts for D+   evttype 21263010
    'genEff_Ds',  # generator level cuts for Ds+  evttype 23263001
    'genEff_Jpsi',  # generator level cuts for Jpsi evttype 24142001
    'genEff_Lc',  # generator level cuts for Lc+  evttype 25103000
    # =========================================================================
    # Reconstruction&Selection Efficiencies
    # =========================================================================
    'selEff_D0',  # reconstruction & selection efficiency for D0
    'selEff_Dp',  # reconstruction & selection efficiency for D+
    'selEff_Ds',  # reconstruction & selection efficiency for Ds+
    'selEff_Lc',  # reconstruction & selection efficiency for Lambda_c+
    'selEff_Jpsi_3D',  # reconstruction & selection efficiency for J/psi
    'selEff_Jpsi_2D',  # reconstruction & selection efficiency for J/psi
    # =========================================================================
    # Trigger (TOS&TRG) efficiencies
    # =========================================================================
    'tosEff_Jpsi',  # TOS-efficiency for J/psi
    'tosEff_D0',  # TOS-efficiency for D0
    'tosEff_Dp',  # TOS-efficiency for D+
    'tosEff_Ds',  # TOS-efficiency for Ds+
    'tosEff_Lc',  # TOS-efficiency for Lc+
    # =========================================================================
    # Particle identification efficiencies
    #  @thanks Andrew Powell
    # =========================================================================
    'pidEff_Jpsi',  # (di)muon idenitfication efficiency ( per J/psi!)
    'pidEff_pion',  # pion     idenitfication efficiency
    'pidEff_kaon',  # kaon     idenitfication efficiency
    'pidEff_proton',  # proton   idenitfication efficiency
    'pidEff_Jpsi_corr',  # correction factor DATA/MC for ISMUON criteria
    # =========================================================================
    # Correction factor to track recontruction efficiency
    #  @thanks Jeroen van Tilburg
    # =========================================================================
    'trackEff_corr',  # Correction factor to track recontruction efficiency
    # =========================================================================
    # branching  ratios, PDG'2k+10
    # =========================================================================
    'br_Jpsi',  # branching  ratio: J/psi     -> mu+ mu-
    'br_D0',  # branching  ratio: D0        -> K-  pi+
    'br_Dp',  # branching  ratio: D+        -> K-  pi+ pi+
    'br_Ds',  # branching  ratio: Ds+       -> ( phi -> K- K+ ) pi+
    'br_Lc',  # branching  ratio: Lambda_c+ -> p+ K- pi+
    'br_ratio'  # dictionry of branching ratios
    # =========================================================================
    # DCSD stuff
    # =========================================================================
    'r_DCSD',  # DCSD ratio for D0 decays
    'corr_DCSD',  # function for DCSD corrections
    # =========================================================================
    # Luminosity
    # =========================================================================
    'lumi_DB',  # for mDB
    'lumi_EXT',  # extrapolated using prompt charm peaks
    # =========================================================================
    # Global Event Cuts
    'gecEff_Spd_Jpsi',  # efficiency of GEC/#Spd cut for J/psiD events
    'gecEff_Spd_CC',  # efficiency of GEC/#Spd cut for CC     events
    'gecEff_Spd_CCbar',  # efficiency of GEC/#Spd cut for CC-bar events
    'gecEff_OT_Jpsi',  # efficiency of GEC/#OT  cut for J/psiD events
    'gecEff_OT_CC',  # efficiency of GEC/#OT  cut for CC     events
    'gecEff_OT_CCbar',  # efficiency of GEC/#OT  cut for CC-bar events
    # =========================================================================
)
# =============================================================================
from AnalysisPython.PyRoUts import VE
import AnalysisPython.ZipShelve as ZipShelve
import os
# =============================================================================
# Lumi
# =============================================================================
# lumi_DOWN  = VE(199, 7**2) ## form DB
lumi_DB = VE(341, 12 ** 2)  # from DB
lumi_EXT = VE(355, 13 ** 2)  # extrapolated using prompt charm peaks
# =============================================================================
# PDG'2k+10
# =============================================================================
br_Jpsi = VE(5.93, 0.06 ** 2) / 100.0
br_D0 = VE(3.89, 0.05 ** 2) / 100.0
br_Dp = VE(9.4, 0.4 ** 2) / 100.0
br_Ds = VE(2.32, 0.14 ** 2) / 100.0
br_Lc = VE(5.0, 1.3 ** 2) / 100.0
# =============================================================================
br_ratio = {}
br_ratio['J/psi(1S)'] = br_Jpsi
br_ratio['Jpsi'] = br_Jpsi
br_ratio['JPsi'] = br_Jpsi
br_ratio['Psi'] = br_Jpsi
br_ratio['psi'] = br_Jpsi

br_ratio['D0'] = br_D0
br_ratio['d0'] = br_D0

br_ratio['D+'] = br_Dp
br_ratio['Dp'] = br_Dp
br_ratio['d+'] = br_Dp
br_ratio['dp'] = br_Dp

br_ratio['D_s+'] = br_Ds
br_ratio['Ds+'] = br_Ds
br_ratio['Ds'] = br_Ds
br_ratio['d_s+'] = br_Ds
br_ratio['ds+'] = br_Ds
br_ratio['ds'] = br_Ds

br_ratio['Lambda_c+'] = br_Lc
br_ratio['Lambda_c'] = br_Lc
br_ratio['lambda_c+'] = br_Lc
br_ratio['lambda_c'] = br_Lc
br_ratio['Lam_c+'] = br_Lc
br_ratio['Lam_c'] = br_Lc
br_ratio['L_c+'] = br_Lc
br_ratio['L_c'] = br_Lc
br_ratio['Lc+'] = br_Lc
br_ratio['Lc'] = br_Lc

# =============================================================================
# Double Cabibbo-suppressed ratios
# =============================================================================
dcsd_D0 = VE(3.80, 0.18 ** 2) * 1.e-3  # D0 -> K+pi- /D->K-pi+
dcsd_fact = (1 + dcsd_D0 ** 2) ** -0.5
r_DCSD = dcsd_D0
f_DCSD = (1 + r_DCSD ** 2) ** -0.5
r2_DCSD = 2 * dcsd_D0
f2_DCSD = (1 + r2_DCSD ** 2) ** -0.5
# correction for DCSD decays


def corr_DCSD(gamma_0, gamma_1, DCSD):
    """
    Correction for DCSD Decays

    >>> gc0, gc1 = corr_DCSD ( g0 , g1 , r_DSCD )
    
    """
    #
    corr_0 = gamma_0 - gamma_1 * DCSD
    corr_1 = gamma_1 - gamma_0 * DCSD

    fact = (1 - DCSD ** 2) ** -0.5

    corr_0 *= fact
    corr_1 *= fact

    return corr_0, corr_1
# =============================================================================
# Global Event Cuts
# =============================================================================
# efficiency of GEC/#Spd cut for J/psiD events
gecEff_Spd_Jpsi = VE(96.6, 0.5 ** 2) / 100
# efficiency of GEC/#Spd cut for CC     events
gecEff_Spd_CC = VE(91.3, 0.9 ** 2) / 100
# efficiency of GEC/#Spd cut for CC-bar events
gecEff_Spd_CCbar = VE(92.6, 0.9 ** 2) / 100
# efficiency of GEC/#OT  cut for J/psiD events
gecEff_OT_Jpsi = VE(97.1, 0.4 ** 2) / 100
# efficiency of GEC/#OT  cut for CC     events
gecEff_OT_CC = VE(99.3, 0.3 ** 2) / 100
# efficiency of GEC/#OT  cut for CC-bar events
gecEff_OT_CCbar = VE(98.8, 0.2 ** 2) / 100

# =============================================================================
# correction factor DATA/MC for ISMUON criteria : 1.024+-0.011
# =============================================================================
pidEff_Jpsi_corr = VE(1.024, 0.011 ** 2)

# =============================================================================
# get the object from read-only db


def _get_eff_from_db(key, dbase='$DICHARMROOT/data/DiCharm.db'):
    """
    Get object from read-only DB 
    """
    dbase = os.path.expandvars(dbase)
    dbase = os.path.expanduser(dbase)

    if not os.path.exists(dbase):
        db_ = './DiCharm.db'
        if os.path.exists(db_):
            dbase = './DiCharm.db'
            print " Try to use Local-DB:'%s'" % dbase

    db = ZipShelve.open(dbase, 'r')
    obj = db[key]
    db.close()
    return obj

# ========================================================================

# ========================================================================
# Generator acceptance
# ========================================================================
_gen_eff_dstar = _get_eff_from_db('GenCuts: D*+, 27163003')
_gen_eff_d0 = _get_eff_from_db('GenCuts: D0, 27163003')
_gen_eff_dplus = _get_eff_from_db('GenCuts: D+, 21263010')
_gen_eff_ds = _get_eff_from_db('GenCuts: Ds+, 23263001')
_gen_eff_lc = _get_eff_from_db('GenCuts: Lc+, 25103000')
## _gen_eff_psi_1 = _get_eff_from_db ('GenCuts: J/psi, 24142001-1')
## _gen_eff_psi_2 = _get_eff_from_db ('GenCuts: J/psi, 24142001-2')
_gen_eff_psi = _get_eff_from_db('GenCuts: J/psi, 24142001')

# ========================================================================
# get the generator cuts efficiency for D*+
#  @param pt  pt ( in GeV)
#  @param y   rapidity


def genEff_Dstar(pt, y):
    """
    Get the generator cut efficiency for D*+

    >>> eff = genEff_Dstar ( pt , y )
    
    """
    return _gen_eff_dstar(pt, y)
# ========================================================================
# get the generator cuts efficiency for D+
#  @param pt  pt ( in GeV)
#  @param y   rapidity


def genEff_Dplus(pt, y):
    """
    Get the generator cut efficiency for D+

    >>> eff = genEff_Dplus ( pt , y )
    
    """
    return _gen_eff_dplus(pt, y)

# ditto
genEff_Dp = genEff_Dplus

# ========================================================================
# get the generator cuts efficiency for D0
#  @param pt  pt ( in GeV)
#  @param y   rapidity


def genEff_D0(pt, y):
    """
    Get the generator cut efficiency for D0

    >>> eff = genEff_D0 ( pt , y )
    
    """
    return _gen_eff_d0(pt, y)

# ========================================================================
# get the generator cuts efficiency for Ds+
#  @param pt  pt ( in GeV)
#  @param y   rapidity


def genEff_Ds(pt, y):
    """
    Get the generator cut efficiency for Ds+

    >>> eff = genEff_Ds ( pt , y )
    
    """
    return _gen_eff_ds(pt, y)

# ========================================================================
# get the generator cuts efficiency for Lc+
#  @param pt  pt ( in GeV)
#  @param y   rapidity


def genEff_Lc(pt, y):
    """
    Get the generator cut efficiency for Lambda_c+

    >>> eff = genEff_Lc ( pt , y )
    
    """
    return _gen_eff_lc(pt, y)

# ========================================================================
# get the generator cuts efficiency for J/psi
#  @param pt  pt ( in GeV)
#  @param y   rapidity


def genEff_Jpsi(pt, y):
    """
    Get the generator cut efficiency for J/psi

    >>> eff = genEff_Jpsi ( pt , y )
    
    """
    return _gen_eff_psi(pt, y)


# ========================================================================
## Recontruction & Selection
# ========================================================================
_sel_eff_D0 = _get_eff_from_db('EffCuts: D0->K-pi+, Reco&Sel')
_sel_eff_Dp = _get_eff_from_db('EffCuts: D+->K-pi+pi+, Reco&Sel')
_sel_eff_Ds = _get_eff_from_db('EffCuts: Ds+->(K-K+)pi+, Reco&Sel')
_sel_eff_Lc = _get_eff_from_db('EffCuts: Lambda_c+ ->pK-pi+, Reco&Sel')
_sel_eff_psi_2D = _get_eff_from_db('EffCuts: Jpsi->mu+mu-, Reco&Sel-2D')
_sel_eff_psi_3D = _get_eff_from_db('EffCuts: Jpsi->mu+mu-, Reco&Sel-3D')


# ========================================================================
# get the Reconstruction&Selection efficiency for D0
#  @param pt  pt ( in GeV)
#  @param y   rapidity
def selEff_D0(pt, y):
    """
    Get the reconsctruction & selecton efficiency for D0

    >>> eff = sellEff_D0 ( pt , y )
    
    """
    return _sel_eff_D0(pt, y)
# ========================================================================
# get the Reconstruction&Selection efficiency for D+
#  @param pt  pt ( in GeV)
#  @param y   rapidity


def selEff_Dp(pt, y):
    """
    Get the reconsctruction & selecton efficiency for D+

    >>> eff = sellEff_Dp ( pt , y )
    
    """
    return _sel_eff_Dp(pt, y)

# ========================================================================
# get the Reconstruction&Selection efficiency for Ds+
#  @param pt  pt ( in GeV)
#  @param y   rapidity


def selEff_Ds(pt, y):
    """
    Get the reconsctruction & selecton efficiency for Ds+

    >>> eff = sellEff_Ds ( pt , y )
    
    """
    return _sel_eff_Ds(pt, y)

# ========================================================================
# get the Reconstruction&Selection efficiency for Lambda_c+
#  @param pt  pt ( in GeV)
#  @param y   rapidity


def selEff_Lc(pt, y):
    """
    Get the reconsctruction & selecton efficiency for Lambda_c+

    >>> eff = sellEff_Lc ( pt , y )
    
    """
    return _sel_eff_Lc(pt, y)

# ========================================================================
# get the Reconstruction&Selectrion efficiency for J/psi
#  @param pt  pt ( in GeV)
#  @param y   rapidity


def selEff_Jpsi_2D(pt, y):
    """
    Get the reconsctruction & selecton efficiency for J/psi

    >>> eff = selEff_Jpsi_2D ( pt , y )
    
    """
    return _sel_eff_psi_2D(pt, y)


# ========================================================================
# get the Reconstruction&Selectrion efficiency for J/psi
#  @param pt     pt ( in GeV)
#  @param y      rapidity
#  @param lv01   cosine of decay angle
def selEff_Jpsi_3D(pt, y, lv01):
    """
    Get the reconsctruction & selecton efficiency for J/psi
    
    >>> eff = selEff_Jpsi_3D ( pt , y , lv01 )
    
    """
    #
    if lv01 < 0:
        lv01 = abs(lv01)
    #
    ptAxis, yAxis, table = _sel_eff_psi_3D
    #
    iPt = ptAxis.FindBin(pt)
    if not 1 <= iPt <= ptAxis.GetNbins():
        return VE(0, 0)
    #
    iY = yAxis.FindBin(y)
    if not 1 <= iY <= yAxis.GetNbins():
        return VE(0, 0)
    #
    eff = table[iPt][iY]
    #
    v = eff(lv01)
    #

    ### alpha  =  1.0
    ### c      =  3.0/(3+alpha)*(1+alpha*lv01*lv01)
    ### v     /=  c
    #
    return v

# =============================================================================
# trigger efficiencies
# =============================================================================
_trg_eff_tos_ = {}
_trg_eff_trg_ = {}

for p in ('D0', 'Dp', 'Ds', 'Lc', 'J/psi'):

    _trg_eff_tos_[p] = {}
    _trg_eff_trg_[p] = {}

    for t in ('L0', 'L1', 'L2', 'L0xL1', 'L0xL1xL2'):

        key_tos = '%s,TisTos,%s, eTOS,Function' % (p, t)
        key_trg = '%s,TisTos,%s, eTRG,Function' % (p, t)

        _trg_eff_tos_[p][t] = _get_eff_from_db(key_tos)[0]
        _trg_eff_trg_[p][t] = _get_eff_from_db(key_trg)[0]


# ========================================================================
# TOS-efficiency for particle
#  @param pt   transverse momentum ( in GeV)
#  @param y    rapidity
#  @param part the particle
#  @param trig the trigger
def _tosEff(pt, y, part, trig, tos):
    """
    Get TOS-efficienciens 
    """

    if tos:
        return _trg_eff_tos_[part][trig](pt, y)
    else:
        return _trg_eff_trg_[part][trig](pt, y)

# ===============================================================
# TOS-efficiency for J/psi
#  @param pt   transverse momentum ( in GeV)
#  @param y    rapidity
#  @param trig the trigger


def tosEff_Jpsi(pt, y, trig, tos=True):
    """
    TOS-efficiency for J/psi

    >>> eff_L0       = tosEff_Jpsi ( 5 , 2.5 , 'L0'       )
    >>> eff_L0xL1xL2 = tosEff_Jpsi ( 5 , 2.5 , 'L0xL1xL2' )    
    """
    return _tosEff(pt, y, 'J/psi', trig, tos)
# ===============================================================
# TOS-efficiency for D0
#  @param pt   transverse momentum ( in GeV)
#  @param y    rapidity
#  @param trig the trigger


def tosEff_D0(pt, y, trig, tos=True):
    """
    TOS-efficiency for D0

    >>> eff_L0       = tosEff_D0 ( 5 , 2.5 , 'L0'       )
    >>> eff_L0xL1xL2 = tosEff_D0 ( 5 , 2.5 , 'L0xL1xL2' )    
    """
    return _tosEff(pt, y, 'D0', trig, tos)
# ===============================================================
# TOS-efficiency for D+
#  @param pt   transverse momentum ( in GeV)
#  @param y    rapidity
#  @param trig the trigger


def tosEff_Dp(pt, y, trig, tos=True):
    """
    TOS-efficiency for D+

    >>> eff_L0       = tosEff_Dp ( 5 , 2.5 , 'L0'       )
    >>> eff_L0xL1xL2 = tosEff_Dp ( 5 , 2.5 , 'L0xL1xL2' )    
    """
    return _tosEff(pt, y, 'Dp', trig, tos)
# ===============================================================
# TOS-efficiency for Ds+
#  @param pt   transverse momentum ( in GeV)
#  @param y    rapidity
#  @param trig the trigger


def tosEff_Ds(pt, y, trig, tos=True):
    """
    TOS-efficiency for Ds+

    >>> eff_L0       = tosEff_Ds ( 5 , 2.5 , 'L0'       )
    >>> eff_L0xL1xL2 = tosEff_Ds ( 5 , 2.5 , 'L0xL1xL2' )    
    """
    return _tosEff(pt, y, 'Ds', trig, tos)
# ===============================================================
# TOS-efficiency for Lambda_c+
#  @param pt   transverse momentum ( in GeV)
#  @param y    rapidity
#  @param trig the trigger


def tosEff_Lc(pt, y, trig, tos=True):
    """
    TOS-efficiency for Lambda_c+

    >>> eff_L0       = tosEff_Lc ( 5 , 2.5 , 'L0'       )
    >>> eff_L0xL1xL2 = tosEff_Lc ( 5 , 2.5 , 'L0xL1xL2' )    
    """
    return _tosEff(pt, y, 'Lc', trig, tos)


# ========================================================================
# PID efficiencies
# ========================================================================
_pid_eff_psi = _get_eff_from_db('J/psi,MuPID, Function')[0]
_pid_eff_pion = _get_eff_from_db('PionPID, Histos')
_pid_eff_kaon = _get_eff_from_db('KaonPID, Histos')
_pid_eff_proton = _get_eff_from_db('ProtonPID, Histos')


# =============================================================================
# (di)muon identification efficiency for J/psi
#  @param   pt   transverse momentum of J/psi ( in GeV
#  @param   y    rapidity of J/psi
#  @return  the efficiency for di-muon identification
def pidEff_Jpsi(pt, y):
    """
    (Di) Muon identificantion efficiency
    
    >>> psi_pt = ...
    >>> psi_y = ... 
    >>> eff_Jpsi_muPID   = pidEff_Jpsi ( psi_pt , psi_y )
    
    """
    return _pid_eff_psi(pt, y) * pidEff_Jpsi_corr
    # return VE ( _pid_eff_psi ( pt , y ).value() , 0 )  * pidEff_Jpsi_corr
    # return _pid_eff_psi ( pt , y ) * VE ( pidEff_Jpsi_corr.value() , 0 )

# =============================================================================
# hadron   identification efficiency


def _pidEffHadron(p, eta, ntrk, data):
    #
    ntrk_axis = data[0]
    pid = data[1]
    #
    ntrk_bin = ntrk_axis.FindBin(ntrk)
    #
    min_bin = 1
    max_bin = ntrk_axis.GetNbins()
    #
    if min_bin > ntrk_bin:
        return pid[min_bin][3](p * 1000.0, eta)
    if max_bin < ntrk_bin:
        return pid[max_bin][3](p * 1000.0, eta)
    #
    bin_center = ntrk_axis.GetBinCenter(ntrk_bin)
    #
    if min_bin == ntrk_bin and ntrk <= bin_center:
        return pid[min_bin][3](p * 1000.0, eta)
    elif max_bin == ntrk_bin and ntrk >= bin_center:
        return pid[max_bin][3](p * 1000.0, eta)
    #
    bin1 = ntrk_bin + 1 if ntrk > bin_center else ntrk_bin - 1
    #
    eff = pid[ntrk_bin][3](p * 1000.0, eta)
    eff1 = pid[bin1][3](p * 1000.0, eta)
    #
    bc1 = ntrk_axis.GetBinCenter(bin1)
    bc0 = bin_center
    #
    c1 = (ntrk - bc0) / (bc1 - bc0)
    c0 = (ntrk - bc1) / (bc0 - bc1)
    #
    return eff * c0 + eff1 * c1

# =============================================================================
# Pion identification efficiency
#  @param p    pion momentum (in Gev/c)
#  @param eta  pion pseudorapidity
#  @param ntrk numbef of tracks
#  @return pion identification efficiency
#  @thanks Andrew Powell


def pidEff_pion(p, eta, ntrk):
    """
    Pion identificantion efficiency

    >>> p    = ...
    >>> eta  = ...
    >>> ntrk = ... 
    >>> eff_pionPID = pidEff_pion ( p , eta , ntrk ) 
    
    """
    return _pidEffHadron(p, eta, ntrk, _pid_eff_pion)

# =============================================================================
# Kaon identification efficiency
#  @param p    kaon momentum (in Gev/c)
#  @param eta  kaon pseudorapidity
#  @param ntrk numbef of tracks
#  @return kaon identification efficiency
#  @thanks Andrew Powell


def pidEff_kaon(p, eta, ntrk):
    """
    Kaon identificantion efficiency

    >>> p    = ...
    >>> eta  = ...
    >>> ntrk = ... 
    >>> eff_kaonPID = pidEff_kaon ( p , eta , ntrk ) 
    """
    return _pidEffHadron(p, eta, ntrk, _pid_eff_kaon)


# =============================================================================
# Proton identification efficiency
#  @param p    proton momentum (in GeV/c)
#  @param eta  proton pseudorapidity
#  @param ntrk numbef of tracks
#  @return proton identification efficiency
#  @thanks Andrew Powell
def pidEff_proton(p, eta, ntrk):
    """
    Proton identificantion efficiency

    >>> p    = ...
    >>> eta  = ...
    >>> ntrk = ... 
    >>> eff_protonPID = pidEff_proton ( p , eta , ntrk  )
    """
    return _pidEffHadron(p, eta, ntrk, _pid_eff_proton)


# ========================================================================
# Track reconstruction
# ========================================================================
_trk_eff_corr = _get_eff_from_db('TrackEff: Data/MC-ratio')
_trk_eff_p_min_ = _trk_eff_corr.GetXaxis().GetXmin() + 0.001
_trk_eff_p_max_ = _trk_eff_corr.GetXaxis().GetXmax() - 0.001
_trk_eff_eta_min_ = _trk_eff_corr.GetYaxis().GetXmin() + 0.001
_trk_eff_eta_max_ = _trk_eff_corr.GetYaxis().GetXmax() - 0.001


# ============================================================================
# DATA/MC correction for the track-reconstruction efficiency
#  @param p    track momentum (in GeV/c)
#  @param eta  track pseudorapidity
#  @param ntrk track multiplicity
#  @return correction factor for the track reconstruction efficiency
#  @thanks Jeroen van Tilburg
def trackEff_corr(p, eta, ntrk=None):
    """
    DATA/MC correction factor for the track recontruction efficiency

    >>> 
    
    >>> p    = ...
    >>> eta  = ...
    >>> corr_factor = trackEff_corr ( p , eta )
    """
    p_ = max(min(p, _trk_eff_p_max_), _trk_eff_p_min_)
    eta_ = max(min(eta, _trk_eff_eta_max_), _trk_eff_eta_min_)

    v = _trk_eff_corr(p_, eta_)
    if 0 == v.value():
        print ' zero ???', p, eta, p_, v
        return VE(1, 0)

    return v

# ========================================================================
# Z0 reconstruction
# ========================================================================


def _readEfficiency(eta, charge, name):
    if 2 > eta or 4.5 < eta:
        return VE(0, 0)  # efficiencies are not valid there.
    charge = {1: 'muplus',
              -1: "muminus"}[charge]
    h = _get_eff_from_db(name)[charge]
    return VE(h.GetBinContent(h.FindBin(eta)),
              h.GetBinError(h.FindBin(eta)))


pidEff_muon = lambda eta, charge: _readEfficiency(eta, charge, "MuonPID")
trackEff_muon = lambda eta, charge: _readEfficiency(
    eta, charge, "MuonTracking")
triggerEff_muon = lambda eta, charge: _readEfficiency(
    eta, charge, "MuonTrigger")


def eff_GEC(nPV):
    h = _get_eff_from_db("GEC,nPV")
    b = h.FindBin(nPV)
    if b == 0:
        b = 1
    if b > h.GetNbinsX():
        b = h.GetNbinsX()
    return VE(h.GetBinContent(b), h.GetBinError(b))


EffOR = lambda e1, e2: e1 + e2 - e1 * e2
EffAND = lambda e1, e2: e1 * e2

pidEff_Z = lambda etapos, etaneg: EffAND(pidEff_muon(etapos, +1),
                                         pidEff_muon(etaneg, -1))
trackEff_Z = lambda etapos, etaneg: EffAND(trackEff_muon(etapos, +1),
                                           pidEff_muon(etaneg, -1))
triggerEff_Z = lambda etapos, etaneg: EffOR(trackEff_muon(etapos, +1),
                                            pidEff_muon(etaneg, -1))

eff_Z = lambda etapos, etaneg, nPV: reduce(EffAND,
                                           [pidEff_Z(etapos, etaneg),
                                            trackEff_Z(etapos, etaneg),
                                            triggerEff_Z(etapos, etaneg),
                                            eff_GEC(nPV)])


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
