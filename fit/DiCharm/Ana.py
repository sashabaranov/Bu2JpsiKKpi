#!/usr/bin/env python
# =============================================================================
# @file DiCharm/Ana.py
#
# Set of helper functions to calculate efficiency base on n-tuple items
#
# @author Vanya BELYAEV Ivan.Belyaev@cern.ch
# @date   2011-06-21
#
#                   $Revision$
# Last modification $Date$
#                by $Author$
# =============================================================================
"""
Decorate the algorithm with proper methods for Pids
"""
# =============================================================================
__version__ = "$Revision: 124897 $"
__author__ = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__ = "2011-06-07"
# =============================================================================
__all__ = (
    # =========================================================================
    # cumulative PID efficiencies
    # =========================================================================
    'pidEff_pions',  # cumulative identification effiicency for pions
    'pidEff_kaons',  # cumulative identification effiicency for kaons
    'pidEff_protons',  # cumulative identification effiicency for protons
    # (di)muon identification efficiency for J/psi
    'pidEff_Jpsi',  # (di)muon identification efficiency for J/psi
    # =========================================================================
    # cumulative MC/DATA correction to track reconstruction efficiency
    # =========================================================================
    'recEff_tracks',  # cumulative MC/DATA correction
    # =========================================================================
    # Acc&Rec&Sel efficiencies
    # =========================================================================
    # (Acc&Rec&Sel&Pid) efficiency for Jpis
    'eff_Jpsi',  # (Acc&Rec&Sel&Pid) efficiency for Jpis
    # (Acc&Rec&Sel)     efficiency for D0
    'eff_D0',  # (Acc&Rec&Sel)     efficiency for D0
    # (Acc&Rec&Sel)     efficiency for D+
    'eff_Dp',  # (Acc&Rec&Sel)     efficiency for D+
    # (Acc&Rec&Sel)     efficiency for Ds+
    'eff_Ds',  # (Acc&Rec&Sel)     efficiency for Ds+
    # (Acc&Rec&Sel)     efficiency for Lambda_c+
    'eff_Lc',  # (Acc&Rec&Sel)     efficiency for Lambda_c+
    # =========================================================================
    # Check if the particle ``is-TOS'' with respect to certain trigger
    # =========================================================================
    'isTos',  # Check if the particle    ``is-TOS''
    'isTis',  # Check if the particle    ``is-TIS''
    'isExTos',  # Check if the particle `` Tos & ~Tis''
    'isExTis',  # Check if the particle ``~Tos &  Tis''
    'isTisTos',  # Check if the particle `` Tis &  Tos''
    'isTob',  # Check if the particle ``~Tis & ~Tos''
    # =========================================================================
    # trigger efficiencies  (tos&trg)
    # =========================================================================
    'trgEff_Jpsi',  # trigger efficiencies  (tos&trg) for J/psi
    'trgEff_D0',  # trigger efficiencies  (tos&trg) for D0
    'trgEff_Dp',  # trigger efficiencies  (tos&trg) for D+
    'trgEff_Ds',  # trigger efficiencies  (tos&trg) for Ds+
    'trgEff_Lc',  # trigger efficiencies  (tos&trg) for Lambda_c+
    # =========================================================================
)
# ========================================================================
import DiCharm.Efficiency as Eff
from AnalysisPython.PyRoUts import VE
# ========================================================================


_eff_Jpsi_ = Eff.pidEff_Jpsi
_eff_pion_ = Eff.pidEff_pion
_eff_kaon_ = Eff.pidEff_kaon
_eff_proton_ = Eff.pidEff_proton
_eff_track_ = Eff.trackEff_corr

# =============================================================================


nBest = lambda s: s.nBest_rec

# =============================================================================
# get the accumulated identification efficiency for all pions


def pidEff_pions(item, suffix=''):
    """
    Get the accumulated identification efficiency for all pions 
    """
    n_pion = getattr(item, 'n_pion' + suffix)

    n_best = nBest(item)

    eff = VE(1, 0)
    for pion in range(0, n_pion):
        #
        p = getattr(item,   'p_pion' + suffix)[pion]
        eta = getattr(item, 'eta_pion' + suffix)[pion]
        #
        eff *= _eff_pion_(p, eta, n_best)

    return eff

# =============================================================================
# get the accumulated identification efficiency for all kaons


def pidEff_kaons(item, suffix=''):
    """
    Get the accumulated identification efficiency for all kaons
    """
    n_kaon = getattr(item, 'n_kaon' + suffix)

    n_best = nBest(item)

    eff = VE(1, 0)
    for kaon in range(0, n_kaon):
        #
        p = getattr(item,   'p_kaon' + suffix)[kaon]
        eta = getattr(item, 'eta_kaon' + suffix)[kaon]
        #
        eff *= _eff_kaon_(p, eta, n_best)

    return eff

# =============================================================================
# get the accumulated identification efficiency for all protons


def pidEff_protons(item, suffix=''):
    """
    Get the accumulated identification efficiency for all protons
    """
    n_proton = getattr(item, 'n_proton' + suffix)

    n_best = nBest(item)

    eff = VE(1, 0)
    for proton in range(0, n_proton):
        #
        p = getattr(item,   'p_proton' + suffix)[proton]
        eta = getattr(item, 'eta_proton' + suffix)[proton]
        #
        eff *= _eff_proton_(p, eta, n_best)

    return eff

# =============================================================================
# get the (di)muon ideintification efficieincy for J/psi


def pidEff_Jpsi(item, suffix='_psi'):
    """
    Get the (di)muon ideintification efficieincy for J/psi 
    """

    pt = getattr(item,   'pt' + suffix)
    y = getattr(item,    'y' + suffix)

    return _eff_Jpsi_(pt, y)


# =============================================================================
# get accumulated  MC/DATA correction to track  recontruction efficiency for
#  all tracks
def recEff_tracks(item, suffix=''):
    """
    Get accumulated  MC/DATA correction to track  recontruction
    efficiency for all tracks 
    """
    n_tracks = getattr(item, 'n_track' + suffix)

    n_best = nBest(item)

    eff = VE(1, 0)
    for track in range(0, n_tracks):
        #
        p = getattr(item,   'p_track' + suffix)[track]
        eta = getattr(item, 'eta_track' + suffix)[track]
        #
        eff *= _eff_track_(p, eta, n_best)

    return eff

    ## n_hadron  = getattr ( item , 'n_pion'   + suffix )
    ## n_hadron += getattr ( item , 'n_kaon'   + suffix )
    ## n_hadron += getattr ( item , 'n_proton' + suffix )

    # eH = VE ( 1 , ( n_hadron * 2.0/100 )**2 ) ## hadronic uncertainties
    # eA = VE ( 1 , ( n_tracks * 0.7/100 )**2 ) ## correlated error per track

    # return eff * eA * eH

# =============================================================================
# get (Acc&Rec&Sel&Pid) efficiency for Jpis


def eff_Jpsi(item, suffix='_psi'):
    """
    Get (Acc&Rec&Sel&Pid) efficiency for J/psi
    """

    pt = getattr(item,   'pt' + suffix)
    y = getattr(item,    'y' + suffix)
    lv01 = getattr(item, 'lv01' + suffix)

    eff = VE(1, 0)
    eff *= Eff.genEff_Jpsi(pt, y)  # acceptance
    #
    eff *= Eff.selEff_Jpsi_3D(pt, y, lv01)  # reco & selection
    # eff *= Eff.selEff_Jpsi_2D ( pt , y ) ## reco & selection
    #
    eff *= Eff.pidEff_Jpsi(pt, y)  # pid

    return eff

# =============================================================================
# get (Acc&Rec&Sel) efficiency for D0


def eff_D0(item, suffix='_D0'):
    """
    Get (Acc&Rec&Sel) efficiency for D0
    
    >>> eff = efF_D0 ( item )
    
    """

    pt = getattr(item,   'pt' + suffix)
    y = getattr(item,    'y' + suffix)

    eff = VE(1, 0)
    eff *= Eff.genEff_D0(pt, y)  # acceptance
    eff *= Eff.selEff_D0(pt, y)  # reco & selection

    return eff

# =============================================================================
# get (Acc&Rec&Sel) efficiency for D+


def eff_Dp(item, suffix='_Dp'):
    """
    Get (Acc&Rec&Sel) efficiency for D+
    
    >>> eff = efF_Dp ( item )
    
    """

    pt = getattr(item,   'pt' + suffix)
    y = getattr(item,    'y' + suffix)

    eff = VE(1, 0)
    eff *= Eff.genEff_Dp(pt, y)  # acceptance
    eff *= Eff.selEff_Dp(pt, y)  # reco & selection

    return eff

# =============================================================================
# get (Acc&Rec&Sel) efficiency for Ds+


def eff_Ds(item, suffix='_Ds'):
    """
    Get (Acc&Rec&Sel) efficiency for Ds+
    
    >>> eff = efF_Ds ( item )
    
    """

    pt = getattr(item,   'pt' + suffix)
    y = getattr(item,    'y' + suffix)

    eff = VE(1, 0)
    eff *= Eff.genEff_Ds(pt, y)  # acceptance
    eff *= Eff.selEff_Ds(pt, y)  # reco & selection

    return eff

# =============================================================================
# get (Acc&Rec&Sel) efficiency for Lambda_c+


def eff_Lc(item, suffix='_Lc'):
    """
    Get (Acc&Rec&Sel) efficiency for Lambda_c+

    >>> eff = efF_Lc ( item )
    
    """

    pt = getattr(item,   'pt' + suffix)
    y = getattr(item,    'y' + suffix)

    eff = VE(1, 0)
    eff *= Eff.genEff_Lc(pt, y)  # acceptance
    eff *= Eff.selEff_Lc(pt, y)  # reco & selection

    return eff


# ========================================================================
# Check if the particle is TOS with respect to certain trigger
#
def isTos(item, part, L0=True, L1=True, L2=True):
    """
    Check if the particle is ``TOS''
    
    >>> tos = isTos ( item , 'psi_' )
    
    """
    #
    tos = True
    #
    if L0:
        tos = tos and 1 == getattr(item, part + 'l0tos_1')
    if L1:
        tos = tos and 1 == getattr(item, part + 'l1tos_1')
    if L2:
        tos = tos and 1 == getattr(item, part + 'l2tos_1')
    #
    return tos

# ========================================================================
# Check if the particle is TIS with respect to certain trigger


def isTis(item, part, L0=True, L1=True, L2=True):
    """
    Check if the particle is ``TIS''
    
    >>> tis = isTis ( item , 'psi_' )
    
    """
    #
    tis = True
    #
    if L0:
        tis = tis and 1 == getattr(item, part + 'l0tis_2')
    if L1:
        tis = tis and 1 == getattr(item, part + 'l1tis_2')
    if L2:
        tis = tis and 1 == getattr(item, part + 'l2tis_2')
    #
    return tis

# ========================================================================
# Check if the particle is ``Ex-TOS'' with respect to certain trigger


def isExTos(item, part, L0=True, L1=True, L2=True):
    """
    Check if the particle is ``Ex-TOS''
    
    >>> exTos = isExTos ( item , 'psi_' )
    
    """
    #
    return isTos(item, part, L0, L1, L2) and not isTis(item, part, L0, L1, L2)

# ========================================================================
# Check if the particle is ``Ex-TIS'' with respect to certain trigger


def isExTis(item, part, L0=True, L1=True, L2=True):
    """
    Check if the particle is ``Ex-Tis''
    
    >>> exTis = isExTis ( item , 'psi_' )
    
    """
    #
    return isTis(item, part, L0, L1, L2) and not isTos(item, part, L0, L1, L2)

# ========================================================================
# Check if the particle is ``Tis&Tos'' with respect to certain trigger


def isTisTos(item, part, L0=True, L1=True, L2=True):
    """
    Check if the particle is ``Tis&Tos''
    
    >>> tistos = isTisTos ( item , 'psi_' )
    
    """
    #
    return isTis(item, part, L0, L1, L2) and isTos(item, part, L0, L1, L2)

# =============================================================================
# Check if the particle is ``!Tis&!Tos'' with respect to certain trigger


def isTob(item, part, L0=True, L1=True, L2=True):
    """
    Check if the particle is ``!Tis&!Tos''
    
    >>> tob = isTos ( item , 'psi_' )
    
    """
    #
    return not isTis(item, part, L0, L1, L2) and not isTos(item, part, L0, L1, L2)

# =============================================================================


# =============================================================================
# get trigger efficiency
def _trgEff_(item,
             suffix,
             data,
             tos,
             L0,
             L1,
             L2):
    """
    Get the trigger efficiency for 
    """
    pt = getattr(item,   'pt' + suffix)
    y = getattr(item,    'y' + suffix)

    if L0 and L1 and L2:
        return data(pt, y, 'L0xL1xL2', tos)
    elif L0 and L1:
        return data(pt, y, 'L0xL1', tos)
    elif L0:
        return data(pt, y, 'L0', tos)
    elif L1:
        return data(pt, y, 'L1', tos)
    elif L2:
        return data(pt, y, 'L2', tos)

    #
    return VE(0, 0)

# =============================================================================
# get trigger efficiency for J/psi


def trgEff_Jpsi(item,
                suffix='_psi',
                tos=True,
                L0=True,
                L1=True,
                L2=True):
    """
    Get the trigger efficiency for Jpsi

    >>> eff = trgEff_Jpsi ( item , L0 = True , L1 = False , L2 = False )
    
    """
    return _trgEff_(item, suffix, Eff.tosEff_Jpsi, tos, L0, L1, L2)

# =============================================================================
# get trigger efficiency for D0


def trgEff_D0(item,
              suffix='_D0',
              tos=True,
              L0=True,
              L1=True,
              L2=True):
    """
    Get the trigger efficiency for D0

    >>> eff = trgEff_D0 ( item , L0 = True , L1 = False , L2 = False )
    
    """
    return _trgEff_(item, suffix, Eff.tosEff_D0, tos, L0, L1, L2)


# =============================================================================
# get trigger efficiency for D+
def trgEff_Dp(item,
              suffix='_Dp',
              tos=True,
              L0=True,
              L1=True,
              L2=True):
    """
    Get the trigger efficiency for D+

    >>> eff = trgEff_Dp ( item , L0 = True , L1 = False , L2 = False )
    
    """
    return _trgEff_(item, suffix, Eff.tosEff_Dp, tos, L0, L1, L2)


# =============================================================================
# get trigger efficiency for Ds+
def trgEff_Ds(item,
              suffix='_Ds',
              tos=True,
              L0=True,
              L1=True,
              L2=True):
    """
    Get the trigger efficiency for Ds+

    >>> eff = trgEff_Ds ( item , L0 = True , L1 = False , L2 = False )
    
    """
    return _trgEff_(item, suffix, Eff.tosEff_Ds, tos, L0, L1, L2)

# =============================================================================
# get trigger efficiency for Lambda_c+


def trgEff_Lc(item,
              suffix='_Lc',
              tos=True,
              L0=True,
              L1=True,
              L2=True):
    """
    Get the trigger efficiency for Lambda_c+

    >>> eff = trgEff_Lc ( item , L0 = True , L1 = False , L2 = False )
    
    """
    return _trgEff_(item, suffix, Eff.tosEff_Lc, tos, L0, L1, L2)


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
