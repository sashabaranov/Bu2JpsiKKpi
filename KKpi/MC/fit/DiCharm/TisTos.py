#!/usr/bin/env python
# =============================================================================
# @file DiCharm/TisTos.py
#
# Helper module for TisTos'ing
#
# @author Vanya BELYAEV Ivan.Belyaev@cern.ch
# @date   2011-06-21
#
#                   $Revision$
# Last modification $Date$
#                by $Author$
# =============================================================================
"""
Helper module for TisTos'ing 
"""
# =============================================================================
__version__ = "$Revision: 124897 $"
__author__ = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__ = "2011-06-07"
# =============================================================================
__all__ = (
    'lines',  # the lines used for TosTos
)
# ========================================================================
from LoKiCore.basic import cpp
# ========================================================================
lines = {}
lines["D0"] = {}
lines["D*+"] = {}
lines["D+"] = {}
lines["Ds+"] = {}
lines["Lc+"] = {}
lines["psi"] = {}
lines["psi1"] = {}  # clean-selection: dimuon
lines["psi2"] = {}  # clean-selection: detached
lines["psi3"] = {}  # clean-selection: unbiased
lines["psi'"] = {}
lines["Upsilon"] = {}
lines["Upsilon1"] = {}  # clean selection
lines["W+"] = {}
lines["Z0"] = {}
#
# D0:
#
lines["D0"]['L0TOS'] = 'L0HadronDecision'
lines["D0"]['L0TIS'] = 'L0(Hadron|DiMuon|Muon|Electron|Photon)Decision'
lines["D0"]['Hlt1TOS'] = 'Hlt1TrackAllL0Decision'
lines["D0"]['Hlt1TIS'] = 'Hlt1(Track|Muon|SingleMuon|DiMion).*Decision'
lines["D0"]['Hlt2TOS'] = 'Hlt2CharmHadD02HH.*Decision'
lines["D0"]['Hlt2TIS'] = 'Hlt2(Charm|Topo|DiMuon|Single).*Decision'
#
# D*+
#
lines["D*+"]['L0TOS'] = lines['D0']['L0TOS']
lines["D*+"]['L0TIS'] = lines['D0']['L0TIS']
lines["D*+"]['Hlt1TOS'] = lines['D0']['Hlt1TOS']
lines["D*+"]['Hlt1TIS'] = lines['D0']['Hlt1TIS']
lines["D*+"]['Hlt2TOS'] = 'Hlt2CharmHadD.*Decision'  # ATTENTION
lines["D*+"]['Hlt2TIS'] = lines['D0']['Hlt2TIS']
#
# D+
#
lines["D+"]['L0TOS'] = lines['D0']['L0TOS']
lines["D+"]['L0TIS'] = lines['D0']['L0TIS']
lines["D+"]['Hlt1TOS'] = lines['D0']['Hlt1TOS']
lines["D+"]['Hlt1TIS'] = lines['D0']['Hlt1TIS']
lines["D+"]['Hlt2TOS'] = 'Hlt2CharmHadD.*Decision'  # ATTENTION
lines["D+"]['Hlt2TIS'] = lines['D0']['Hlt2TIS']
#
# Ds+
#
lines["Ds+"]['L0TOS'] = lines['D0']['L0TOS']
lines["Ds+"]['L0TIS'] = lines['D0']['L0TIS']
lines["Ds+"]['Hlt1TOS'] = lines['D0']['Hlt1TOS']
lines["Ds+"]['Hlt1TIS'] = lines['D0']['Hlt1TIS']
lines["Ds+"]['Hlt2TOS'] = lines['D+']['Hlt2TOS']  # ATTENTION
lines["Ds+"]['Hlt2TIS'] = lines['D0']['Hlt2TIS']
#
# Lambda_c+
#
lines["Lc+"]['L0TOS'] = lines['D0']['L0TOS']
lines["Lc+"]['L0TIS'] = lines['D0']['L0TIS']
lines["Lc+"]['Hlt1TOS'] = lines['D0']['Hlt1TOS']
lines["Lc+"]['Hlt1TIS'] = lines['D0']['Hlt1TIS']
lines["Lc+"]['Hlt2TOS'] = lines['D+']['Hlt2TOS']  # ATTENTION
lines["Lc+"]['Hlt2TIS'] = lines['D0']['Hlt2TIS']
#
# J/psi
#
lines["psi"]['L0TOS'] = 'L0(DiMuon|Muon)Decision'
lines["psi"]['L0TIS'] = 'L0(Hadron|DiMuon|Muon|Electron|Photon)Decision'
lines["psi"]['Hlt1TOS'] = 'Hlt1(DiMuon|SingleMuon|TrackMuon).*Decision'
lines["psi"]['Hlt1TIS'] = 'Hlt1(DiMuon|SingleMuon|Track).*Decision'
lines["psi"]['Hlt2TOS'] = 'Hlt2(DiMuon|ExpressJPsi|SingleMuon).*Decision'
lines["psi"]['Hlt2TIS'] = 'Hlt2(Charm|Topo|DiMuon|Single).*Decision'
#
## J/psi ("clean-1")
#
lines["psi1"]['L0TOS'] = 'L0(DiMuon|Muon)Decision'
lines["psi1"]['L0TIS'] = 'L0(Hadron|DiMuon|Muon|Electron|Photon)Decision'
lines["psi1"]['Hlt1TOS'] = 'Hlt1(DiMuon|TrackMuon).*Decision'
lines["psi1"]['Hlt1TIS'] = 'Hlt1(DiMuon|SingleMuon|Track).*Decision'
lines["psi1"]['Hlt2TOS'] = 'Hlt2DiMuon.*Decision'
lines["psi1"]['Hlt2TIS'] = 'Hlt2(Charm|Topo|DiMuon|Single).*Decision'
#
## J/psi ("clean-2")
#
lines["psi2"]['L0TOS'] = 'L0(DiMuon|Muon)Decision'
lines["psi2"]['L0TIS'] = 'L0(Hadron|DiMuon|Muon|Electron|Photon)Decision'
lines["psi2"]['Hlt1TOS'] = 'Hlt1DiMuonHighMass.*Decision'
lines["psi2"]['Hlt1TIS'] = 'Hlt1(DiMuon|SingleMuon|Track).*Decision'
lines["psi2"]['Hlt2TOS'] = 'Hlt2DiMuonDetached(Heavy|JPsi)Decision'
lines["psi2"]['Hlt2TIS'] = 'Hlt2(Charm|Topo|DiMuon|Single).*Decision'
#
## J/psi ("unbiased")
#
lines["psi3"]['L0TOS'] = 'L0(DiMuon|Muon)Decision'
lines["psi3"]['L0TIS'] = 'L0(Hadron|DiMuon|Muon|Electron|Photon)Decision'
lines["psi3"]['Hlt1TOS'] = 'Hlt1DiMuonHighMass.*Decision'
lines["psi3"]['Hlt1TIS'] = 'Hlt1(DiMuon|SingleMuon|Track).*Decision'
lines["psi3"]['Hlt2TOS'] = 'Hlt2DiMuonJPsiHighPT.*Decision'
lines["psi3"]['Hlt2TIS'] = 'Hlt2(Charm|Topo|DiMuon|Single).*Decision'

#
#
# Psi'
#
import copy
lines["psi'"] = copy.deepcopy(lines["psi"])
lines["psi'"]['Hlt2TOS'] = 'Hlt2(DiMuon|SingleMuon).*Decision'

#
# Upsilon
#
lines["Upsilon"]['L0TOS'] = 'L0(DiMuon|Muon)Decision'
lines["Upsilon"]['L0TIS'] = 'L0(Hadron|DiMuon|Muon|Electron|Photon)Decision'
lines["Upsilon"]['Hlt1TOS'] = 'Hlt1(DiMuon|SingleMuon|TrackMuon).*Decision'
lines["Upsilon"]['Hlt1TIS'] = 'Hlt1(DiMuon|SingleMuon|Track).*Decision'
lines["Upsilon"]['Hlt2TOS'] = 'Hlt2(DiMuon|SingleMuonHighPT).*Decision'
lines["Upsilon"]['Hlt2TIS'] = 'Hlt2(Charm|Topo|DiMuon|Single).*Decision'
#
# ``clean-Y''
lines["Upsilon1"]['L0TOS'] = 'L0DiMuonDecision'
lines["Upsilon1"][
    'L0TIS'] = 'L0(Hadron|DiMuon|Muon|Electron|Photon)Decision'
lines["Upsilon1"]['Hlt1TOS'] = 'Hlt1DiMuonHighMassDecision'
lines["Upsilon1"]['Hlt1TIS'] = 'Hlt1(DiMuon|SingleMuon|Track).*Decision'
lines["Upsilon1"]['Hlt2TOS'] = 'Hlt2DiMuonBDecision'
lines["Upsilon1"]['Hlt2TIS'] = 'Hlt2(Charm|Topo|DiMuon|Single).*Decision'

#
# Z0
#
lines["Z0"]['L0TOS'] = 'L0MuonDecision'
lines["Z0"]['L0TIS'] = 'L0(Hadron|DiMuon|Muon|Electron|Photon)Decision'
lines["Z0"]['Hlt1TOS'] = 'Hlt1SingleMuonHighPTDecision'
lines["Z0"]['Hlt1TIS'] = 'Hlt1(DiMuon|SingleMuon|Track).*Decision'
lines["Z0"]['Hlt2TOS'] = 'Hlt2SingleMuonHighPTDecision'
lines["Z0"]['Hlt2TIS'] = 'Hlt2(Charm|Topo|DiMuon|Single).*Decision'


#
# W+/W-
#
lines["W+"]['L0TOS'] = 'L0MuonDecision'
lines["W+"]['L0TIS'] = 'L0(Hadron|DiMuon|Muon|Electron|Photon)Decision'
lines["W+"]['Hlt1TOS'] = 'Hlt1SingleMuonHighPTDecision'
lines["W+"]['Hlt1TIS'] = 'Hlt1(DiMuon|SingleMuon|Track).*Decision'
lines["W+"]['Hlt2TOS'] = 'Hlt2SingleMuon(HighPT|VHighPT)Decision'
lines["W+"]['Hlt2TIS'] = 'Hlt2(Charm|Topo|DiMuon|Single).*Decision'

lines["W-"] = copy.deepcopy(lines["W+"])

#
# Subset of J/psi-detached
#
lines["psi_clean"] = {}
lines["psi_clean"]['L0TOS'] = 'L0(DiMuon|Muon)Decision'
lines["psi_clean"][
    'L0TIS'] = 'L0(Hadron|DiMuon|Muon|Electron|Photon)Decision'
lines["psi_clean"]['Hlt1TOS'] = 'Hlt1DiMuonHighMassDecision'
lines["psi_clean"]['Hlt1TIS'] = 'Hlt1(DiMuon|SingleMuon|Track).*Decision'
lines["psi_clean"]['Hlt2TOS'] = 'Hlt2DiMuonDetached(Heavy|JPsi)Decision'
lines["psi_clean"]['Hlt2TIS'] = 'Hlt2(Charm|Topo|DiMuon|Single).*Decision'


#
# =============================================================================
# The actual job steering
if '__main__' == __name__:

    print 80 * '*'
    print __doc__
    print ' Author  : ', __author__
    print ' Version : ', __version__
    print ' Date    : ', __date__
    print 80 * '*'


# =============================================================================
# The END
# =============================================================================
