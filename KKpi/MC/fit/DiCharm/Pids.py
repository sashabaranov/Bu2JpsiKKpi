#!/usr/bin/env python
# =============================================================================
# @file DiCharm/Pids.py
#
# Decorate the algorithm with proper methods for Pids
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
    'treatPions',  # add information about pions
    'treatKaons',  # add information about kaons
    'treatProtons',  # add information about protons
    'treatMuons',  # add information about muons
    'treatTracks',  # add information about the tracks
)
# ========================================================================
from Bender.Main import LoKi, SUCCESS
from LoKiPhys.decorators import *
from LoKiTracks.decorators import TrKEY
from GaudiKernel.SystemOfUnits import GeV, MeV, mm

# ========================================================================
# add pion information into n-tuple
#  @param tup   n-tuple
#  @param p     the particle


def treatPions(self,
               tup,
               p,
               suffix=''):

    #
    if hasattr(p, 'particle'):
        p = p.particle()
    #
    # get all pions from decay tree
    #
    good = LHCb.Particle.ConstVector()
    p.children(self._pions, good)
    #
    sc = tup.column('mindll_piK' + suffix, self._minPi(p))
    #
    tup.fArrayP('phi_pion' + suffix, PHI, LHCb.Particle.Range(good),
                'n_pion' + suffix, 10)

    return tup.fArrayP('p_pion' + suffix, P / GeV,
                       'pt_pion' + suffix, PT / GeV,
                       'eta_pion' + suffix, ETA,
                       #                         'phi_pion' + suffix          , PHI          ,
                       'pid_pion' + suffix, PIDpi - PIDK,
                       LHCb.Particle.Range(good),
                       'n_pion' + suffix, 10)


# ========================================================================
# add kaon information into n-tuple
#  @param tup   n-tuple
#  @param p     the particle
def treatKaons(self,
               tup,
               p,
               suffix=''):

    #
    if hasattr(p, 'particle'):
        p = p.particle()
    #
    # get all kaons form decay tree
    #
    good = LHCb.Particle.ConstVector()
    p.children(self._kaons, good)
    #
    sc = tup.column('mindll_K' + suffix, self._minK(p))
    #
    tup.fArrayP('phi_kaon' + suffix, PHI, LHCb.Particle.Range(good),
                'n_kaon' + suffix, 10)
    return tup.fArrayP('p_kaon' + suffix, P / GeV,
                       'pt_kaon' + suffix, PT / GeV,
                       'eta_kaon' + suffix, ETA,
                       'pid_kaon' + suffix, PIDK - PIDpi,
                       LHCb.Particle.Range(good),
                       'n_kaon' + suffix, 10)


# ========================================================================
# add proton information into n-tuple
#  @param tup   n-tuple
#  @param p     the particle
def treatProtons(self,
                 tup,
                 p,
                 suffix=''):

    #
    if hasattr(p, 'particle'):
        p = p.particle()
    #
    # get all protons form decay tree
    #
    good = LHCb.Particle.ConstVector()
    p.children(self._protons, good)
    #
    sc = tup.column('mindll_pK' + suffix, self._minPK(p))
    sc = tup.column('mindll_ppi' + suffix, self._minPpi(p))
    #
    tup.fArrayP('p_proton' + suffix, P / GeV,
                'pt_proton' + suffix, PT / GeV,
                'eta_proton' + suffix, ETA,
                LHCb.Particle.Range(good),
                'n_proton' + suffix, 10)
    return tup.fArrayP('pid_proton_pi' + suffix, PIDp - PIDpi,
                       'pid_proton_K' + suffix, PIDp - PIDK,
                       LHCb.Particle.Range(good),
                       'n_proton' + suffix, 10)


# ========================================================================
# add muon information into n-tuple
#  @param tup   n-tuple
#  @param p     the particle
def treatMuons(self,
               tup,
               p,
               suffix=''):

    #
    if hasattr(p, 'particle'):
        p = p.particle()
    #
    # get all muons from decay tree
    #
    good = LHCb.Particle.ConstVector()
    p.children(self._muons, good)
    #
    sc = tup.column('mindll_mu' + suffix, self._minMu(p))
    sc = tup.column('minpt_mu' + suffix, self._minPt(p))
    sc = tup.column('maxEtC_mu' + suffix, self._maxEtC(p))
    sc = tup.column('maxPtC_mu' + suffix, self._maxPtC(p))
    sc = tup.column('maxEcalE_mu' + suffix, self._maxEcalE(p))
    sc = tup.column('maxHcalE_mu' + suffix, self._maxHcalE(p))
    #
    sc = tup.fArrayP('p_mu' + suffix, P / GeV,
                     'pt_mu' + suffix, PT / GeV,
                     'eta_mu' + suffix, ETA,
                     'pid_mu' + suffix, PIDmu - PIDpi,
                     LHCb.Particle.Range(good),
                     'n_muon' + suffix,  10)
    return tup.fArrayP('ptC_mu' + suffix, self._EtC,
                       'etC_mu' + suffix, self._PtC,
                       'eEcal_mu' + suffix, self._EcalE,
                       'eHcal_mu' + suffix, self._HcalE,
                       LHCb.Particle.Range(good),
                       'n_muon' + suffix, 10)


# ========================================================================
# add tracks information into n-tuple
#  @param tup   n-tuple
#  @param p     the particle
def treatTracks(self,
                tup,
                p,
                suffix=''):

    #
    if hasattr(p, 'particle'):
        p = p.particle()
    #
    # get all tracks from decay tree
    #
    from math import sqrt, acos

    good = LHCb.Particle.ConstVector()
    p.children(self._tracks, good)
    #
    m2min = 2.0
    for i in range(0, good.size()):
        p_i = good[i]
        q_i = Q(p_i)

        for j in range(i + 1, good.size()):
            p_j = good[j]
            q_j = Q(p_j)

            if q_i * q_j < 0:
                continue

            m2 = self._delta_m2(p_i, p_j)
            m2min = min(m2min, m2)
    sc = SUCCESS
    sc = tup.column('m2min_track' + suffix, m2min)
    sc = tup.column('minPt_track' + suffix, self._minPt(p))
    sc = tup.column('minEta_track' + suffix, self._minEta(p))
    sc = tup.column('maxEta_track' + suffix, self._maxEta(p))
    sc = tup.column('maxChi2_track' + suffix, self._maxTrChi2(p))
    sc = tup.column('maxTrGh_track' + suffix, self._maxTrGhost(p))

    sc = tup.fArrayP('p_track' + suffix, P / GeV,
                     'pt_track' + suffix, PT / GeV,
                     'eta_track' + suffix, ETA,
                     #                         'chi1_track' + suffix        , TRCHI2DOF  ,
                     'PChi2_track' + suffix, TRPCHI2,
                     LHCb.Particle.Range(good),
                     'n_track' + suffix, 10)
    sc = tup.fArrayP('chi2_track' + suffix, TRCHI2DOF,
                     LHCb.Particle.Range(good),
                     'n_track' + suffix, 10)
    return sc

# =============================================================================
# add basic kinematical information into N-tuple


def treatKine(self,
              tup,
              p,
              suffix,
              good=None):
    """
    Add basic kinematical information into N-tuple
    """
    if hasattr(p, 'particle'):
        p = p.particle()
    #
    tup.column('pid' + suffix, int(ID(p)))
    tup.column('pt' + suffix,       PT(p) / GeV)
    tup.column('m' + suffix,        M(p) / GeV)
    tup.column('y' + suffix,        Y(p))
    tup.column('lv01' + suffix, self._lv01(p))
    tup.column('phi' + suffix,      PHI(p))
    tup.column('p4' + suffix,                 p.momentum() / GeV)
    tup.column('ctau' + suffix, self._ctau(p))
    tup.column('dtf' + suffix, self._dtfchi2(p))
    tup.column('c2ip' + suffix, self._ipchi2(p))
    tup.column('dls' + suffix, self._dls(p))
    tup.column('vchi2' + suffix, self._vchi2(p))
    #
    ok = True
    if good:
        ok = good(p)

    tup.column('good' + suffix,  ok, 0, 1)

    return SUCCESS

# =============================================================================
# add some event summary information


def addRecSummary(self,
                  tup,
                  data):
    """
    Add event summary information
    """

    #
    nPV = int(data(LHCb.RecSummary.nPVs))
    nLong = int(data(LHCb.RecSummary.nLongTracks))
    nDown = int(data(LHCb.RecSummary.nDownstreamTracks))
    nUp = int(data(LHCb.RecSummary.nUpstreamTracks))
    nVelo = int(data(LHCb.RecSummary.nVeloTracks))
    nTT = int(data(LHCb.RecSummary.nTTracks))
    nBack = int(data(LHCb.RecSummary.nBackTracks))
    #
    hRich1 = int(data(LHCb.RecSummary.nRich1Hits))
    hRich2 = int(data(LHCb.RecSummary.nRich2Hits))
    hVelo = int(data(LHCb.RecSummary.nVeloClusters))
    hIT = int(data(LHCb.RecSummary.nITClusters))
    hTT = int(data(LHCb.RecSummary.nTTClusters))
    hOT = int(data(LHCb.RecSummary.nOTClusters))
    hSPD = int(data(LHCb.RecSummary.nSPDhits))

    tup.column('nPV', nPV, 0,    30)
    tup.column('nLong', nLong, 0,  1000)
    tup.column('nDown', nDown, 0,   500)
    tup.column('nUp', nUp, 0,   500)
    tup.column('nVelo', nVelo, 0,  1000)
    tup.column('nTT', nTT, 0,   500)
    tup.column('nBack', nBack, 0,   200)

    tup.column('hSPD', hSPD, 0,  2000)
    tup.column('hRich1', hRich1, 0, 30000)
    tup.column('hRich2', hRich2, 0, 30000)
    tup.column('hVelo', hVelo, 0, 60000)
    tup.column('hIT', hIT, 0, 30000)
    tup.column('hTT', hTT, 0, 30000)
    tup.column('hOT', hOT, 0, 30000)

    return SUCCESS

# =============================================================================
# add some event summary information


def addGecInfo(self,
               tup,
               data):
    """
    Add global event counters 
    """
    for key in data:

        v = int(data[key])

        tup.column(key + '_gec', v)

    return SUCCESS

# =============================================================================
# initialize pid/tracks machinery


def _pid_initialize(self):
    """
    Initialize Pid/Tracks machinery 
    """
    self._pions = LoKi.Child.Selector('pi+' == ABSID)
    self._kaons = LoKi.Child.Selector('K+' == ABSID)
    self._protons = LoKi.Child.Selector('p+' == ABSID)
    self._muons = LoKi.Child.Selector('mu+' == ABSID)
    self._tracks = LoKi.Child.Selector(HASTRACK)
    self._basic = LoKi.Child.Selector(ISBASIC)
    #
    from LoKiCore.functions import switch
    self._ctau = switch(ISBASIC, -100, BPVLTIME(9) * c_light)
    self._lv01 = switch(ISBASIC, -100, LV01)
    self._vchi2 = switch(ISBASIC, -100, VFASPF(VCHI2))
    self._dtfchi2 = switch(
        ISBASIC,  BPVIPCHI2() / 2, DTF_CHI2NDOF(True))
    self._ipchi2 = BPVIPCHI2()
    self._dls = switch(
        ISBASIC, -100, LoKi.Particles.DecayLengthSignificanceDV())

    self._minK = MINTREE('K+' == ABSID, PIDK - PIDpi)
    self._minPi = MINTREE('pi+' == ABSID, PIDpi - PIDK)
    self._minPK = MINTREE('p+' == ABSID, PIDp - PIDK)
    self._minPpi = MINTREE('p+' == ABSID, PIDp - PIDpi)
    self._minMu = MINTREE('mu+' == ABSID, PIDmu - PIDpi)
    self._minPt = MINTREE(ISBASIC & HASTRACK, PT) / GeV
    self._minEta = MINTREE(ISBASIC & HASTRACK, ETA)
    self._maxEta = MAXTREE(ISBASIC & HASTRACK, ETA)
    #
    self._maxTrChi2 = MAXTREE(ISBASIC & HASTRACK, TRCHI2DOF)
    self._maxTrGhost = MAXTREE(ISBASIC & HASTRACK, TRGHOSTPROB)
    self._minTrIPchi2 = MINTREE(ISBASIC & HASTRACK, BPVIPCHI2())
    #
    self._trIpChi = switch(ISBASIC & HASTRACK, BPVIPCHI2(), -1000)
    #
    #
    self._EtC = PINFO(55001, -100 * GeV)
    self._PtC = PINFO(55002, -100 * GeV)
    self._maxEtC = MAXTREE('mu+' == ABSID, self._EtC)
    self._maxPtC = MAXTREE('mu+' == ABSID, self._PtC)
    self._EcalE = PPINFO(LHCb.ProtoParticle.CaloEcalE, -100 * GeV)
    self._HcalE = PPINFO(LHCb.ProtoParticle.CaloHcalE, -100 * GeV)
    self._maxEcalE = MAXTREE('mu+' == ABSID, self._EcalE)
    self._maxHcalE = MAXTREE('mu+' == ABSID, self._HcalE)
    #
    self._delta_m2 = LoKi.PhysKinematics.deltaM2
    #
    _pid = LoKi.Particles.pidFromName
    self._dimu_pids = [
        (ADMASS('J/psi(1S)'),  _pid('J/psi(1S)')),
        (ADMASS('psi(2S)'),  _pid('psi(2S)')),
        (ADMASS('Upsilon(1S)'),  _pid('Upsilon(1S)')),
        (ADMASS('Upsilon(2S)'),  _pid('Upsilon(2S)')),
        (ADMASS('Upsilon(3S)'),  _pid('Upsilon(3S)')),
    ]

    return SUCCESS

# =============================================================================
# finalize pid/tracks machinery


def _pid_finalize(self):
    """
    Initialize Pid/Tracks machinery 
    """
    self._pions = None
    self._kaons = None
    self._protons = None
    self._muons = None
    self._tracks = None
    self._basic = None
    #
    self._ctau = None
    self._lv01 = None
    self._vchi2 = None
    self._dtfchi2 = None
    self._minK = None
    self._minPi = None
    self._minPK = None
    self._minPpi = None
    self._minMu = None
    self._minPt = None
    self._maxEta = None
    self._minEta = None
    #
    self._maxTrChi2 = None
    self._maxTrGhost = None
    self._maxTrIpChi2 = None
    self._trIpChi = None
    self._ipchi2 = None
    self._dls = None
    #
    self._EtC = None
    self._PtC = None
    self._maxEtC = None
    self._maxPtC = None
    self._EcalE = None
    self._HcalE = None
    self._maxEcalE = None
    self._maxHcalE = None

    self._dimu_pids = None

    return SUCCESS

# set new particle ID for dimuons,
#  accoring to nearest resonance


def _dimuon_ID_(self, dimu):
    """
    Set new particle ID for dimuons,
    accoring to nearest dimuon resonance     
    """

    pid = LHCb.ParticleID(0)
    m_ = -1
    for i in self._dimu_pids:
        dm = i[0](dimu)
        if dm < m_ or m_ < 0 or 0 == pid.pid():
            m_ = dm
            pid = i[1]
            dimu.setParticleID(pid)

    return pid


# =============================================================================
LoKi.Algo.treatPions = treatPions
LoKi.Algo.treatKaons = treatKaons
LoKi.Algo.treatProtons = treatProtons
LoKi.Algo.treatMuons = treatMuons
LoKi.Algo.treatTracks = treatTracks
LoKi.Algo.treatKine = treatKine

LoKi.Algo.addRecSummary = addRecSummary
LoKi.Algo.addGecInfo = addGecInfo

LoKi.Algo.pid_initialize = _pid_initialize
LoKi.Algo.pid_finalize = _pid_finalize
LoKi.Algo.dimuon_ID = _dimuon_ID_
# =============================================================================


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
