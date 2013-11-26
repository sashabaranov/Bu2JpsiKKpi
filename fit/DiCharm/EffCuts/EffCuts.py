#!/usr/bin/env ipython
# ========================================================================
# @file
#
#  Simple algorithm to get reconstrcution &
#   selection efficiencies for prompt charm decays
#
#  This file is a part of
#  <a href="http://cern.ch/lhcb-comp/Analysis/Bender/index.html">Bender project</a>
#  <b>``Python-based Interactive Environment for Smart and Friendly
#   Physics Analysis''</b>
#
#  The package has been designed with the kind help from
#  Pere MATO and Andrey TSAREGORODTSEV.
#  And it is based on the
#  <a href="http://cern.ch/lhcb-comp/Analysis/LoKi/index.html">LoKi project:</a>
#  ``C++ ToolKit for Smart and Friendly Physics Analysis''
#
#  By usage of this code one clearly states the disagreement
#  with the smear campaign of Dr.O.Callot et al.:
#  ``No Vanya's lines are allowed in LHCb/Gaudi software.''
#
#  @date   2011-06-21
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#
#                    $Revision$
#  Last modification $Date$
#                 by $Author$
# =============================================================================
"""

Simple algorithm to get reconstrcution &
selection efficiencies for prompt charm decays 

This file is a part of BENDER project:
   ``Python-based Interactive Environment for Smart and Friendly Physics Analysis''

The project has been designed with the kind help from Pere MATO and Andrey TSAREGORODTSEV. 

And it is based on the LoKi project:
    ``C++ ToolKit for Smart and Friendly Physics Analysis''

By usage of this code one clearly states the disagreement 
with the smear campaign of Dr.O.Callot et al.: 
    ``No Vanya's lines are allowed in LHCb/Gaudi software.''

"""
# =============================================================================
__author__ = 'Vanya BELYAEV  Ivan.Belyaev@cern.ch'
__date__ = "2011-06-21"
__version__ = '$Revision$'
# =============================================================================
## needed to produce/visualize the histograms
import ROOT
from Bender.MainMC import *  # import all bender goodies
## easy access to various geometry routines
import LHCbMath.Types
from Gaudi.Configuration import *  # needed for job configuration

from GaudiKernel.SystemOfUnits import GeV, MeV, mm
from GaudiKernel.PhysicalConstants import c_light
# =============================================================================
hasB = MCINANCESTORS(BEAUTY)
# =============================================================================
import DiCharm.Pids  # add methods for Pid/Track information
#AlgoMC.pid_initialize = Algo.pid_initialize
#AlgoMC.pid_finalize   = Algo.pid_finalize
# =============================================================================


def bpid(part):
    if BEAUTY(part):
        return part.particleID().pid()
    mother = part.mother()
    while mother:
        if BEAUTY(mother):
            return mother.particleID().pid()
        mother = mother.mother()
    return 0
# ========================================================================
# @class MCD0
#  Simple algorithm to get MC-efficiencies  for D0->K-pi+ events
#  @date   2011-05-27
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class MCD0(AlgoMC):

    """
    Simple algorithm to get MC-efficiencies for D0->K-pi+ events 
    """
    # the only one esential method:

    def initialize(self):

        sc = AlgoMC .     initialize(self)
        if sc.isFailure():
            return sc

        sc = self   . pid_initialize()
        if sc.isFailure():
            return sc

        from DiCharm.GoodParticles import (basicD0,
                                           basicJpsi,
                                           basicDp,
                                           basicDs,
                                           basicLc,
                                           prompt)

        self._goodD0 = in_range(2, Y, 4.5) & basicD0() & prompt()
        self._goodJpsi = in_range(2, Y, 4.5) & basicJpsi() & prompt()
        self._goodDp = in_range(2, Y, 4.5) & basicDp() & prompt()
        self._goodDs = in_range(2, Y, 4.5) & basicDs() & prompt()
        self._goodLc = in_range(2, Y, 4.5) & basicLc() & prompt()

        self.k_plus = LoKi.MCChild.Selector(
            ' Charm ==> ^K+  K-  ( pi+ | pi- ) ')
        self.k_minus = LoKi.MCChild.Selector(
            ' Charm ==>  K+ ^K-  ( pi+ | pi- ) ')

        self._ctau = BPVLTIME(9) * c_light
        self._ipchi2 = BPVIPCHI2()
        self._vxchi2 = VFASPF(VCHI2)

        return SUCCESS

    # finalize it!
    def finalize(self):

        self.dumpHistos()

        self._goodD0 = None
        self._goodJpsi = None
        self._goodDp = None
        self._goodDs = None
        self._goodLc = None

        self.k_plus = None
        self.k_minus = None

        self . pid_finalize()

        self._ctau = None
        self._ipchi2 = None
        self._vxchi2 = None

        return AlgoMC.finalize(self)

    # ========================================================================
    def treatRC(self,
                tup,
                parts,
                suffix=''):
        """
        Add basic kinematical information into N-tuple
        """
        tup.fArrayP(
            'minK' + suffix, self._minK,
            'minPi' + suffix, self._minPi,
            'minPK' + suffix, self._minPK,
            'minPpi' + suffix, self._minPpi,
            parts,
            'n_rc' + suffix, 10)
        tup.fArrayP(
            'mass' + suffix, M / GeV,
            'minPT' + suffix, self._minPt,
            'minEta' + suffix, self._minEta,
            'maxEta' + suffix, self._maxEta,
            parts,
            'n_rc' + suffix, 10)
        tup.fArrayP(
            'ctau' + suffix, self._ctau,
            'dtf' + suffix, self._dtfchi2,
            'ipchi2' + suffix, self._ipchi2,
            'vxchi2' + suffix, self._vxchi2,
            parts,
            'n_rc' + suffix, 10)
        return tup.fArrayP(
            'lv01' + suffix, LV01,
            parts,
            'n_rc' + suffix, 10)

    # the only one esential method:
    def analyse(self):
        """
        The only one essential method
        """
        beauty = self.mcselect('Beauty',  BEAUTY)
        beauty_g = self. gselect('gBeauty', GBEAUTY)
        # if not beauty.empty() :
        #    return self.Warning('Skip event with beauty-production', SUCCESS )

        # get true MC-decays
        mcD0 = self.mcselect(
            'mcD0', '[ D*(2010)+ -> ^( D0 => K-  pi+ )  pi+]CC')
        if mcD0   .empty():
            return self.Warning('No True MC-decays are found', SUCCESS)

        goodMC = self.mcselect('goodMC',
                               mcD0,
                               in_range(1.95, MCY, 4.55) & (MCPT > 1.95 * GeV))

        if goodMC.empty():
            return self.Warning('No Good MC-decays are found', SUCCESS)

        if 1 != len(goodMC):
            self.Warning('Too many MC-decays are found %s' %
                         goodMC.size(), SUCCESS)

        # keep only the first ``good'' particle
        D0mc = goodMC[0]

        mcK_ = self.mcselect(
            'mcK_', '[ D*(2010)+ ->  ( D0 =>  ^K- pi+ )  pi+]CC')
        if mcK_   .empty():
            return self.Warning('No True MC-decays (K) are found', SUCCESS)
        mcK = self.mcselect('mcK', mcK_,  FROMMCTREE(D0mc))
        if mcK_   .empty():
            return self.Warning('No Good MC-decays (K) are found', SUCCESS)
        #
        mcPi_ = self.mcselect(
            'mcPi_', '[ D*(2010)+ ->  ( D0 =>   K- ^pi+ )  pi+]CC')
        if mcPi_  .empty():
            return self.Warning('No True MC-decays (pi) are found', SUCCESS)
        mcPi = self.mcselect('mcPi', mcPi_, FROMMCTREE(D0mc))
        if mcPi   .empty():
            return self.Warning('No Good MC-decays (pi) are found', SUCCESS)

        if not beauty.empty() or hasB(D0mc):
            mother = D0mc.mother()
            while mother:
                if BEAUTY(mother):
                    break
                mother = mother.mother()
            if mother and not hasB(D0mc):
                self.Error('crazyyyyy', SUCCESS, 100)
                print 'found beauty mother: ', mother.decay(), hasB(D0mc)

        tup = self.nTuple('MCD0')

        tup.column('hasBevt', not beauty   .empty(), 0, 1)
        tup.column('hasBevtG', not beauty_g .empty(), 0, 1)
        tup.column('hasB', hasB(D0mc), 0, 1)

        tup.column('mcpt', MCPT(D0mc) / GeV)
        tup.column('mcy', MCY(D0mc))

        # create the MC-truth matching functors:
        truePi = MCTRUTH(self.mcTruth(), mcPi)
        trueK = MCTRUTH(self.mcTruth(), mcK)
        trueD0 = MCTRUTH(self.mcTruth(), D0mc)

        # get reconstructed pions
        pions = self.select('pion', ('pi+' == ABSID) & truePi)
        # get reconstructed kaons
        kaons = self.select('kaon', ('K+' == ABSID) & trueK)

        d0 = self.loop('kaon pion', 'D0')
        for d in d0:

            k = d(1)
            p = d(2)

            if Q(k) * Q(p) > 0:
                continue

            if not 1.7 * GeV < d.mass(1, 2) < 2.0 * GeV:
                continue

            if not 0 <= VCHI2(d) < 1000:
                continue

            if not trueD0(d):
                continue

            if Q(k) < 0:
                d.setPID('D0')
            else:
                d.setPID('D~0')

            self.plot(M(d) / GeV, 'mass D0->K-pi+ ', 1.8, 1.9)

            d0.save('d0')

        drec = self.selected('d0')
        dsel = self.select('D0', drec, self._goodD0)
        self.treatRC(tup, drec)

        tup.column('nDrec', drec . size())
        tup.column('nDsel', dsel . size())

        tup.write()

        self.setFilterPassed(not dsel.empty())

        return SUCCESS

# ========================================================================
# @class MCPsi
#  Simple algorithm to get MC-efficiencies  for J/psi ->mu+ mu- events
#  @date   2011-05-27
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class MCPsi(MCD0):

    """
    Simple algorithm to get MC-efficiencies for J/psi -> mu+ mu- events 
    """
    # the only one esential method:

    def analyse(self):
        """
        The only one essential method
        """
        beauty = self.mcselect('Beauty',  BEAUTY)
        beauty_g = self. gselect('gBeauty', GBEAUTY)
        # if not beauty.empty() :
        #    return self.Warning('Skip event with beauty-production', SUCCESS )

        # get true MC-decays
        mcPsi = self.mcselect('mcPsi', 'J/psi(1S) => mu+ mu-')
        if mcPsi  .empty():
            return self.Warning('No True MC-decays are found', SUCCESS)

        goodMC = self.mcselect('goodMC',
                               mcPsi,
                               in_range(1.95, MCY, 4.55))

        if goodMC.empty():
            return self.Warning('No Good MC-decays are found', SUCCESS)

        if 1 != len(goodMC):
            self.Warning('Too many MC-decays are found %s' %
                         goodMC.size(), SUCCESS)

        # keep only the first ``good'' particle
        PSImc = goodMC[0]

        mcMuP_ = self.mcselect('mcMuP_', 'J/psi(1S) => ^mu+  mu-')
        if mcMuP_  .empty():
            return self.Warning('No True MC-decays (mu+) are found', SUCCESS)
        mcMuP = self.mcselect('mcMuP', mcMuP_,  FROMMCTREE(PSImc))
        if mcMuP   .empty():
            return self.Warning('No Good MC-decays (mu+) are found', SUCCESS)

        mcMuM_ = self.mcselect('mcMuM_', 'J/psi(1S) =>  mu+  ^mu-')
        if mcMuM_  .empty():
            return self.Warning('No True MC-decays (mu0) are found', SUCCESS)
        mcMuM = self.mcselect('mcMuM', mcMuM_,  FROMMCTREE(PSImc))
        if mcMuM   .empty():
            return self.Warning('No Good MC-decays (mu-) are found', SUCCESS)

        tup = self.nTuple('MCpsi')

        tup.column('hasBevt', not beauty   .empty(), 0, 1)
        tup.column('hasBevtG', not beauty_g .empty(), 0, 1)
        tup.column('hasB', hasB(PSImc), 0, 1)

        mup = mcMuP[0].momentum()
        mum = mcMuM[0].momentum()
        mc_lv01 = LoKi.Kinematics.decayAngle(mup, mup + mum)
        tup.column('mcpt', MCPT(PSImc) / GeV)
        tup.column('mcy', MCY(PSImc))
        tup.column('mclv01', mc_lv01)

        # create the MC-truth matching functors:
        trueMuP = MCTRUTH(self.mcTruth(), mcMuP)
        trueMuM = MCTRUTH(self.mcTruth(), mcMuM)
        truePsi = MCTRUTH(self.mcTruth(), PSImc)

        # get reconstructed muons
        mu_plus = self.select('Mu+', ('mu+' == ID) & trueMuP)
        mu_minus = self.select('Mu-', ('mu-' == ID) & trueMuM)

        psi = self.loop('Mu+ Mu-', 'J/psi(1S)')
        for j in psi:

            if not 2.9 * GeV < j.mass(1, 2) < 3.3 * GeV:
                continue

            if not 0 <= VCHI2(j) < 1000:
                continue

            if not truePsi(j):
                continue

            self.plot(M(j) / GeV, 'mass Jpsi->mu+mu-', 3.0, 3.2)

            j.save('Psi')

        jrec = self.selected('Psi')
        jsel = self.select('PSI', jrec, self._goodJpsi)
        self.treatRC(tup, jrec)

        tup.column('nJrec', jrec . size())

        tup.fArrayP('m_rc',
                    M / GeV,
                    'pt_rc',
                    PT / GeV,
                    'y_rc',
                    Y,
                    'lv01',
                    LV01,
                    jsel, 'nJsel', 100)

        tup.write()

        self.setFilterPassed(not jsel.empty())

        return SUCCESS


# ========================================================================
# @class MCDplus
#  Simple algorithm to get MC-efficiencies  for D+ ->K-pi+pi+ events
#  @date   2011-06-26
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
class MCDplus(MCD0):

    """
    Simple algorithm to get MC-efficiencies for D+ ->K-pi+pi+ events 
    """
    # the only one esential method:

    def analyse(self):
        """
        The only one essential method
        """
        beauty = self.mcselect('Beauty',  BEAUTY)
        beauty_g = self. gselect('gBeauty', GBEAUTY)
        # if not beauty.empty() :
        #    return self.Warning('Skip event with beauty-production', SUCCESS )

        # get true MC-decays
        mcDp = self.mcselect('mcDp', '[ D+ ==> K- pi+ pi+]CC')
        if mcDp .empty():
            return self.Warning('No True MC-decays are found', SUCCESS)

        goodMC = self.mcselect('goodMC',
                               mcDp,
                               in_range(1.95, MCY, 4.55) & (MCPT > 1.95 * GeV))

        if goodMC.empty():
            return self.Warning('No Good MC-decays are found', SUCCESS)

        if 1 != len(goodMC):
            self.Warning('Too many MC-decays are found %s' %
                         goodMC.size(), SUCCESS)

        # keep only the first ``good'' particle
        DPmc = goodMC[0]

        mcK_ = self.mcselect('mcK_', '[ D+ ==> ^K- pi+ pi+]CC')
        if mcK_   .empty():
            return self.Warning('No True MC-decays (K) are found', SUCCESS)
        mcK = self.mcselect('mcK', mcK_,  FROMMCTREE(DPmc))
        if mcK    .empty():
            return self.Warning('No Good MC-decays (K) are found', SUCCESS)
        #
        mcPi_ = self.mcselect('mcPi_', '[ D+ ==> K- ^pi+ ^pi+]CC')
        if mcPi_  .empty():
            return self.Warning('No True MC-decays (pi) are found', SUCCESS)
        mcPi = self.mcselect('mcPi', mcPi_, FROMMCTREE(DPmc))
        if mcPi   .empty():
            return self.Warning('No Good MC-decays (pi) are found', SUCCESS)

        tup = self.nTuple('MCDplus')

        tup.column('hasBevt', not beauty  .empty(), 0, 1)
        tup.column('hasBevtG', not beauty_g.empty(), 0, 1)
        tup.column('hasB', hasB(DPmc), 0, 1)

        tup.column('mcpt', MCPT(DPmc) / GeV)
        tup.column('mcy', MCY(DPmc))

        # create the MC-truth matching functors:
        truePi = MCTRUTH(self.mcTruth(), mcPi)
        trueK = MCTRUTH(self.mcTruth(), mcK)
        trueDP = MCTRUTH(self.mcTruth(), DPmc)

        # get reconstructed pions
        pions = self.select('pion', ('pi+' == ABSID) & truePi)
        # get reconstructed kaons
        kaons = self.select('kaon', ('K+' == ABSID) & trueK)

        dp = self.loop('kaon pion pion', 'D+')
        for d in dp:

            k = d(1)
            p1 = d(2)
            p2 = d(3)

            if Q(k) * Q(p1) > 0:
                continue
            if Q(k) * Q(p2) > 0:
                continue

            if not 1.8 * GeV < d.mass(1, 2, 3) < 2.0 * GeV:
                continue

            if not 0 <= VCHI2(d) < 1000:
                continue

            if not trueDP(d):
                continue

            if Q(k) < 0:
                d.setPID('D+')
            else:
                d.setPID('D-')

            self.plot(M(d) / GeV, 'mass D+ ->K-pi+pi+ ', 1.8, 1.95, 150)

            d.save('d+')

        drec = self.selected('d+')
        dsel = self.select('D+', drec, self._goodDp)

        tup.column('nDrec', drec . size())
        tup.column('nDsel', dsel . size())
        self.treatRC(tup, drec)

        tup.write()

        self.setFilterPassed(not dsel.empty())

        return SUCCESS


# ========================================================================
# @class MCDs
#  Simple algorithm to get MC-efficiencies  for Ds+ ->(K-K+)pi+ events
#  @date   2011-06-26
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
class MCDs(MCD0):

    """
    Simple algorithm to get MC-efficiencies for Ds+ ->K-K+pi+ events 
    """
    # the only one esential method:

    def analyse(self):
        """
        The only one essential method
        """
        beauty = self.mcselect('Beauty',  BEAUTY)
        beauty_g = self. gselect('gBeauty', GBEAUTY)
        # if not beauty.empty() :
        #    return self.Warning('Skip event with beauty-production', SUCCESS )

        # get true MC-decays
        mcDs = self.mcselect('mcDs', '[ D_s+ ==> K- K+ pi+]CC')
        if mcDs .empty():
            return self.Warning('No True MC-decays are found', SUCCESS)

        goodMC = self.mcselect('goodMC',
                               mcDs,
                               in_range(1.95, MCY, 4.55) & (MCPT > 1.95 * GeV))

        if goodMC.empty():
            return self.Warning('No Good MC-decays are found', SUCCESS)

        if 1 != len(goodMC):
            self.Warning('Too many MC-decays are found %s' %
                         goodMC.size(), SUCCESS)

        DSmc = []
        for d in goodMC:
            k1 = self.k_minus.child(d)
            if not k1:
                continue
            k2 = self.k_plus.child(d)
            if not k2:
                continue
            phi = k1.momentum() + k2.momentum()
            phiM = phi.M()
            if phiM > 1040 * MeV:
                continue
            #
            DSmc += [d]

        if not DSmc:
            return self.Warning(' phi is not found %s', SUCCESS)

        # keep only the first ``good'' particle
        DSmc = DSmc[0]

        mcK_ = self.mcselect('mcK_', '[ D_s+ ==> ^K- ^K+ pi+]CC')
        if mcK_   .empty():
            return self.Warning('No True MC-decays (K) are found', SUCCESS)
        mcK = self.mcselect('mcK', mcK_,  FROMMCTREE(DSmc))
        if mcK    .empty():
            return self.Warning('No Good MC-decays (K) are found', SUCCESS)
        #
        mcPi_ = self.mcselect('mcPi_', '[ D_s+ ==> K- K+ ^pi+]CC')
        if mcPi_  .empty():
            return self.Warning('No True MC-decays (pi) are found', SUCCESS)
        mcPi = self.mcselect('mcPi', mcPi_, FROMMCTREE(DSmc))
        if mcPi   .empty():
            return self.Warning('No Good MC-decays (pi) are found', SUCCESS)

        tup = self.nTuple('MCDs')

        tup.column('hasBevt', not beauty  .empty(), 0, 1)
        tup.column('hasBevtG', not beauty_g.empty(), 0, 1)
        tup.column('hasB', hasB(DSmc), 0, 1)

        tup.column('mcpt', MCPT(DSmc) / GeV)
        tup.column('mcy', MCY(DSmc))

        # create the MC-truth matching functors:
        truePi = MCTRUTH(self.mcTruth(), mcPi)
        trueK = MCTRUTH(self.mcTruth(), mcK)
        trueDS = MCTRUTH(self.mcTruth(), DSmc)

        # get reconstructed pions
        pions = self.select('pion', ('pi+' == ABSID) & truePi)
        # get reconstructed kaons
        kaons = self.select('kaon', ('K+' == ABSID) & trueK)

        dp = self.loop('kaon kaon pion', 'D_s+')
        for d in dp:

            k1 = d(1)
            k2 = d(2)
            pion = d(3)

            if Q(k1) * Q(k2) > 0:
                continue

            if not 1.8 * GeV < d.mass(1, 2, 3) < 2.0 * GeV:
                continue

            if not 0 <= VCHI2(d) < 1000:
                continue

            if not trueDS(d):
                continue

            if Q(pion) > 0:
                d.setPID('D_s+')
            else:
                d.setPID('D_s-')

            self.plot(M(d) / GeV, 'mass Ds+ ->K-K+pi+ ', 1.8, 1.95, 150)

            d.save('ds+')

        drec = self.selected('ds+')
        dsel = self.select('Ds+', drec, self._goodDs)
        self.treatRC(tup, drec)

        tup.column('nDsrec', drec . size())
        tup.column('nDssel', dsel . size())

        tup.write()

        self.setFilterPassed(not dsel.empty())

        return SUCCESS


# ========================================================================
# @class MCLc
#  Simple algorithm to get MC-efficiencies  for Lambda_c+ ->p K- pi+ events
#  @date   2011-06-26
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
class MCLc(MCD0):

    """
    Simple algorithm to get MC-efficiencies for D+ ->K-pi+pi+ events 
    """

    # the only one esential method:
    def analyse(self):
        """
        The only one essential method
        """

        beauty = self.mcselect('Beauty',  BEAUTY)
        beauty_g = self. gselect('gBeauty', GBEAUTY)
        # if not beauty.empty() :
        #    return self.Warning('Skip event with beauty-production', SUCCESS )

        # get true MC-decays
        mcLc = self.mcselect('mcLc', '[ Lambda_c+ ==> p+ K- pi+]CC')
        if mcLc .empty():
            return self.Warning('No True MC-decays are found', SUCCESS)

        goodMC = self.mcselect('goodMC',
                               mcLc,
                               in_range(1.95, MCY, 4.55) & (MCPT > 1.95 * GeV))

        if goodMC.empty():
            return self.Warning('No Good MC-decays are found', SUCCESS)

        if 1 != len(goodMC):
            self.Warning('Too many MC-decays are found %s' %
                         goodMC.size(), SUCCESS)

        # keep only the first ``good'' particle
        LCmc = goodMC[0]

        mcP_ = self.mcselect('mcP_', '[ Lambda_c+ ==> ^p+ K- pi+ ]CC')
        if mcP_   .empty():
            return self.Warning('No True MC-decays (p) are found', SUCCESS)
        mcP = self.mcselect('mcP', mcP_,  FROMMCTREE(LCmc))
        if mcP    .empty():
            return self.Warning('No Good MC-decays (p) are found', SUCCESS)
        #
        mcK_ = self.mcselect('mcK_', '[ Lambda_c+ ==>  p+ ^K- pi+ ]CC')
        if mcK_   .empty():
            return self.Warning('No True MC-decays (K) are found', SUCCESS)
        mcK = self.mcselect('mcK', mcK_,  FROMMCTREE(LCmc))
        if mcK    .empty():
            return self.Warning('No Good MC-decays (K) are found', SUCCESS)
        #
        mcPi_ = self.mcselect('mcPi_', '[ Lambda_c+ ==>  p+ K- ^pi+ ]CC')
        if mcPi_  .empty():
            return self.Warning('No True MC-decays (pi) are found', SUCCESS)
        mcPi = self.mcselect('mcPi', mcPi_,  FROMMCTREE(LCmc))
        if mcPi   .empty():
            return self.Warning('No Good MC-decays (pi) are found', SUCCESS)
        #

        tup = self.nTuple('MCLc')

        tup.column('hasBevt', not beauty   .empty(), 0, 1)
        tup.column('hasBevtG', not beauty_g .empty(), 0, 1)
        tup.column('hasB', hasB(LCmc), 0, 1)
        tup.column('bpid', bpid(LCmc))

        tup.column('mcpt', MCPT(LCmc) / GeV)
        tup.column('mcy', MCY(LCmc))

        # create the MC-truth matching functors:
        trueP = MCTRUTH(self.mcTruth(), mcP)
        trueK = MCTRUTH(self.mcTruth(), mcK)
        truePi = MCTRUTH(self.mcTruth(), mcPi)
        trueLC = MCTRUTH(self.mcTruth(), LCmc)

        # get reconstructed pions
        pions = self.select('pion', ('pi+' == ABSID) & truePi)
        # get reconstructed protons
        protons = self.select('proton', ('p+' == ABSID) & trueP)
        # get reconstructed kaons
        kaons = self.select('kaon', ('K+' == ABSID) & trueK)

        lc = self.loop('proton kaon pion', 'Lambda_c+')
        for l in lc:

            p = l(1)
            k = l(2)
            pion = l(3)

            if Q(k) * Q(p) > 0:
                continue
            if Q(k) * Q(pion) > 0:
                continue

            if not 2.2 * GeV < l.mass(1, 2, 3) < 2.4 * GeV:
                continue

            if not 0 <= VCHI2(l) < 1000:
                continue

            if not trueLC(l):
                continue

            if Q(p) > 0:
                l.setPID('Lambda_c+')
            else:
                l.setPID('Lambda_c~-')

            self.plot(M(l) / GeV, 'mass Lc+ -> pKpi ', 2.2, 2.35, 150)

            l.save('lc+')

        lrec = self.selected('lc+')
        lsel = self.select('Lc+', lrec, self._goodLc)
        self.treatRC(tup, lrec)

        tup.column('nLrec', lrec . size())
        tup.column('nLsel', lsel . size())

        tup.write()

        self.setFilterPassed(not lsel.empty())

        return SUCCESS


# =============================================================================
# common configuration part
def common_configure(datafiles, catalogs=[]):
    """
    Common configuration part
    """

    ## needed for job configuration
    from Configurables import DaVinci
    ## needed for job configuration
    from Configurables import EventSelector
    ## needed for job configuration
    from GaudiConf.Configuration import FileCatalog
    ## needed for job configuration
    from GaudiConf.Configuration import NTupleSvc

    from PhysConf.Filters import LoKi_Filters

    davinci = DaVinci(
        DataType='2010',
        ## InputType     = 'MDST' ,
        Simulation=True,
        PrintFreq=100,
        EvtMax=-1,
        #
        HistogramFile='MC_Histos.root',
        TupleFile='MC.root',
        #
        Lumi=True,
        #
        # OLD:
        ## DDDBtag   = "head-20101026" ,
        ## CondDBtag = "head-20101112"
        # NEW:
        ## DDDBtag   = "head-20110303" ,
        ## CondDBtag = "head-20110407"
    )

    from Configurables import CondDB
    CondDB().UseLatestTags = ["2010"]

    #
    # come back to Bender
    #

    setData(datafiles, catalogs)

    return SUCCESS

# =============================================================================
# configure the job


def configure_D0(datafiles, catalogs=[]):
    """
    Job configuration 
    """

    common_configure(datafiles, catalogs)

    #

    gaudi = appMgr()

    alg = MCD0(
        'MCD0',  # Algorithm name ,
        Inputs=['Phys/StdNoPIDsKaons/Particles',
                 'Phys/StdNoPIDsPions/Particles'],
        PP2MCs=['Relations/Rec/ProtoP/Charged']
    )
    #
    mainSeq = gaudi.algorithm('GaudiSequencer/DaVinciUserSequence', True)
    mainSeq.Members += [alg.name()]

    return SUCCESS

# =============================================================================
# configure the job


def configure_Psi(datafiles, catalogs=[]):
    """
    Job configuration 
    """

    common_configure(datafiles, catalogs)

    #

    gaudi = appMgr()

    alg = MCPsi(
        'MCPsi',  # Algorithm name ,
        Inputs=['Phys/StdLooseMuons/Particles'],
        PP2MCs=['Relations/Rec/ProtoP/Charged']
    )
    #
    mainSeq = gaudi.algorithm('GaudiSequencer/DaVinciUserSequence', True)
    mainSeq.Members += [alg.name()]

    return SUCCESS

# =============================================================================
# configure the job


def configure_Dplus(datafiles, catalogs=[]):
    """
    Job configuration 
    """

    common_configure(datafiles, catalogs)

    #

    gaudi = appMgr()

    alg = MCDplus(
        'MCplus',  # Algorithm name ,
        Inputs=['Phys/StdNoPIDsKaons/Particles',
                   'Phys/StdNoPIDsPions/Particles'],
        PP2MCs=['Relations/Rec/ProtoP/Charged']
    )
    #
    mainSeq = gaudi.algorithm('GaudiSequencer/DaVinciUserSequence', True)
    mainSeq.Members += [alg.name()]

    return SUCCESS

# =============================================================================
# configure the job


def configure_Ds(datafiles, catalogs=[]):
    """
    Job configuration 
    """

    common_configure(datafiles, catalogs)

    #

    gaudi = appMgr()

    alg = MCDs(
        'MCDs',  # Algorithm name ,
        Inputs=['Phys/StdNoPIDsKaons/Particles',
                 'Phys/StdNoPIDsPions/Particles'],
        PP2MCs=['Relations/Rec/ProtoP/Charged']
    )
    #
    mainSeq = gaudi.algorithm('GaudiSequencer/DaVinciUserSequence', True)
    mainSeq.Members += [alg.name()]

    return SUCCESS

# =============================================================================
# configure the job


def configure_Lc(datafiles, catalogs=[]):
    """
    Job configuration 
    """

    common_configure(datafiles, catalogs)

    #

    gaudi = appMgr()

    alg = MCLc(
        'MCLc',  # Algorithm name ,
        Inputs=['Phys/StdNoPIDsKaons/Particles',
                 'Phys/StdNoPIDsPions/Particles',
                 'Phys/StdNoPIDsProtons/Particles'],
        PP2MCs=['Relations/Rec/ProtoP/Charged']
    )
    #
    mainSeq = gaudi.algorithm('GaudiSequencer/DaVinciUserSequence', True)
    mainSeq.Members += [alg.name()]

    return SUCCESS

# =============================================================================
# configure the job


def configure(datafiles, catalogs=[]):
    """
    Job configuration 
    """
    # return configure_Psi ( datafiles , catalogs )
    # return configure_D0 ( datafiles , catalogs )
    # return configure_Dplus ( datafiles , catalogs )
    return configure_Ds(datafiles, catalogs)
    # return configure_Lc ( datafiles , catalogs )


# =============================================================================
# The actual job steering
if '__main__' == __name__:

    files = [
        # J/psi->mu+mu-
        #"/castor/cern.ch/grid/lhcb/MC/MC10/ALLSTREAMS.DST/00009117/0000/00009117_00000%03d_1.allstreams.dst" % i for i in range(1,720)
        # [ D*+ -> ( D0 => K- pi+ ) pi+ ]CC
        #"/castor/cern.ch/grid/lhcb/MC/MC10/ALLSTREAMS.DST/00008554/0000/00008554_00000%03d_1.allstreams.dst" % i for i in range(1,300)
        # D+ -> K- pi+ pi+  21263010
        #"/castor/cern.ch/grid/lhcb/MC/MC10/ALLSTREAMS.DST/00009404/0000/00009404_00000%03d_1.allstreams.dst" % i for i in range(1,400)
        # D_s+ -> K+ K- pi+   23263001
        #"/castor/cern.ch/grid/lhcb/MC/MC10/ALLSTREAMS.DST/00011124/0000/00011124_00000%03d_1.allstreams.dst" % i for i in range(1,600)
        # Lambda_c+ -> p K- pi+ 25103000
        "/castor/cern.ch/grid/lhcb/MC/MC10/ALLSTREAMS.DST/00011118/0000/00011118_00000%03d_1.allstreams.dst" % i for i in range(1, 255)
    ]

    configure(files)

    run(1000)

    gaudi = appMgr()


# =============================================================================
# The END
# =============================================================================
