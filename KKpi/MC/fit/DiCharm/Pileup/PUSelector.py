#!/usr/bin/env ipython
# ==========================================================================================#
# $Id:$
# ========================================================================
# @file  DiCharm/Pileup/Selectors.py
#
#  Various selectors for 2xCharm analysis
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
#  @date   2011-07-22
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#
#                    $Revision$
#  Last modification $Date$
#                 by $Author$
# =============================================================================
"""

Various selectors for 2xCharm analysis

This file is a part of BENDER project:
    ``Python-based Interactive Environment for Smart and Friendly Physics Analysis''

The project has been designed with the kind help from Pere MATO and Andrey TSAREGORODTSEV. 

And it is based on the LoKi project:
    ``C++ ToolKit for Smart and Friendly Physics Analysis''

By usage of this code one clearly states the disagreement 
with the smear campain of Dr.O.Callot et al.: 
    ``No Vanya's lines are allowed in LHCb/Gaudi software.''

"""
# =============================================================================
__author__ = 'Vanya BELYAEV  Ivan.Belyaev@cern.ch'
__date__ = "2011-09-10"
__version__ = '$Revision$'
# =============================================================================
import ROOT
from AnalysisPython.PyRoUts import hID, VE, SE
from AnalysisPython.PySelector import Selector
# =============================================================================
from math import atan2, pi
# =============================================================================

import DiCharm.Efficiency as Eff

# ========================================================================
# Reweight generator-level monte carlo
#  @date   2011-07-22
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class PUSel(Selector):

    """
    Reweight generator-level monte carlo 
    """

    def __init__(self,
                 scale=VE(1, 0)):

        Selector.__init__(self, None, self)  # initialize the base

        self.scale = scale
        self.dz_0 = ROOT.TH1F(hID(), '', 500, -250, 250)
        self.dz_1 = ROOT.TH1F(hID(), '', 100,  -10,  10)
        self.pid_0 = ROOT.TH1F(hID(), '', 204,  -51,  51)  # true pileup
        ## pileup with |dz|
        self.pid_1 = ROOT.TH1F(hID(), '', 204,  -51,  51)
        self.pid_2 = ROOT.TH1F(hID(), '', 204,  -51,  51)  # B-feeddown
        self.pid_3 = ROOT.TH1F(hID(), '', 204,  -51,  51)  # signal

        self.dz_0  . Sumw2()
        self.dz_1  . Sumw2()
        self.pid_0 . Sumw2()
        self.pid_1 . Sumw2()
        self.pid_2 . Sumw2()
        self.pid_3 . Sumw2()

        self.se_0 = SE()
        self.se_1 = SE()
        self.se_2 = SE()
        self.se_3 = SE()

        self._events = 0

    def eff(self, pid, pt, y, pt_min=2):

        if not 2 <= y < 4:
            return VE(0, 0)
        if 12 <= pt:
            return VE(0, 0)

        if 443 == pid:
            return Eff.selEff_Jpsi_2D(pt, y)
        elif 421 == abs(pid):
            if pt < pt_min:
                return VE(0, 0)
            return Eff.selEff_D0(pt, y)
        elif 411 == abs(pid):
            if pt < pt_min:
                return VE(0, 0)
            return Eff.selEff_Dp(pt, y)
        elif 431 == abs(pid):
            if pt < pt_min:
                return VE(0, 0)
            return Eff.selEff_Ds(pt, y)
        elif 4000 < abs(pid):
            if pt < pt_min:
                return VE(0, 0)
            return Eff.selEff_Lc(pt, y)
        #
        return VE(0.0, 0.0)

    def trg(self, pid, pt, y):

        if 443 == pid:
            return Eff.selEff_Jpsi_2D(pt, y)
        elif 421 == abs(pid):
            return Eff.selEff_D0(pt, y)
        elif 411 == abs(pid):
            return Eff.selEff_Dp(pt, y)
        elif 431 == abs(pid):
            return Eff.selEff_Ds(pt, y)
        elif 4000 < abs(pid):
            return Eff.selEff_Lc(pt, y)
        #
        return VE(0, 0)

    def br(self, pid):

        if 443 == pid:
            return Eff.br__Jpsi
        elif 421 == abs(pid):
            return Eff.br_D0
        elif 411 == abs(pid):
            return Eff.br_Dp
        elif 431 == abs(pid):
            return Eff.br_Ds
        elif 4000 < abs(pid):
            return Eff.br_Lc
        #
        return VE(0.0, 0.0)

    def same(self, bamboo, i, j):

        if abs(bamboo.s_px[i] - bamboo.c_px[j]) < 0.001:
            return True
        if abs(bamboo.s_py[i] - bamboo.c_py[j]) < 0.001:
            return True
        if abs(bamboo.s_pz[i] - bamboo.c_pz[j]) < 0.001:
            return True
        if abs(bamboo.s_pt[i] - bamboo.c_pt[j]) < 0.001:
            return True
        if abs(bamboo.s_y[i] - bamboo.c_y[j]) < 0.001:
            return True
        #
        return False

    # the only one important method
    def Process(self, entry):
        """
        Fills data set 
        """
        #
        # == getting the next entry from the tree
        #
        if self.GetEntry(entry) <= 0:
            return 0  # RETURN

        self._events += 1
        if 0 == self._events % 100000:
            print self._events, self.dz_0.Integral(), self.pid_1.Integral() + self.pid_2.Integral()
        #
        # == for more convenience
        #
        bamboo = self.fChain

        nS = bamboo.nSignal
        nC = bamboo.nC

        for i in range(0, nS):

            pid_s = bamboo.s_id[i]
            pt_s = bamboo.s_pt[i] / 1000
            y_s = bamboo.s_y[i]
            b_s = bamboo.s_fromB[i]
            pvk_s = bamboo.s_pvk[i]

            eff_s = self.eff(pid_s, pt_s, y_s)
            if 0 >= eff_s.value():
                continue
            trg_s = self.trg(pid_s, pt_s, y_s)
            if 0 >= trg_s.value():
                continue

            for j in range(0, nC):

                pid_c = bamboo.c_id[j]
                pt_c = bamboo.c_pt[j] / 1000
                y_c = bamboo.c_y[j]
                b_c = bamboo.c_fromB[j]
                pvk_c = bamboo.c_pvk[j]

                if pid_s == pid_c and self.same(bamboo, i, j):
                    continue

                eff_c = self.eff(pid_c, pt_c, y_c)
                eff_c *= self.br(pid_c)
                if 0 >= eff_c.value():
                    continue

                trg_c = self.trg(pid_c, pt_c, y_c)
                if 0 >= trg_c.value():
                    continue

                # total efficiency
                eff_0 = eff_s * eff_c
                trg_0 = VE(1, 0)

                # true pileup (different PVs )
                pileup = pvk_s != pvk_c
                fromB = b_s or b_c

                if 443 == pid_s:
                    trg_0 = trg_s
                else:
                    trg_0 = 1 - (1 - trg_s) * (1 - trg_c)

                # correct for lower B-efficiencies
                if b_s and 443 == pid_s:
                    eff_0 /= 10
                elif b_s:
                    eff_0 /= 2 * 2

                if b_c:
                    eff_0 /= 2

                weight = eff_0 * trg_0 * self.scale
                w = weight.value()

                dz = bamboo.s_pvz[i] - bamboo.c_pvz[j]
                dz = min(dz,  249)
                dz = max(dz, -249)

                pid = abs(pid_c)
                pid = min(pid, 450)
                pid -= 400

                if 443 != pid_s and pid_c * pid_s < 0:
                    pid *= -1

                if pileup:
                    self.dz_0.Fill(dz, w)
                    self.pid_0.Fill(pid, w)
                    self.se_0 += w
                    if abs(dz) < 1.5:
                        self.dz_1.Fill(dz, w)
                        self.pid_1.Fill(pid, w)
                        self.se_1 += w
                elif fromB:
                    # print 'fromB:' , pid
                    if 443 == pid_s:
                        w /= (300. / 5)
                    self.pid_2.Fill(pid, w)
                    self.se_2 += w
                else:
                    # print 'true:' , pid
                    self.pid_3.Fill(pid, w)
                    self.se_3 += w

        return 1


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
