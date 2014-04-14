#!/usr/bin/env ipython
# ==========================================================================================#
# $Id:$
# ========================================================================
# @file  DiCharm/Analysis/Slices.py
#
#  Helper structure for fits in slices
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
#  @date   2011-11-14
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#
#                    $Revision$
#  Last modification $Date$
#                 by $Author$
# =============================================================================
"""

Helper structure for fits in slices 

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
__date__ = "2011-11-14"
__version__ = '$Revision$'
# =============================================================================
import ROOT
from AnalysisPython.PyRoUts import hID, SE
# =============================================================================
# get the average binwidth for axis


def binw(axis):
    """
    Get binwidth 
    """
    return float(axis.GetXmax() - axis.GetXmin()) / axis.GetNbins()

ROOT.TAxis.binw = binw
ROOT.TH1D .binw = lambda s: s.GetXaxis().binw()
ROOT.TH1F .binw = lambda s: s.GetXaxis().binw()


# =============================================================================
# @class Slices
#  helper structure to keep the slices
#  @date   2011-11-14
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
class Slices(object):

    """
    Helper structure for fits in slices 
    """
    # constructor fomr two historgam templates

    def __init__(self, h1, h2):

        h = h1.Clone(hID())
        h.Reset()
        if not h.GetSumw2():
            h.Sumw2()
        self.h1_ = h

        h = h2.Clone(hID())
        h.Reset()
        if not h.GetSumw2():
            h.Sumw2()
        self.h2_ = h

        self.axis_1 = self.h1_.GetXaxis()
        self.axis_2 = self.h2_.GetXaxis()

        self.histos_1 = {}
        self.histos_2 = {}

        for ibin in self.axis_2:
            h = self.h1_.Clone(hID())
            self.histos_1[ibin] = h

        for ibin in self.axis_1:
            h = self.h2_.Clone(hID())
            self.histos_2[ibin] = h

    # fill the structure
    def fill(self, m1, m2, weight=1):

        i1 = self.axis_1.FindBin(m1)
        if not i1 in self.axis_1:
            return False

        i2 = self.axis_2.FindBin(m2)
        if not i2 in self.axis_2:
            return False

        self.h1_             . Fill(m1, weight)
        self.h2_             . Fill(m2, weight)

        self.histos_1[i2] . Fill(m1, weight)
        self.histos_2[i1] . Fill(m2, weight)

        return True

    # fit the slices:
    def fit(self,
            fun1,
            fun2):

        h1 = self.h1_.Clone(hID())
        h2 = self.h2_.Clone(hID())

        mnerr = 1.e+9

        n1, r1 = fun1(self.h1_)
        for k in self.histos_1:
            h = self.histos_1[k]
            n, r = fun1(h)
            if r and 0 == r.Status() and 0 != n.value() and 0 < n.error():
                mnerr = min(mnerr, n.error())
            h2[k] = abs(n)

        n2, r2 = fun2(self.h2_)
        for k in self.histos_2:
            h = self.histos_2[k]
            n, r = fun2(h)
            if r and 0 == r.Status() and 0 != n.value() and 0 < n.error():
                mnerr = min(mnerr, n.error())
            h1[k] = abs(n)

        # adjust null entries (== empty histograms)
        for h in (h1, h2):
            for i in h:
                n = h[i]
                if n.error() < mnerr:
                    n.setError(mnerr)
                    h[i] = n

        n1, r1 = fun1(h1)
        n2, r2 = fun2(h2)

        return (n1, r1, h1), (n2, r2, h2)


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
