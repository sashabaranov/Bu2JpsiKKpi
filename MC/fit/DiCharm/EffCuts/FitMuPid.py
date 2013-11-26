#!/usr/bin/env ipython
# =============================================================================
# $Id:$
# =============================================================================
# @file  DiCharm/EffCuts/FitMuPid.py
#
#  Calculate mu-PID efficiency for J/psi
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
#  @date   2011-07-04
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#
#                    $Revision$
#  Last modification $Date$
#                 by $Author$
# =============================================================================
"""

Calculate mu-PID efficiency for J/psi 


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
__date__ = "2011-07-04"
__version__ = '$Revision$'
# =============================================================================

import ROOT
from AnalysisPython.PyRoUts import VE


def fitMuPid(func,
             fitfunc,
             key,
             dbase='$DICHARMROOT/data/DiCharm.db'):

    import AnalysisPython.ZipShelve as ZipShelve
    db = ZipShelve.open(dbase, 'r')
    histos = db[key]
    db.close()

    pids = {}
    keys = histos[0].keys()
    keys.sort()
    keys.reverse()

    nTotal = VE(0, 0)
    nCorr = VE(0, 0)

    for key in keys:

        hists = histos[0][key]

        nums = []
        m0_ = VE(-1000, 0)
        s0_ = VE(-1000, 0)
        b0_ = VE(-1000, 0)

        for h in hists:
            if m0_.value() > 0 and s0_.value() > 0:
                n, r = func(h,
                            fitfunc,
                            fixMass=m0_.value(),
                            fixSigma=s0_.value())
            else:
                n, r = func(h, fitfunc)
                if hasattr(r, 'Status') and 0 == r.Status():
                    m0_ = r[1]
                    s0_ = r[2]
                    b0_ = r[3]

            if not hasattr(r, 'Status'):
                print 'WARNING(1): ', key, r
            elif 0 != r.Status():
                print 'WARNING(2): ', key, r

            nums += [n]

        if abs(nums[2].value()) < 0.5:
            print 'WARNING(3): adjust ', nums[2]
            nums[2] = VE(1.1, 1)

        ePID = nums[1].frac(nums[2])

        nTotal += nums[0]

        pids[key] = ePID, nums

        print ' key, nums : ', key, ePID, nums, m0_, s0_, b0_

    nTT = VE(0, 0)
    nTA = VE(0, 0)
    nTR = VE(0, 0)
    for k in pids:
        e = pids[k]
        nTT += e[1][0]
        nTA += e[1][1]
        nTR += e[1][2]

    epid = nTA.frac(nTR)

    print ' #total ', nTT, nTA, nTR, epid

    return pids, histos[1], histos[2]

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
