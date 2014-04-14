#!/usr/bin/env python
# =============================================================================
# @file
#
# Simple script to get the generator cuts acceptance for D+
#
# The actual code has been kindly provided by Greig Cowan
#
# @author Vanya  BELYAEV Ivan.Belyaev@cern.ch
# @date   2011-06-18
#
#                   $Revision$
# Last modification $Date$
#                by $Author$
# =============================================================================
"""

Simple script to get the generator cuts acceptance for D+

The actual code has been kindly provided by Greig Cowan

"""
# =============================================================================
__author__ = " Vanya BELYAEV Ivan.Belyaev@cern.ch "
__date__ = " 201-06-18"
__version__ = " $Revision: 1.2 $ "
# =============================================================================
##from DiCharm.GenCuts.CheckGen import Dplus, configure_alg, run
from DiCharm.CheckGen import Dplus, configure_alg, run

# =============================================================================
# the only one important function


def configure(datafiles, catalogs=[]):
    """
    configure the job
    """
    return configure_alg(Dplus, datafiles, catalogs)

# =============================================================================
# job steering
if __name__ == '__main__':

    # make printout of the own documentations
    print '*' * 120
    print __doc__
    print ' Author  : %s ' % __author__
    print ' Version : %s ' % __version__
    print ' Date    : %s ' % __date__
    print '*' * 120

    configure([
        '/castor/cern.ch/grid/lhcb/user/i/ibelyaev/2011_06/21695/21695982/Gauss-21263010-1000ev-20110617.sim'
    ]
    )

    from Bender.MainMC import run
    run(1000)

# =============================================================================
# The END
# =============================================================================
