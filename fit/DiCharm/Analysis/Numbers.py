#!/usr/bin/env python
# =============================================================================
# @file DiCharm/Numbers.py
#
# Helper function to get certain numbers
#
# @author Vanya BELYAEV Ivan.Belyaev@cern.ch
# @date   2011-06-19
#
#                   $Revision$
# Last modification $Date$
#                by $Author$
# =============================================================================
"""
Helper function to get certain numbers 
"""
# =============================================================================
__version__ = "$Revision: 124897 $"
__author__ = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__ = "2011-06-07"
# =============================================================================
__all__ = (
    'cross_N',
    'cross_T',
)
# =============================================================================
import ROOT
from AnalysisPython.PyRoUts import hID
import DiCharm.Efficiency as Eff


# =============================================================================
# calculate the cross-section in nanobarns
#  @param N      efficiency corrected event yield
#  @param first  the first particle
#  @param second the second particle
#  @param the efficiency of global event cuts
#  @return the cross-section [in nb] and "irreducible systematics"
def cross_N(N,
            first,
            second,
            gec):
    """
    calculate the cross-section in nanobarns
    
    >>> N = ...
    >>> cross_section, syst = cross_N ( N , 'D0' , 'D0' )
    
    """

    br1 = Eff.br_ratio[first]
    br2 = Eff.br_ratio[second]

    factor = Eff.lumi_EXT * br1 * br2 * gec

    val1 = N.value() / factor / 1000  # nb
    val2 = N / factor.value() / 1000  # nb

    return val2, val1.error()

# =============================================================================
# calculate the cross-section in nanobarns
#  @param tree (or dataset) the tree
#  @param first  the first particle
#  @param second the second particle
#  @param the efficiency of global event cuts
#  @return the cross-section [in nb] and "irreducible systematics"


def cross_T(tree,
            first,
            second,
            gec):
    """
    calculate the cross-section in nanobarns
    
    >>> tree = ...
    >>> cross_section, syst = cross_T ( tree , 'D0' , 'D0' )
    
    """
    if issubclass(tree .__class__,  (ROOT.RooAbsData, ROOT.RooDataSet)):
        tree = tree.store().tree()

    h0 = ROOT.TH1F(hID(), '', 10, -1, 4)
    h0.Sumw2()
    tree.Project(h0.GetName(), '1', 'S1S2_sw')
    n_raw = h0.accumulate()
    tree.Project(h0.GetName(), '1', 'S1S2_sw*weight')
    n_corr = h0.accumulate()

    return n_raw, n_corr, cross_N(n_corr, first, second, gec)

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
