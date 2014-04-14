#!/usr/bin/env python
# =============================================================================
# $Id:$
# =============================================================================
# @file DiCharm/Utils.py
# Set of small techical utilities for 2xCharm analysis
#
# @author Vanya BELYAEV Ivan.Belyaev Ivan.BElyaev@cern.ch
# @date   2011-07-05
#
#                   $Revision$
# Last Modification $Date$
#                by $Author$
# =============================================================================
"""
Set of small techical utilities for 2xCharm analysis
"""
# =============================================================================
__version__ = "$Revision:$"
__author__ = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__ = "2011-06-07"
# =============================================================================
__all__ = (
    'lfn2disk',
    'LHCbStyle',
    'getLumi'
)
# =============================================================================
import ROOT
# =============================================================================
# load LHCb-style file


def LHCbStyle(stylefile='$HOME/public/lhcbStyle.C'):
    """
    Load LHCb-Style fiel for plots 
    """

    if ROOT.gROOT.IsBatch():
        ROOT.gROOT.SetBatch(False)

    import os
    if not os.path.exists(stylefile):
        stylefile = os.path.expandvars(stylefile)
    if not os.path.exists(stylefile):
        stylefile = os.path.expanduser(stylefile)
    if not os.path.exists(stylefile):
        stylefile = '/afs/cern.ch/user/i/ibelyaev/public/lhcbStyle.C'

    if os.path.exists(stylefile):
        ROOT.gROOT.LoadMacro(stylefile)
        ROOT.LHCb.setLHCbStyle()
    else:
        print "WARNING: can't find/load LHCb-style file"

    #
    # put some minor polishing atop of the style
    #
    # ROOT.gStyle.SetOptStat (0)

    ROOT.gStyle.SetOptFit(0)
    ROOT.gStyle.SetOptFit(1111)
    # ROOT.gStyle.SetOptFit  (0)

    if ROOT.gROOT.IsBatch():
        ROOT.gROOT.SetBatch(False)

    return ROOT.gROOT.IsBatch()

# =============================================================================
# copy castor files into some DISK location


def lfn2disk(lst,
             lfndir='/afs/cern.ch/project/lbcern/vol7/ibelyaev/DATA/LFN'):
    """
    Copy castor file (by LFN) into disk  and return the list of
    names of disk
    """

    import os
    if not lfndir:
        lfndir = '/tmp'
        print " DIRECTORY is not specified , use '%s'" % lfndir

    if isinstance(lst, str):
        lst = [lst]

    result = []
    size = 0
    import os
    import tempfile
    for l in lst:
        #
        # already real file name ?
        #
        if os.path.exists(l):
            result += [l]
            size += os.path.getsize(l)
            continue
        #
        # valid directory?
        #
        if not os.path.exists(lfndir) or not os.path.isdir(lfndir):
            lfndir = os.path.expandvars(lfndir)
            lfndir = os.path.expandvars(lfndir)
            lfndir = os.path.expandvars(lfndir)
            lfndir = os.path.abspath(lfndir)
        #
        if not os.path.exists(lfndir) or not os.path.isdir(lfndir):
            # try to create it:
            print ' Try to create LFNDIR: ', lfndir
            os.makedirs(lfndir)
        #
        if not os.path.exists(lfndir) or not os.path.isdir(lfndir):
            raise TypeError, 'invalid "LFNDIR"'
        #
        #
        f = lfndir + l
        #
        # already copied to disk ?
        #
        if os.path.exists(f):
            result += [f]
            size += os.path.getsize(f)
            continue
        tmp = lfndir + '/_tmp'
        if os.path.exists(tmp):
            os.remove(tmp)
        #
        # copy from castor:
        #
        cmd = 'rfcp %s %s ' % ('/castor/cern.ch/grid' + l, tmp)
        cin, cout, cerr = os.popen3(cmd)
        for a in cerr:
            print a
        # for a in cout : print a
        if os.path.exists(tmp):
            os.renames(tmp, f)
        if os.path.exists(f):
            result += [f]
            size += os.path.getsize(f)
            print 'Disk size: ', l, " %.1fkB " % int(os.path.getsize(f) / 1000.0), len(result), len(lst)
        if os.path.exists(tmp):
            os.remove(tmp)
    if result:
        print 'DIRECTORY: ', lfndir, ' #FILES:', len(result), " SIZE: %.1fkB" % int(size / 1000.0)
    return result

# =============================================================================
# tiny decoration for shelve module
# =============================================================================
import shelve
_old_shelve_open_ = shelve.open


def _new_shelve_open_(filename, *kargs, **kwargs):
    """
    A bit extended version of open
    """
    import os
    filename = os.path.expandvars(filename)
    filename = os.path.expandvars(filename)
    filename = os.path.expandvars(filename)
    filename = os.path.abspath(filename)
    #
    return _old_shelve_open_(filename, *kargs, **kwargs)

_new_shelve_open_ .__doc__ += '\n' + _old_shelve_open_ .__doc__

shelve.open = _new_shelve_open_

# =============================================================================
# Real decoration
# =============================================================================
# get the lumi


def getLumi(tree):
    """
    Get Lumi 
    """

    import ROOT
    from AnalysisPython.PyRoUts import hID, VE
    #
    hL = ROOT.TH1F(hID(), '', 10, -1, 4)
    hL.Sumw2()
    tree.Project(hL.GetName(), '1', 'IntegratedLuminosity')
    l0 = hL.accumulate()
    #
    tree.Project(hL.GetName(), '1',
                 'IntegratedLuminosity+IntegratedLuminosityErr')
    l1 = hL.accumulate()
    #
    # tree.Project ( hL.GetName() , 'IntegratedLuminosity+IntegratedLuminosityErr')
    # l2 = hL.accumulate()
    #
    return VE(l0.value(), (l1.value() - l0.value()) ** 2)


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
