#!/usr/bin/env python
# =============================================================================
# $Id:$
# =============================================================================
# @file DiCharm/GenCuts/GenDs.py
#
# Simple script to get the generator acceptance for Ds -> phi pi+
# Event type: 23103002
#
# The actual code has been kindly provided by Greig COWAN
#
# The algorithm gets the generator-files (without cuts)
# and (re)applies the GenCut-tool
#
# @attention one needs to build package Gen/GenCuts!
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
# @author Vanya  BELYAEV Ivan.Belyaev@cern.ch
# @thanks Greig  COWAN
# @date   2013-09-03
#
#                   $Revision: 155091 $
# Last modification $Date: 2013-04-16 19:47:24 +0200 (Tue, 16 Apr 2013) $
#                by $Author: ibelyaev $
# =============================================================================
"""

Simple script to get the generator acceptance for Ds -> phi pi+
Event type: 23103002

The actual code has been kindly provided by Greig COWAN

The algorithm gets the generator-files (without cuts)
and (re)applies the GenCut-tool

ATTENTION: one needs to build package Gen/GenCuts! 

This file is a part of BENDER project:
``Python-based Interactive Environment for Smart and Friendly Physics Analysis''

The project has been designed with the kind help from
Pere MATO and Andrey TSAREGORODTSEV. 

And it is based on the 
LoKi project: ``C++ ToolKit for Smart and Friendly Physics Analysis''

By usage of this code one clearly states the disagreement 
with the smear campaign of Dr.O.Callot et al.: 
``No Vanya's lines are allowed in LHCb/Gaudi software.''


"""
# =============================================================================
__author__ = " Vanya BELYAEV Ivan.Belyaev@cern.ch "
__date__ = " 2013-09-03"
__version__ = " $Revision: 155091 $ "
# =============================================================================
# import everything from bender
from Bender.MainMC import *
from GaudiKernel.SystemOfUnits import GeV
from GaudiKernel.PhysicalConstants import c_light
# =============================================================================
# logging
# =============================================================================
from Bender.Logger import getLogger
logger = getLogger(__name__)
# =============================================================================

# =============================================================================
# configure it


def configure(datafiles,
              catalogs=[],
              castor=False):

    #
    #
    #
    from Configurables import LoKi__GenCutTool
    tightCut = LoKi__GenCutTool('TightCut')
    tightCut.PropertiesPrint = True
    tightCut.Decay = '^[D_s+ ==> ^K- ^K+ ^pi+]CC'
    tightCut.Preambulo += [
        'from GaudiKernel.SystemOfUnits import millimeter, micrometer,MeV,GeV',
        'GY           =  LoKi.GenParticles.Rapidity () ## to be sure       ',
        'inAcc        =  in_range ( 0.005 , GTHETA , 0.400 )               ',
        'inEta        =  in_range ( 1.95  , GETA   , 5.050 )               ',
        'fastTrack    =  ( GPT > 220 * MeV ) & ( GP  > 3.0 * GeV )         ',
        'goodTrack    =  inAcc & inEta & fastTrack                         ',
        'longLived    =  75 * micrometer < GTIME                           ',
        'inY          =  in_range ( 1.9   , GY     , 4.6   )               ',
        'goodDs       =  inY & longLived     & ( GPT > 0.9 * GeV )         ',
        'Bancestors   =  GNINTREE ( GBEAUTY , HepMC.ancestors )            ',
        'notFromB     =  0 == Bancestors                                   '
    ]
    tightCut.Cuts = {
        '[D_s+]cc': 'goodDs & notFromB ',
        '[K+]cc': 'goodTrack ',
        '[pi+]cc': 'goodTrack '
    }

    from BenderTools.GenCuts import configure_decay
    return configure_decay(
        'Alg',
        "( [D_s+ ==> ( phi(1020) => K- K+ ) pi+]CC ) && ~( Charm  --> pi0 ...  ) ",
        datafiles,
        catalogs,
        castor,
        tool='LoKi::GenCutTool/TightCut:PUBLIC')


# ========================================================================
#
if '__main__' == __name__:

    input = [
        '/lhcb/user/i/ibelyaev/2013_09/57296/57296391/Gauss-23103002-20ev-20130901.gen',
        '/lhcb/user/i/ibelyaev/2013_09/57296/57296394/Gauss-23103002-20ev-20130901.gen',
        '/lhcb/user/i/ibelyaev/2013_09/57296/57296397/Gauss-23103002-20ev-20130901.gen',
        '/lhcb/user/i/ibelyaev/2013_09/57296/57296404/Gauss-23103002-20ev-20130901.gen',
        '/lhcb/user/i/ibelyaev/2013_09/57296/57296414/Gauss-23103002-20ev-20130901.gen',
        '/lhcb/user/i/ibelyaev/2013_09/57296/57296424/Gauss-23103002-20ev-20130901.gen',
        '/lhcb/user/i/ibelyaev/2013_09/57296/57296432/Gauss-23103002-20ev-20130901.gen',
        '/lhcb/user/i/ibelyaev/2013_09/57296/57296439/Gauss-23103002-20ev-20130901.gen',
        '/lhcb/user/i/ibelyaev/2013_09/57296/57296444/Gauss-23103002-20ev-20130901.gen',
        '/lhcb/user/i/ibelyaev/2013_09/57296/57296455/Gauss-23103002-20ev-20130901.gen'
    ]

    configure(input, catalogs=[], castor=True)

    run(10)

# ========================================================================
# The END
# ========================================================================
