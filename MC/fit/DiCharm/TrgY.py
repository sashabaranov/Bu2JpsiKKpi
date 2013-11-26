#!/usr/bin/env python
# ========================================================================
# @file DiCharm/TrgY.py
#
#  Get trigger efficiency for Y
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
#  @date   2011-08-18
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#
#                    $Revision$
#  Last modification $Date$
#                 by $Author$
# =============================================================================
"""

Get trigger efficiency for Y

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
__date__ = "2011-07-18"
__version__ = '$Revision$'
# =============================================================================
## needed to produce/visualize the histograms
import ROOT
from Bender.Main import *  # import all bender goodies
## easy access to various geometry routines
import LHCbMath.Types
from Gaudi.Configuration import *  # needed for job configuration
# =============================================================================
from GaudiKernel.SystemOfUnits import GeV, MeV, mm
from GaudiKernel.PhysicalConstants import c_light
# =============================================================================
# logging
# =============================================================================
from Bender.Logger import getLogger
logger = getLogger(__name__)
# =============================================================================
import BenderTools.TisTos  # add methods for TisTos
import BenderTools.Fill  # add methods to fill n-tuple info
# =============================================================================

# =============================================================================
# @class TrgY
#  Simple algorithm to studty trigger efficiency for Y->mumu
#  @date   2011-05-27
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class TrgY (Algo):

    """
    Simple algorithm 
    """

    def initialize(self):

        sc = Algo.       initialize(self)
        if sc.isFailure():
            return sc

        triggers = {}
        from DiCharm.TisTos import lines

        sc = self.tisTos_initialize(triggers, lines)
        if sc.isFailure():
            return sc

        sc = self.  fill_initialize()
        if sc.isFailure():
            return sc

        self.mfit = DTF_FUN(M, True)
        self.c2dtf = DTF_CHI2NDOF(True)
        self.ctau = DTF_CTAU(0, True)

        return SUCCESS

    # finalize the algorithm
    def finalize(self):
        """
        Finalize the action 
        """
        self.mfit = None
        self.c2dtf = None
        self.ctau = None

        self.  fill_finalize()
        self.tisTos_finalize()

        return Algo.finalize(self)

    # the only one esential method:
    def analyse(self):
        """
        The only one essential method
        """

        primaries = self.vselect('PVs', ISPRIMARY)

        if primaries.empty():
            return self.Warning('No primary vertices are found', SUCCESS)

        #
        # rec summary
        rc_summary = self.get('/Event/Rec/Summary').summaryData()
        #
        ups = self.select('ups', ' Meson -> mu+ mu-')

        tup = self.nTuple('Y')
        for y in ups:
            #
            m = M(y) / GeV
            if not 8 < m < 15:
                continue
            #
            mu1 = y(1)
            mu2 = y(2)

            self.treatKine(tup, y, '_Y')

            tup.column_float('mass', self.mfit(y) / GeV)
            tup.column_float('c2dtf', self.c2dtf(y))

            # add the information needed for TisTos
            self.tisTos(y, tup, 'Y_',
                        self.lines['Upsilon'],
                        self.l0tistos,
                        self.tistos)

            # add the information needed for TisTos
            self.tisTos(y, tup, 'Y1_',
                        self.lines['Upsilon1'],
                        self.l0tistos,
                        self.tistos)

            # add the information for Pid efficiency correction
            self.treatMuons(tup, y)
            # add the infomation needed for track efficiency correction
            self.treatTracks(tup, y)

            # add some reco-summary information
            self.addRecSummary(tup, rc_summary)
            #
            tup.write()

        return SUCCESS

# =============================================================================
# configure the job


def configure(datafiles,
              catalogs=[],
              castor=False,
              params={}):
    """
    Job configuration 
    """

    ## needed for job configuration
    from Configurables import DaVinci

    the_year = "2012"

    from BenderTools.Parser import hasInFile

    if params:
        the_year = params['Year']
        logger.info('Year is set from params to be %s ' % the_year)
    else:
        if hasInFile(datafiles, 'Collision11'):
            the_year = '2011'
        elif hasInFile(datafiles, 'Collision12'):
            the_year = '2012'
        elif hasInFile(datafiles, 'Collision13'):
            the_year = '2013'
        logger.info('Year is set from files  to be %s ' % the_year)

    #
    # check
    #
    if '2011' == the_year and hasInFile(datafiles, 'Collision12'):
        raise AttributeError, 'Invalid Year %s ' % the_year
    if '2012' == the_year and hasInFile(datafiles, 'Collision11'):
        raise AttributeError, 'Invalid Year %s ' % the_year

    logger.info('Use the Year = %s ' % the_year)

    rootInTES = '/Event/BOTTOM'

    dv = DaVinci(
        DataType=the_year,
        InputType='DST',
        RootInTES=rootInTES,
        Simulation=False,
        PrintFreq=1000,
        EvtMax=-1,
        #
        HistogramFile='Y_Histos.root',
        TupleFile='Y.root',
        #
        Lumi=True,
        #
    )

    from Configurables import CondDB
    CondDB(LatestGlobalTagByDataType=the_year)

    dv.UserAlgorithms = ['Y']

    from BenderTools.Utils import silence
    silence()

    #
    # come back to Bender
    #
    setData(datafiles, catalogs, castor)

    #
    # start Gaudi
    #
    gaudi = appMgr()

    alg = TrgY(
        'Y',
        RootInTES=rootInTES,
        Inputs=['Phys/Upsilon/Particles']
    )

    return SUCCESS

# =============================================================================
# The actual job steering
if '__main__' == __name__:

    input = [
        '/lhcb/LHCb/Collision12/BOTTOM.MDST/00024670/0000/00024670_000000%02d_1.bottom.mdst' % i for i in range(1, 100)
    ]

    configure(input, castor=True)

    run(5000)

    gaudi = appMgr()

# =============================================================================
# The END
# =============================================================================
