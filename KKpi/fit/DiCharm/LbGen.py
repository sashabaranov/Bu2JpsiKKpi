#!/usr/bin/env ipython
# ========================================================================
# @file DiCharm/LbGen.oy
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
from Bender.MainMC import *  # import all bender goodies
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


# =============================================================================
# @class LbGen
#  Simple template algorithm to study Lb
#  @date   2011-05-27
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
class LbGen (AlgoMC):

    """
    Simple algorithm to study Lb
    """

    # the only one esential method:
    def analyse(self):

        lb = self.gselect(
            'lb',
            '[ Lambda_b0 --> Lambda_c+ D_s- ... ]CC')

        if lb.empty():
            self.Warning(" No good candidates!")
            lb = self.gselect(
                'all', 'Lambda_b0' == GABSID)
            print lb
            return SUCCESS

        tup = self.nTuple('Lb')

        lc = LoKi.GenChild.Selector('Lambda_c+' == GABSID)
        ds = LoKi.GenChild.Selector('D_s+' == GABSID)

        for b in lb:

            tup.column('pt', GPT(b) / GeV)
            tup.column('p', GP(b) / GeV)

            p_lc = lc .child(b)
            p_ds = ds .child(b)

            m_lc = p_lc . momentum()
            m_ds = p_ds . momentum()

            mm = LoKi.Kinematics.mass(m_lc, m_ds) / GeV

            tup.column('mm', mm)

            tup.write()

        return SUCCESS

# =============================================================================
# configure the job


def configure(datafiles,
              catalogs=[],
              castor=False):
    """
    Job configuration 
    """

    from Configurables import DaVinci  # needed for job configuration

    the_year = '2012'

    davinci = DaVinci(
        DataType=the_year,
        InputType='DST',
        Simulation=True,
        PrintFreq=1000,
        EvtMax=-1,
        #
        HistogramFile='LcGen_Histos.root',
        TupleFile='LcGen-7.root',
        #
        Lumi=True,
        #
    )

    from Configurables import CondDB
    CondDB(LatestGlobalTagByDataType=the_year)

    # from Configurables import Gaudi__IODataManager as IODataManager
    # IODataManager().AgeLimit = 2

    # come back to Bender
    setData(datafiles, catalogs, castor)

    #
    # start Gaudi
    #
    gaudi = appMgr()

    alg = LbGen(
        'LbGen',  # Algorithm name ,
        PP2MCPs=[],
    )
    #

    gaudi.setAlgorithms([alg])

    return SUCCESS


# =============================================================================
# The actual job steering
if '__main__' == __name__:

    input = [
        '/afs/cern.ch/user/i/ibelyaev/cmtuser/Gauss_v46r2p1/Gauss-7-10000ev-20130702.gen',
    ]

    configure(input, castor=True)

    run(-1)

    gaudi = appMgr()

# =============================================================================
# The END
# =============================================================================
