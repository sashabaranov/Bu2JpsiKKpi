#!/usr/bin/env ipython
# ========================================================================
# @file DiCharm/Alg
#
#  Simple template algorithm for Charm+X
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

Simple template algorithm for 2xCharm

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
# @class BcGen
#  Simple template algorithm to study Bc
#  @date   2011-05-27
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
class BcGen (AlgoMC):

    """
    Simple algorithm to study Charm & Z production
    """

    # the only one esential method:
    def analyse(self):

        bc = self.gselect(
            'bc',
            '[ B_c+ ==> ( J/psi(1S) => mu+ mu- ) ( D_s+ ==> K+ K- pi+ ) ]CC')

        if bc.empty():
            self.Warning(" No good candidates!")
            bc = self.gselect(
                'bcc', 'B_c+' == GABSID)
            print bc
            return SUCCESS

        tup = self.nTuple('Bc')

        psi = LoKi.GenChild.Selector(
            '[ B_c+ ==> ^( J/psi(1S) => mu+ mu- )  ( D_s+ ==> K+ K- pi+ ) ]CC')
        ds = LoKi.GenChild.Selector(
            '[ B_c+ ==>  ( J/psi(1S) => mu+ mu- ) ^( D_s+ ==> K+ K- pi+ ) ]CC')

        for b in bc:

            tup.column('pt', GPT(b) / GeV)
            tup.column('p', GP(b) / GeV)

            p_psi = psi.child(b)
            p_ds = ds .child(b)

            m_psi = p_psi. momentum()
            m_ds = p_ds . momentum()

            mm = LoKi.Kinematics.mass(m_psi, m_ds) / GeV

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

    the_year = '2011'

    davinci = DaVinci(
        DataType=the_year,
        InputType='DST',
        Simulation=True,
        PrintFreq=1000,
        EvtMax=-1,
        #
        HistogramFile='BcGen_Histos.root',
        TupleFile='BcGen.root',
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

    alg = BcGen(
        'BcGen',  # Algorithm name ,
    )
    #

    gaudi.setAlgorithms([alg])

    return SUCCESS


# =============================================================================
# The actual job steering
if '__main__' == __name__:

    input = [
        '/afs/cern.ch/user/i/ibelyaev/cmtuser/Gauss_v42r1/Gauss-14145019-1000ev-20120911.gen',
        '/afs/cern.ch/user/i/ibelyaev/cmtuser/Gauss_v42r1/Gauss-14145019-10000ev-20120911.gen']

    configure(input, castor=True)

    run(1000)

    gaudi = appMgr()

# =============================================================================
# The END
# =============================================================================
