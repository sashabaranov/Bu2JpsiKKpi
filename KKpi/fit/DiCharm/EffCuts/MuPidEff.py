#!/usr/bin/env ipython
# ========================================================================
# @file
#
#  Simple algorithm to get muID efficiency for J/psi
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

Simple algorithm to get muID efficiency for J/psi 

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
# ========================================================================
## needed to produce/visualize the histograms
import ROOT
from Bender.Main import *  # import all bender goodies
## easy access to various geometry routines
import LHCbMath.Types
from Gaudi.Configuration import *  # needed for job configuration

from GaudiKernel.SystemOfUnits import GeV, MeV, mm
from GaudiKernel.PhysicalConstants import c_light
# ========================================================================

# ========================================================================
# @class MuPidEff
#  Simple algorithm to get muID efficiency for J/psi
#  @date   2011-07-16
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class MuPidEff(Algo):

    """
    Simple algorithm to study trigger efficiencies 
    """

    def initialize(self):
        """
        Initialization
        """
        sc = Algo.initialize(self)
        if sc.isFailure():
            return sc

        import DiCharm.GoodParticles as GP

        self.goodPsi = GP.basicJpsi() & GP.prompt()

        return SUCCESS

    # the only one esential method:
    def analyse(self):
        """
        The only one essential method
        """
        primaries = self.vselect('PVs', ISPRIMARY)
        if primaries.empty():
            self.Warning('No primary vertices are found', SUCCESS)

        psi_ = self.select('psi_', 'J/psi(1S) -> mu+ mu- ')
        psi = self.select('psi', psi_, self.goodPsi)

        tPsi = self.nTuple('Psi')

        for j in psi:

            tPsi.column('mass',  M(j) / GeV)
            tPsi.column('pt', PT(j) / GeV)
            tPsi.column('y',  Y(j))
            tPsi.column('lv01', LV01(j))

            mu1 = j.child(1)
            mu2 = j.child(2)

            tPsi.column('DLLmu1', PIDmu(mu1))
            tPsi.column('DLLmu2', PIDmu(mu2))
            tPsi.column('pt1', PT(mu1) / GeV)
            tPsi.column('pt2', PT(mu2) / GeV)

            tPsi.write()

        return SUCCESS

    def finalize(self):

        self.goodPsi = None

        return Algo.finalize(self)

# =============================================================================
# configure the job


def configure(datafiles, catalogs=[]):
    """
    Job configuration 
    """

    ## needed for job configuration
    from Configurables import DaVinci
    ## needed for job configuration
    from Configurables import EventSelector
    ## needed for job configuration
    from GaudiConf.Configuration import FileCatalog
    ## needed for job configuration
    from GaudiConf.Configuration import NTupleSvc

    from PhysConf.Filters import LoKi_Filters

    fltrs = LoKi_Filters(
        STRIP_Code="""
        HLT_PASS_RE ( 'Stripping.*DiMuonInc.*Decision'  )
        """,
        VOID_Code="""
        0.5 < CONTAINS ('/Event/Leptonic/Phys/MicroDSTDiMuonDiMuonIncLine/Particles' )
        """
    )
    filters = fltrs.filters('Filters')
    filters.reverse()

    davinci = DaVinci(
        DataType='2011',
        InputType='MDST',
        Simulation=False,
        PrintFreq=1000,
        EventPreFilters=filters,
        EvtMax=-1,
        #
        HistogramFile='MuPidEff_Histos.root',
        TupleFile='MuPidEff.root',
        #
        Lumi=True,
        #
    )

    from Configurables import CondDB
    CondDB(UseLatestTags=["2011"],
           UseOracle=True
           )

    # from Configurables import Gaudi__IODataManager as IODataManager
    # IODataManager().AgeLimit = 2

    # configure TIS-TOS
    rootInTES = '/Event/Leptonic'

    # ------- decoding set-up start ----------

    ## from MicroDSTConf.TriggerConfUtils import configureL0AndHltDecoding
    # configureL0AndHltDecoding(rootInTES)

    # ------- decoding set-up end  -----------

    # come back to Bender

    setData(datafiles, catalogs)

    gaudi = appMgr()

    alg = MuPidEff(
        'MuPidEff',  # Algorithm name ,
        RootInTES=rootInTES,
        Inputs=['Phys/MicroDSTDiMuonDiMuonIncLine/Particles']
    )
    #
    mainSeq = gaudi.algorithm('GaudiSequencer/DaVinciUserSequence', True)
    mainSeq.Members += [alg.name()]

    return SUCCESS


# =============================================================================
# The actual job steering
if '__main__' == __name__:

    print 80 * '*'
    print __doc__
    print ' Author  : ', __author__
    print ' Version : ', __version__
    print ' Date    : ', __date__
    print 80 * '*'

    data1_3 = [
        # Leptonic microDST
        "/castor/cern.ch/grid/lhcb/LHCb/Collision11/LEPTONIC.MDST/00010831/0000/00010831_0000%04d_1.leptonic.mdst" % i for i in range(1, 1800)
    ]

    files = data1_3

    configure(files)

    run(100)

    gaudi = appMgr()


# =============================================================================
# The END
# =============================================================================
