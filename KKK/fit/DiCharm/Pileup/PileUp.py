#!/usr/bin/env ipython
# =============================================================================
# @file
#
#  Simple algorithm to check the pileup for 2xCharm production
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
#  @date   2011-07-11
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#
#                    $Revision$
#  Last modification $Date$
#                 by $Author$
# =============================================================================
"""

Simple algorithm to check the pileup for 2xCharm production 

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
__date__ = "2011-07-11"
__version__ = '$Revision$'
# =============================================================================
## needed to produce/visualize the histograms
import ROOT
from Bender.MainMC import *  # import all bender goodies
## easy access to various geometry routines
import LHCbMath.Types
from Gaudi.Configuration import *  # needed for job configuration

from GaudiKernel.SystemOfUnits import GeV, MeV, mm
from GaudiKernel.PhysicalConstants import c_light
# =============================================================================

# =============================================================================
# @class PsiPU
#  check for psi + charm pileup
#  @date   2011-07-11
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class PsiPU (AlgoMC):

    """
    Check for psi+``D'' pileup 
    """

    # collect charm hadrons
    def getOpenCharm(self):

        return self.mcselect('charm',
                             in_range(2.0, MCY, 4.0) &
                             in_range(2 * GeV, MCPT, 12 * GeV) &
                             (('D0' == MCABSID) |
                              ('D+' == MCABSID) |
                                 ('D_s+' == MCABSID) |
                              ('Lambda_c+' == MCABSID))
                             )

    def analyse(self):

        # collect open charm:
        charm = self.getOpenCharm()

        if charm.empty():
            return self.Warning("No charm is found", SUCCESS)

        # collect Jpsi:
        psis__ = self.mcselect('psi__', 'J/psi(1S) => mu+ mu-')
        psis = self.mcselect('psi', psis__,
                             in_range(2.0, MCY, 4.5) &
                             (MCPT < 12 * GeV))

        if psis.empty():
            return self.Warning("No good J/psi is found", SUCCESS)

        self.fillTuple(psis, charm, 'mcJpsiC')

        self.setFilterPassed(True)

        return SUCCESS

    # fill n-tuple
    def fillTuple(self, signal, charm, tname):

        beauty = self.mcselect('b',  BEAUTY)

        fromB = MCINANCESTORS(BEAUTY)

        pv_keys = set()
        for p in signal:
            pv = p.primaryVertex()
            if pv:
                pv_keys.add(pv.key())
        for p in charm:
            pv = p.primaryVertex()
            if pv:
                pv_keys.add(pv.key())

        tup0 = self.nTuple(tname)

        tup0.fArrayMCP('mcid', MCID,  beauty, 'nBmc', 100)
        tup0.column('npv', len(pv_keys), 0, 100)
        tup0.fArrayMCP(
            's_id', MCID,
            's_pt', MCPT,
            's_y', MCY,
            's_fromB', switch(fromB, 1, 0),
            signal, 'nSignal', 100)
        tup0.fArrayMCP(
            'c_id', MCID,
            'c_pt', MCPT,
            'c_y', MCY,
            'c_fromB', switch(fromB, 1, 0),
            charm, 'nC', 100)
        tup0.fArrayMCP(
            's_px', MCPX,
            's_py', MCPY,
            's_pz', MCPZ,
            signal, 'nSignal', 100)
        tup0.fArrayMCP(
            'c_px', MCPX,
            'c_py', MCPY,
            'c_pz', MCPZ,
            charm, 'nC', 100)
        tup0.fArrayMCP(
            's_vx', MCVFASPF(MCVX),
            's_vy', MCVFASPF(MCVY),
            's_vz', MCVFASPF(MCVZ),
            's_vk', MCVFASPF(MCVKEY),
            signal, 'nSignal', 100)
        tup0.fArrayMCP(
            'c_vx', MCVFASPF(MCVX),
            'c_vy', MCVFASPF(MCVY),
            'c_vz', MCVFASPF(MCVZ),
            'c_vk', MCVFASPF(MCVKEY),
            charm, 'nC', 100)
        tup0.fArrayMCP(
            's_pvx', MCVPXFUN(MCVX),
            's_pvy', MCVPXFUN(MCVY),
            's_pvz', MCVPXFUN(MCVZ),
            's_pvk', MCVPXFUN(MCVKEY),
            signal, 'nSignal', 100)
        tup0.fArrayMCP(
            'c_pvx', MCVPXFUN(MCVX),
            'c_pvy', MCVPXFUN(MCVY),
            'c_pvz', MCVPXFUN(MCVZ),
            'c_pvk', MCVPXFUN(MCVKEY),
            charm, 'nC', 100)

        tup0.write()

        return SUCCESS


# ============================================================================
# @class D0PU
#  check for D0 + charm pileup
#  @date   2011-07-11
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
class D0PU (PsiPU):

    """
    Check for D0+``D'' pileup 
    """

    def analyse(self):

        d0_ = self.mcselect('d0_', '[D0 => K- pi+]CC')
        d0s = self.mcselect('d0', d0_,
                            in_range(2.0, MCY, 4.0) &
                            in_range(2 * GeV, MCPT, 12 * GeV)
                            )

        if d0s.empty():
            return self.Warning("No good D0 is found", SUCCESS)

        charm = self.getOpenCharm()

        self.fillTuple(d0s, charm, 'mcD0C')

        self.setFilterPassed(not charm.empty())

        return SUCCESS

# ============================================================================
# @class DpPU
#  check for D+ + charm pileup
#  @date   2011-07-11
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class DpPU (PsiPU):

    """
    Check for D+ + ``D'' pileup 
    """

    def analyse(self):

        dp_ = self.mcselect('dp_', '[D+ --> K- pi+ pi+]CC')
        dps = self.mcselect('dp', dp_,
                            in_range(2.0, MCY, 4.0) &
                            in_range(2 * GeV, MCPT, 12 * GeV)
                            )

        if dps.empty():
            return self.Warning("No good D+ is found", SUCCESS)

        charm = self.getOpenCharm()

        self.fillTuple(dps, charm, 'mcDpC')

        self.setFilterPassed(not charm.empty())

        return SUCCESS

# ============================================================================
# @class DsPU
#  check for Ds+ + charm pileup
#  @date   2011-07-11
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class DsPU (PsiPU):

    """
    Check for Ds+ + ``D'' pileup 
    """

    def analyse(self):

        ds_ = self.mcselect('ds_', '[ D_s+ --> K- K+ pi+]CC')
        dss = self.mcselect('ds', ds_,
                            in_range(2.0, MCY, 4.0) &
                            in_range(2 * GeV, MCPT, 12 * GeV)
                            )

        if dss.empty():
            return self.Warning("No good D_s+ is found", SUCCESS)

        charm = self.getOpenCharm()

        self.fillTuple(dss, charm, 'mcDsC')

        self.setFilterPassed(not charm.empty())

        return SUCCESS

# ============================================================================
# @class LcPU
#  check for Lc+ + charm pileup
#  @date   2011-07-11
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class LcPU (PsiPU):

    """
    Check for Lc+ + ``D'' pileup 
    """

    def analyse(self):

        lc_ = self.mcselect('lc_', '[Lambda_c+ --> p+ K- pi+]CC')
        lcs = self.mcselect('lc', lc_,
                            in_range(2.0, MCY, 4.0) &
                            in_range(2 * GeV, MCPT, 12 * GeV)
                            )

        if lcs.empty():
            return self.Warning("No good Lc+ is found", SUCCESS)

        charm = self.getOpenCharm()

        self.fillTuple(lcs, charm, 'mcLcC')

        self.setFilterPassed(not charm.empty())

        return SUCCESS

# =============================================================================
# configure the job


def configure_common(alg_type, datafiles, catalogs=[]):
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

    davinci = DaVinci(
        DataType='2010',
        InputType='DST',
        Simulation=True,
        PrintFreq=1000,
        ## EventPreFilters = filters ,
        EvtMax=-1,
        #
        HistogramFile='MCPU_2xCharm_Histos.root',
        TupleFile='MCPU_2xCharm.root',
        #
        Lumi=True,
        #
    )

    from Configurables import CondDB
    CondDB(IgnoreHeartBeat=True)

    setData(datafiles, catalogs)

    gaudi = appMgr()

    alg = alg_type(
        'PU',  # Algorithm name ,
        Inputs=[],
        PP2MCs=[]
    )
    #
    gaudi.setAlgorithms([alg])

    return SUCCESS


# =============================================================================
# configure the job
def configure(datafiles, catalogs=[]):
    """
    Job configuration 
    """

    return configure_common(PsiPU, datafiles,  catalogs)
    # return configure_common ( D0PU , datafiles ,  catalogs )
    # return configure_common ( DpPU , datafiles ,  catalogs )
    # return configure_common ( LcPU , datafiles ,  catalogs )
    # return configure_common ( DsPU , datafiles ,  catalogs )

# =============================================================================
# The actual job steering
if '__main__' == __name__:

    data1_3 = [
        #
        # jpsi
        #'/castor/cern.ch/grid/lhcb/MC/MC10/ALLSTREAMS.DST/00008941/0000/00008941_00000%03d_1.allstreams.dst' % i for i in range ( 1 , 90 )
        #
        # [ D*+ -> ( D0 => K- pi+ ) pi+ ]CC
        #"/castor/cern.ch/grid/lhcb/MC/MC10/ALLSTREAMS.DST/00008554/0000/00008554_00000%03d_1.allstreams.dst" % i for i in range(1,300)
        # D+ -> K- pi+ pi+  21263010
        #'/castor/cern.ch/grid/lhcb/MC/MC10/ALLSTREAMS.DST/00009404/0000/00009404_00000%03d_1.allstreams.dst' % i for i in range(1,400)
        # Lambda_c+ -> p K- pi+ 25103000
        # "/castor/cern.ch/grid/lhcb/MC/MC10/ALLSTREAMS.DST/00011118/0000/00011118_00000%03d_1.allstreams.dst" % i for i in range(1,255)
        # D_s+ -> K+ K- pi+   23263001
        "/castor/cern.ch/grid/lhcb/MC/MC10/ALLSTREAMS.DST/00011124/0000/00011124_00000%03d_1.allstreams.dst" % i for i in range(1, 600)
    ]

    files = data1_3

    configure(files)

    run(100)

    gaudi = appMgr()


# =============================================================================
# The END
# =============================================================================
