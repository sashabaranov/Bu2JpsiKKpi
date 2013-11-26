#!/usr/bin/env ipython
# ========================================================================
# @file
#
#  Simple algorithm to study trigger effiicency for double charm
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

Simple algorithm to study trigger efficiency for double charm 

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

import DiCharm.TisTos

# ========================================================================
# @class TrgCheck
#  Simple algorithm to study trigger efficiencies
#  @date   2010-12-01
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class TrgCheck(Algo):

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

        self.l0tistos = self.tool(
            cpp.ITriggerTisTos, 'L0TriggerTisTos', parent=self)
        self.tistos = self.tool(
            cpp.ITriggerTisTos,   'TriggerTisTos', parent=self)

        self.tistos   .setOfflineInput()
        self.l0tistos .setOfflineInput()

        self.l0trigs = {}
        self.trigs = {}

        self.l0trigs['TIS'] = {}
        self.l0trigs['TOS'] = {}

        self.  trigs['TIS'] = {}
        self.  trigs['TOS'] = {}

        self.L0TOSlines = 'L0(Muon|DiMuon)Decision'
        self.L0TISlines = 'L0(Muon|DiMuon|Electron|Hadron|Photon)Decision'

        self.Hlt1TOSlines = 'Hlt1(DiMuon|SingleMuon|MuTrack).*Decision'
        self.Hlt1TISlines = 'Hlt1(DiMuon|SingleMuon|MuTrack|DiHadron).*Decision'

        self.Hlt2TOSlines = 'Hlt2(DiMuon|SingleMuon|MuTrack).*Decision'
        self.Hlt2TISlines = self.Hlt2TOSlines

        self.triggers = {}
        self.triggers['psi'] = {}
        self.triggers['D0'] = {}
        self.triggers['D+'] = {}
        self.triggers['Ds+'] = {}
        self.triggers['D*+'] = {}
        self.triggers['Lc+'] = {}

        return SUCCESS

    # the only one esential method:
    def analyse(self):
        """
        The only one essential method
        """

        primaries = self.vselect('PVs', ISPRIMARY)
        if primaries.empty():
            self.Warning('No primary vertices are found', SUCCESS)

        psi = self.select(
            'psi', '  Meson -> ^( J/psi(1S) -> mu+ mu-   ) Charm    ')
        d0 = self.select(
            'd0', '[ Meson -> ^( D0        -> K-  pi+   ) Charm]CC ')
        dp = self.select(
            'd+', '[ Meson -> ^( D+   -> K- pi+ pi+     ) Charm]CC ')
        ds = self.select(
            'ds', '[ Meson -> ^( D_s+ -> K- K+  pi+     ) Charm]CC ')
        lc = self.select(
            'Lc', '[ Meson -> ^( Lambda_c+ -> p+ K- pi+ ) Charm]CC ')

        # attention!!!
        dst = self.select(
            'dst', '[ Meson ->  ( D*(2010)+ -> ^( D0 -> K- pi+ ) pi+ ) Charm]CC ')

        for d in psi:
            self.decisions(
                d, self.triggers['psi'], self.l0tistos, self.tistos)
        for d in d0:
            self.decisions(
                d, self.triggers['D0'], self.l0tistos, self.tistos)
        for d in dp:
            self.decisions(
                d, self.triggers['D+'], self.l0tistos, self.tistos)
        for d in ds:
            self.decisions(
                d, self.triggers['Ds+'], self.l0tistos, self.tistos)
        for d in lc:
            self.decisions(
                d, self.triggers['Lc+'], self.l0tistos, self.tistos)
        for d in dst:
            self.decisions(
                d, self.triggers['D*+'], self.l0tistos, self.tistos)

        if psi.empty() and d0.empty() and dp.empty() and ds.empty() and dst.empty() and lc.empty():
            return self.Warning('No interesting particles are found', SUCCESS)

        return SUCCESS

    def finalize(self):

        histos = self.Histos()
        for key in histos:
            h = histos[key]
            if hasattr(h, 'dump'):
                print h.dump(50, 30, True)

        self.trgDecs(self.triggers)

        self.l0tistos = None
        self.tistos = None

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
        HLT_PASS_RE ( 'Stripping.*CharmAndDiMuon.*Decision'  ) |
        HLT_PASS_RE ( 'Stripping.*DiCharm.*Decision'         )
        """ ,
        VOID_Code="""
        ( 0.5 < CONTAINS('/Event/Charm/Phys/CharmAndDiMuonForPromptCharm/Particles' ) )
        |
        ( 0.5 < CONTAINS('/Event/Charm/Phys/DiCharmForPromptCharm/Particles' ) )        
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
        HistogramFile='Charm_Histos.root',
        TupleFile='Charm.root',
        #
        Lumi=True,
        #
    )

    from Configurables import CondDB
    CondDB().UseLatestTags = ["2011"]

    # from Configurables import Gaudi__IODataManager as IODataManager
    # IODataManager().AgeLimit = 2

    # configure TIS-TOS
    rootInTES = '/Event/Charm'

    # ------- decoding set-up start ----------

    from MicroDSTConf.TriggerConfUtils import configureL0AndHltDecoding
    configureL0AndHltDecoding(rootInTES)

    # ------- decoding set-up end  -----------

    # come back to Bender
    setData(datafiles, catalogs)

    gaudi = appMgr()

    # Set properties of the TisTosTools
    for t in (gaudi.tool('TrgCheck.L0TriggerTisTos'),
              gaudi.tool('TrgCheck.TriggerTisTos')):

        t . UseParticle2LHCbIDsMap = 2
        t . PropertiesPrint = True

    alg = TrgCheck(
        'TrgCheck',  # Algorithm name ,
        RootInTES=rootInTES,
        Inputs=['Phys/CharmAndDiMuonForPromptCharm/Particles',
                'Phys/DiCharmForPromptCharm/Particles']

    )
    #
    mainSeq = gaudi.algorithm('GaudiSequencer/DaVinciUserSequence', True)
    mainSeq.Members += [alg.name()]

    return SUCCESS


# =============================================================================
# The actual job steering
if '__main__' == __name__:

    data1_3 = [
        '/castor/cern.ch/grid/lhcb/user/i/ibelyaev/2011_07/22614/22614404/MicroDST.DiCharmLine.mdst',
        '/castor/cern.ch/grid/lhcb/user/i/ibelyaev/2011_07/22614/22614427/MicroDST.DiCharmLine.mdst',
        ##'/castor/cern.ch/grid/lhcb/user/i/ibelyaev/2011_07/22614/22614427/MicroDST.CharmAndDiMuonLine.mdst' ,
    ]

    files = data1_3

    configure(files)

    run(100)

    gaudi = appMgr()


# =============================================================================
# The END
# =============================================================================
