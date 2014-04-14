#!/usr/bin/env ipython
# ========================================================================
# @file
#
#  Simple algorithm to get trigger efficiencies for charm particles
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

Simple algorithm to study trigger efficiencies for charm particles 

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
import BenderTools.TisTos  # add methods for tistos
import DiCharm.Pids  # add methods for Pid/Track information
# ========================================================================

# ========================================================================
# @class TrgEff
#  Simple algorithm to study trigger efficiencies for charm particles
#  @date   2010-12-01
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class TrgEff(Algo):

    """
    Simple algorithm to study trigger efficiencies 
    """
    # constructor

    def __init__(self,
                 name='TrgEff',
                 **args):
        """
        Constructor
        """
        Algo.__init__(self, name, **args)

    def initialize(self):
        """
        Initialization
        """
        sc = Algo.initialize(self)
        if sc.isFailure():
            return sc

        #
        # trigger lists
        #
        triggers = {}
        triggers["D0"] = {}
        triggers["D*+"] = {}
        triggers["D+"] = {}
        triggers["Ds+"] = {}
        triggers["Lc+"] = {}
        triggers["psi"] = {}
        triggers["psi'"] = {}
        triggers["Y"] = {}

        #
        # list of trigger lines for TisTos
        #
        from DiCharm.TisTos import lines

        sc = self.tisTos_initialize(triggers=triggers,
                                    lines=lines)

        if sc.isFailure():
            return sc

        sc = self.   pid_initialize()
        if sc.isFailure():
            return sc
        #
        import DiCharm.GoodParticles as GP
        #
        self.goodD0 = GP.basicD0() & GP.prompt(
        ) & GP.pidKaon() & GP.pidPion()
        self.goodDst = GP.basicDstar() & GP.prompt(
        ) & GP.pidKaon() & CHILDCUT(1, GP.pidPion())
        self.goodDp = GP.basicDplus() & GP.prompt(
        ) & GP.pidKaon() & GP.pidPion()
        self.goodDs = GP.basicDs() & GP.prompt(
        ) & GP.pidKaon() & GP.pidPion()
        self.goodLc = GP.basicLamC() & GP.prompt(
        ) & GP.pidKaon() & GP.pidProton() & GP.pidPion()

        self.goodPsi = GP.basicJpsi() & GP.prompt() & GP.pidMuon()
        self.goodPsiPrime = GP.basicPsiPrime() & GP.prompt() & GP.pidMuon()
        self.goodY = GP.basicY() & GP.prompt() & GP.pidMuon()

        return SUCCESS

    # add record to n-tuple
    def addRec(self,
               tup,
               part,
               id,
               trg):

        self.treatKine(tup, part,  '_' + id)

        self.tisTos(part,
                    tup,
                    id + '_',
                    self.  lines[trg],
                    self.l0tistos,
                    self.  tistos,
                    self.triggers[trg], True)

        # add the information for Pid efficiency correction
        self.treatPions(tup, part)
        self.treatKaons(tup, part)
        self.treatProtons(tup, part)
        self.treatMuons(tup, part)
        # add the infomation needed for track efficiency correction
        self.treatTracks(tup,  part)

        return SUCCESS

    # the only one esential method:
    def analyse(self):
        """
        The only one essential method
        """
        primaries = self.vselect('PVs', ISPRIMARY)
        if primaries.empty():
            self.Warning('No primary vertices are found', SUCCESS)

        d0_ = self.select('d0_', '[D0   -> K-  pi+ ]CC ')
        d0 = self.select('d0', d0_, self.goodD0)

        dp_ = self.select('dp_', '[D+   -> K-  pi+ pi+ ]CC ')
        dp = self.select('dp', dp_, self.goodDp)

        ds_ = self.select('ds_', '[D_s+ -> K- K+ pi+ ]CC ')
        ds = self.select('ds', ds_, self.goodDs)

        lc_ = self.select('lc_', '[Lambda_c+ -> p+ K- pi+ ]CC ')
        lc = self.select('lc', lc_, self.goodLc)

        # counters
        ## gec_counters = self.get('/Event/Charm/GlobalEventCounters').counters()

        # rec summary
        rc_summary = self.get('/Event/Charm/Rec/Summary').summaryData()

        tD0 = self.nTuple('D0')
        tDp = self.nTuple('Dplus')
        tDs = self.nTuple('Ds')
        tLc = self.nTuple('Lc')

        for d in d0:
            self.decisions(
                d, self.triggers['D0'], self.l0tistos, self.tistos)
            self.addRec(tD0, d, 'D0', 'D0')
            #
            self.addRecSummary(tD0, rc_summary)
            #
            tD0.write()

        for d in dp:
            self.decisions(
                d, self.triggers['D+'], self.l0tistos, self.tistos)
            self.addRec(tDp, d, 'Dp', 'D+')
            #
            self.addRecSummary(tDp, rc_summary)
            #
            tDp.write()

        for d in ds:
            self.decisions(
                d, self.triggers['Ds+'], self.l0tistos, self.tistos)
            self.addRec(tDs, d, 'Ds', 'Ds+')
            #
            self.addRecSummary(tDs, rc_summary)
            #
            tDs.write()

        for d in lc:
            self.decisions(
                d, self.triggers['Lc+'], self.l0tistos, self.tistos)
            self.addRec(tLc, d, 'Lc', 'Lc+')
            #
            self.addRecSummary(tLc, rc_summary)
            #
            tLc.write()

        return SUCCESS

    def finalize(self):

        histos = self.Histos()
        for key in histos:
            h = histos[key]
            if hasattr(h, 'dump'):
                print h.dump(50, 30, True)

        self.tisTos_finalize()
        self.   pid_finalize()

        self.goodD0 = None
        self.goodDst = None
        self.goodDp = None
        self.goodDs = None
        self.goodLc = None
        self.goodPsi = None
        self.goodPsiPrime = None
        self.goodY = None

        return Algo.finalize(self)

# ========================================================================
# @class TrgPsiEff
#  Simple algorithm to study trigger efficiencies for J/psi
#  @date   2010-12-01
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class TrgPsiEff(TrgEff):

    """
    Simple algorithm to study trigger efficiencies 
    """

    #
    # the only one esential method:
    def analyse(self):
        """
        The only one essential method
        """
        primaries = self.vselect('PVs', ISPRIMARY)
        if primaries.empty():
            self.Warning('No primary vertices are found', SUCCESS)

        psi_ = self.select('psi_', 'Meson -> mu+ mu- ')
        psi = self.select('psi', psi_, self.goodPsi)

        if psi.empty():
            self.Warning("No good J/psi are are found", SUCCESS)

        tPsi = self.nTuple('Psi')

        odin = self.get('/Event/Leptonic/DAQ/ODIN')

        # rec summary
        rc_summary = self.get('/Event/Leptonic/Rec/Summary').summaryData()

        for j in psi:

            self.plot(M(j) / GeV, 'm(J/psi)', 3.0, 3.2, 200)

            self.decisions(
                j, self.triggers['psi'], self.l0tistos, self.tistos)

            self.addRec(tPsi, j, 'psi', 'psi')
            self.addRecSummary(tPsi, rc_summary)

            tPsi.column_aux(odin)
            tPsi.write()

        return SUCCESS

# ========================================================================
# @class TrgPsiPEff
#  Simple algorithm to study trigger efficiencies for psi'
#  @date   2012-05-12
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru


class TrgPsiPEff(TrgEff):

    """
    Simple algorithm to study trigger efficiencies for psi'
    """

    #
    # the only one esential method:
    def analyse(self):
        """
        The only one essential method
        """
        primaries = self.vselect('PVs', ISPRIMARY)
        if primaries.empty():
            self.Warning('No primary vertices are found', SUCCESS)

        psip_ = self.select('psip_', 'Meson -> mu+ mu- ')
        psip = self.select('psip', psip_, self.goodPsiPrime)

        if psip.empty():
            self.Warning("No good psi' are are found", SUCCESS)

        tPsiP = self.nTuple('PsiP')

        # ODIN
        odin = self.get('/Event/Leptonic/DAQ/ODIN')

        # rec summary
        rc_summary = self.get('/Event/Leptonic/Rec/Summary').summaryData()

        for p in psip:
            #
            self.plot(M(p) / GeV, "m(psi_2s)", 3.5, 3.9, 200)
            #
            self.decisions(p,
                           self.triggers["psi'"],
                           self.l0tistos,
                           self.tistos)

            self.addRec(tPsiP,  p, 'psip', "psi'")
            self.addRecSummary(tPsiP, rc_summary)

            tPsiP.column_aux(odin)
            tPsiP.write()

        return SUCCESS

# ========================================================================
# @class TrgUpsEff
#  Simple algorithm to study trigger efficiencies for Upsilons
#  @date   2012-05-12
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru


class TrgUpsEff(TrgEff):

    """
    Simple algorithm to study trigger efficiencies for Upsilons
    """
    #
    # the only one esential method:

    def analyse(self):
        """
        The only one essential method
        """
        primaries = self.vselect('PVs', ISPRIMARY)
        if primaries.empty():
            self.Warning('No primary vertices are found', SUCCESS)

        ups_ = self.select('ups_', 'Meson -> mu+ mu- ')
        ups = self.select('ups', ups_, self.goodY)

        if ups.empty():
            self.Warning("No good Ys are are found", SUCCESS)

        tY = self.nTuple('Y')

        # ODIN
        odin = self.get('/Event/DAQ/ODIN')

        #
        # rec summary, here FullDST is used
        #
        rc_summary = self.get('/Event/Rec/Summary').summaryData()

        for y in ups:

            self.plot(M(y) / GeV, "m(Y)", 8.5, 12.5, 200)
            self.decisions(y,
                           self.triggers["Y"],
                           self.l0tistos,
                           self.tistos)

            self.addRec(tY, y, 'ups', "Y")

            self.addRecSummary(tY, rc_summary)

            tY . column_aux(odin)
            tY . write()

        return SUCCESS


# =============================================================================
# configure the job
def configure_Charm(datafiles, catalogs=[]):
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
        HLT_PASS_RE ( 'Stripping.*ForPromptCharm.*Decision'  ) 
        """,
        VOID_Code="""
        ( 0.5 < CONTAINS ('/Event/Charm/Phys/D02HHForPromptCharm/Particles'   ) )
        |
        ( 0.5 < CONTAINS ('/Event/Charm/Phys/DstarForPromptCharm/Particles'   ) ) 
        |
        ( 0.5 < CONTAINS ('/Event/Charm/Phys/DForPromptCharm/Particles'       ) )
        |
        ( 0.5 < CONTAINS ('/Event/Charm/Phys/DsForPromptCharm/Particles'      ) ) 
        |
        ( 0.5 < CONTAINS ('/Event/Charm/Phys/LambdaCForPromptCharm/Particles' ) )
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
        HistogramFile='TrgEff_Histos.root',
        TupleFile='TrgEff.root',
        #
        Lumi=True,
        #
    )

    from Configurables import CondDB
    CondDB(UseLatestTags=["2011"])

    # from Configurables import Gaudi__IODataManager as IODataManager
    # IODataManager().AgeLimit = 2

    #
    # configure TIS-TOS
    #
    rootInTES = '/Event/Charm'
    from MicroDSTConf.TriggerConfUtils import configureL0AndHltDecoding
    configureL0AndHltDecoding(rootInTES)

    # come back to Bender
    setData(datafiles, catalogs)

    gaudi = appMgr()

    # Set properties of the TisTosTools
    for t in (gaudi.tool('TrgEff.L0TriggerTisTos'),
              gaudi.tool('TrgEff.TriggerTisTos')):

        t . UseParticle2LHCbIDsMap = 2
        t . PropertiesPrint = True

    alg = TrgEff(
        RootInTES=rootInTES,
        Inputs=['Phys/D02HHForPromptCharm/Particles',
                'Phys/DstarForPromptCharm/Particles',
                'Phys/DForPromptCharm/Particles',
                'Phys/DsForPromptCharm/Particles',
                'Phys/LambdaCForPromptCharm/Particles']
    )
    #
    mainSeq = gaudi.algorithm('GaudiSequencer/DaVinciUserSequence', True)
    mainSeq.Members += [alg.name()]

    return SUCCESS

# =============================================================================
# configure the job


def configure_Psi(datafiles, catalogs=[], castor=False):
    """
    Job configuration 
    """

    ## needed for job configuration
    from Configurables import DaVinci

    from PhysConf.Filters import LoKi_Filters

    fltrs = LoKi_Filters(
        STRIP_Code="""
        HLT_PASS_RE ( 'Stripping.*DiMuonInc.*Decision'  )        
        """,
        VOID_Code="""
        ( 0.5 < CONTAINS ( '/Event/Leptonic/Phys/MicroDSTDiMuonDiMuonIncLine/Particles' ) ) & 
        ( 0.5 < CONTAINS ( '/Event/Leptonic/Rec/Vertex/Primary'                         ) )         
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
        HistogramFile='TrgEff_Histos.root',
        TupleFile='TrgEff.root',
        #
        Lumi=True,
        #
    )

    from Configurables import CondDB
    CondDB(UseLatestTags=["2011"])

    #
    # configure micro-DST treatment
    #
    rootInTES = '/Event/Leptonic'
    from BenderTools.MicroDST import uDstConf
    uDstConf(rootInTES)

    #
    # come back to Bender
    #
    setData(datafiles, catalogs, castor)

    gaudi = appMgr()

    alg = TrgPsiEff(
        RootInTES=rootInTES,
        Inputs=['Phys/MicroDSTDiMuonDiMuonIncLine/Particles']
    )
    #
    mainSeq = gaudi.algorithm('GaudiSequencer/DaVinciUserSequence', True)
    mainSeq.Members += [alg.name()]

    return SUCCESS


# =============================================================================
# configure the job
def configure_PsiP(datafiles, catalogs=[], castor=False):
    """
    Job configuration 
    """
    ## needed for job configuration
    from Configurables import DaVinci

    from PhysConf.Filters import LoKi_Filters

    fltrs = LoKi_Filters(
        STRIP_Code="""
        HLT_PASS_RE ( 'Stripping.*DiMuonInc.*Decision'  )        
        """,
        VOID_Code="""
        ( 0.5 < CONTAINS ( '/Event/Leptonic/Phys/MicroDSTDiMuonDiMuonIncLine/Particles' ) ) &
        ( 0.5 < CONTAINS ( '/Event/Leptonic/Rec/Vertex/Primary'                         ) )         
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
        HistogramFile='TrgEff_Histos.root',
        TupleFile='TrgEff.root',
        #
        Lumi=True,
        #
    )

    from Configurables import CondDB
    CondDB(UseLatestTags=["2011"])

    #
    # configure microDST treatment
    #
    rootInTES = '/Event/Leptonic'

    from BenderTools.MicroDST import uDstConf
    uDstConf(rootInTES)

    #
    # come back to Bender
    #
    setData(datafiles, catalogs, castor)

    gaudi = appMgr()

    alg = TrgPsiPEff(
        RootInTES=rootInTES,
        Inputs=['Phys/MicroDSTDiMuonDiMuonIncLine/Particles']
    )
    #
    mainSeq = gaudi.algorithm('GaudiSequencer/DaVinciUserSequence', True)
    mainSeq.Members += [alg.name()]

    return SUCCESS


# =============================================================================
# configure the job
def configure_Ups(datafiles, catalogs=[], castor=False):
    """
    Job configuration 
    """
    ## needed for job configuration
    from Configurables import DaVinci

    from PhysConf.Filters import LoKi_Filters

    fltrs = LoKi_Filters(
        STRIP_Code="""
        HLT_PASS_RE ( 'Stripping.*DiMuonHighMass.*Decision'  )
        """,
        VOID_Code="""
        ( 0.5 < CONTAINS ( '/Event/Dimuon/Phys/FullDSTDiMuonDiMuonHighMassLine/Particles' ) ) & 
        ( 0.5 < CONTAINS ( '/Event/Rec/Vertex/Primary'                                    ) )         
        """
    )
    filters = fltrs.filters('Filters')
    filters.reverse()

    davinci = DaVinci(
        DataType='2011',
        InputType='DST',
        Simulation=False,
        PrintFreq=1000,
        EventPreFilters=filters,
        EvtMax=-1,
        #
        HistogramFile='TrgEff_Histos.root',
        TupleFile='TrgEff.root',
        #
        Lumi=True,
        #
    )

    from Configurables import CondDB
    CondDB(UseLatestTags=["2011"])

    #
    # come back to Bender
    #
    setData(datafiles, catalogs, castor)

    gaudi = appMgr()

    # for t in ( 'TrgEff.L0TriggerTisTos',
    #           'TrgEff.TriggerTisTos'  ) :
    #
    #    t1 = gaudi.tool ( t )
    #    t1.TOSFracMuon     =  0.01
    #    t1.PropertiesPrint = True

    alg = TrgUpsEff(
        'TrgEff',  # Algorithm name ,
        Inputs=[
            '/Event/Dimuon/Phys/FullDSTDiMuonDiMuonHighMassLine/Particles']
    )
    #
    mainSeq = gaudi.algorithm('GaudiSequencer/DaVinciUserSequence', True)
    mainSeq.Members += [alg.name()]

    return SUCCESS

# =============================================================================
# configure the job


def configure(datafiles, catalogs=[], castor=False):
    """
    Job configuration 
    """
    # return configure_Charm ( datafiles , catalogs , castor )
    return configure_Psi(datafiles, catalogs, castor)
    # return configure_PsiP ( datafiles , catalogs , castor )
    # return configure_Ups ( datafiles , catalogs , castor )

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
        # Charm microDST
        # "/castor/cern.ch/grid/lhcb/LHCb/Collision11/CHARM.MDST/00010646/0000/00010646_000000%02d_1.charm.mdst" % i for i in range(1,20)
        # Leptonic microDST, stripping 13
        #"/castor/cern.ch/grid/lhcb/LHCb/Collision11/LEPTONIC.MDST/00010831/0000/00010831_0000%04d_1.leptonic.mdst" % i for i in range(1,1800)
        # Leptonic microDST, stripping 17
        "/castor/cern.ch/grid/lhcb/LHCb/Collision11/LEPTONIC.MDST/00012551/0000/00012551_0000%04d_1.leptonic.mdst" % i for i in range(1, 100)
        # Dimuon DST, stripping 17
        #"/castor/cern.ch/grid/lhcb/LHCb/Collision11/DIMUON.DST/00012549/0000/00012549_0000%04d_1.dimuon.dst" % i for i in range(1,100)
    ]

    files = data1_3

    configure(files)

    run(1000)

    gaudi = appMgr()

    alg = gaudi.algorithm('TrgEff')
    alg.dumpHistos()

# =============================================================================
# The END
# =============================================================================
