#!/usr/bin/env ipython
# ========================================================================
# @file DiCharm/Alg
#
#  Simple template algorithm for 2xCharm
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
import DiCharm.Pids  # add methods for Pid/Track information
# =============================================================================
# @class DiCharmAlg
#  Simple template algorithm to study 2xCharm
#  @date   2011-05-27
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class DiCharmAlg (Algo):

    """
    Simple algorithm to study J/psi+``D'' production
    """

    def initialize(self):
        """
        Initialization
        """
        sc = Algo.       initialize(self)
        if sc.isFailure():
            return sc

        triggers = {}
        from DiCharm.TisTos import lines

        sc = self.tisTos_initialize(triggers, lines)
        if sc.isFailure():
            return sc

        sc = self.   pid_initialize()
        if sc.isFailure():
            return sc

        import DiCharm.GoodParticles as GP

        self.goodD0 = GP.basicD0() & GP.pidPion() & GP.pidKaon()
        self.goodDp = GP.basicDp() & GP.pidPion() & GP.pidKaon()
        self.goodDs = GP.basicDs() & GP.pidPion() & GP.pidKaon()
        self.goodLc = GP.basicLc() & GP.pidPion(
        ) & GP.pidKaon() & GP.pidProton()
        self.goodJpsi = GP.basicJpsi() & GP.pidMuon()
        self.goodDiMu = GP.basicDiMu() & GP.pidMuon()
        self.goodY = GP.basicY() & GP.pidMuon()
        self.goodW = GP.basicW()
        self.goodZ = GP.basicZ()
        #
        self.prompt = GP.prompt()
        self.prompt = PALL

        _pid = LoKi.Particles.pidFromName
        self.dimu_pids = [
            (ADMASS('J/psi(1S)'),  _pid('J/psi(1S)')),
            (ADMASS('psi(2S)'),  _pid('psi(2S)')),
            (ADMASS('Upsilon(1S)'),  _pid('Upsilon(1S)')),
            (ADMASS('Upsilon(2S)'),  _pid('Upsilon(2S)')),
            (ADMASS('Upsilon(3S)'),  _pid('Upsilon(3S)')),
        ]

        return SUCCESS

    # set new particle ID for dimuons,
    #  accoring to nearest resonance
    def dimuon_ID(self, dimu):
        """
        Set new particle ID for dimuons,
        accoring to nearest resonance     
        """

        pid = 0
        m_ = -1
        for i in self.dimu_pids:
            dm = i[0](dimu)
            if dm < m_ or m_ < 0:
                m_ = dm
                pid = i[1]

        dimu.setParticleID(pid)

        return pid

    # fill N-tuple
    def addItem2(self,
                 tup,
                 pair,
                 d1, name1, trg1, fun1,
                 d2, name2, trg2, fun2):

        #
        # the major kinematical properties
        self.treatKine(tup, pair, '_' + '2c')  # di charm
        self.treatKine(tup, d1, '_' + name1, fun1)
        self.treatKine(tup, d2, '_' + name2, fun2)

        #
        # promptness
        tup.column('prompt_' + '2c', self.prompt(pair), 0, 1)
        tup.column('prompt_' + name1, self.prompt(d1), 0, 1)
        tup.column('prompt_' + name2, self.prompt(d2), 0, 1)

        #
        # add the information needed for TisTos
        self.tisTos(d1, tup, name1 + '_',
                    self.lines[trg1], self.l0tistos, self.tistos)
        self.tisTos(d2, tup, name2 + '_',
                    self.lines[trg2], self.l0tistos, self.tistos)

        #
        # add the information for Pid efficiency correction
        self.treatPions(tup, pair)
        self.treatKaons(tup, pair)
        self.treatProtons(tup, pair)
        self.treatMuons(tup, pair)

        #
        # add the infomation needed for track efficiency correction
        self.treatTracks(tup,  pair)

        return SUCCESS

    # fill N-tuple
    def addItem(self,
                tup,
                pair,
                name1, trg1, fun1,
                name2, trg2, fun2):

        d1 = pair.child(1)
        d2 = pair.child(2)

        return self.addItem2(tup, pair,
                             d1, name1, trg1, fun1,
                             d2, name2, trg2, fun2)

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

        # D0
        tup = self.nTuple('PsiD0')
        psiD0 = self.select(
            'PsiD0', '[Meson -> J/psi(1S)  ( D0 -> K- pi+)]CC ')
        for pair in psiD0:
            self.addItem(tup,
                         pair,
                         'psi', 'psi', self.goodJpsi,
                         'D0', 'D0', self.goodD0)
            #
            # add some reco-summary information
            self.addRecSummary(tup, rc_summary)
            #
            tup.write()

        # D+
        tup = self.nTuple('PsiDp')
        psiDp = self.select(
            'PsiDp', '[Meson -> J/psi(1S)  ( D+ -> K- pi+ pi+)]CC ')
        for pair in psiDp:
            self.addItem(tup,
                         pair,
                         'psi', 'psi', self.goodJpsi,
                         'Dp', 'D+', self.goodDp)
            #
            # add some reco-summary information
            self.addRecSummary(tup, rc_summary)
            #
            tup.write()

        # Ds+
        tup = self.nTuple('PsiDs')
        psiDs = self.select(
            'PsiDs', '[Meson -> J/psi(1S)  ( D_s+ -> K- K+ pi+)]CC ')
        for pair in psiDs:
            self.addItem(tup,
                         pair,
                         'psi', 'psi', self.goodJpsi,
                         'Ds', 'Ds+', self.goodDs)
            #
            # add some reco-summary information
            self.addRecSummary(tup, rc_summary)
            #
            tup.write()

        # Lambda_c+
        tup = self.nTuple('PsiLc')
        psiLc = self.select(
            'PsiLc', '[Meson -> J/psi(1S)  ( Lambda_c+ -> p+ K- pi+)]CC ')
        for pair in psiLc:
            self.addItem(tup,
                         pair,
                         'psi', 'psi', self.goodJpsi,
                         'Lc', 'Lc+', self.goodLc)

            #
            # add some reco-summary information
            self.addRecSummary(tup, rc_summary)
            #
            tup.write()

        return SUCCESS

    # finalize it!
    def finalize(self):
        """
        finalize the algorithm
        """
        #
        self.goodD0 = None
        self.goodDp = None
        self.goodDs = None
        self.goodLc = None
        self.goodJpsi = None
        self.prompt = None
        #
        self.goodY = None
        self.goodW = None
        self.goodZ = None
        #
        self . tisTos_finalize()
        self .    pid_finalize()
        #
        self . dumpHistos()
        #
        return Algo.finalize(self)

# =============================================================================
# @class DiDAlg
#  Simple template algorithm to study 2xCharm
#  @date   2011-05-27
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class DiDAlg (DiCharmAlg):

    """
    Simple algorithm to study ``DD'' production
    """
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
        rc_summary = self.get('/Event/Charm/Rec/Summary').summaryData()

        ## D0 &D0
        D0D0 = self.select(
            'D0&D0',
            '[ ( Meson -> ( D0 -> K- pi+ ) ( D0 -> K- pi+ ) ) || ( Meson -> ( D0 -> K- pi+ ) ( D~0 -> K+ pi- ) ) ]CC'
        )
        tup = self.nTuple('D0D0')
        for pair in D0D0:
            self.addItem(tup,
                         pair,
                         'D01', 'D0', self.goodD0,
                         'D02', 'D0', self.goodD0)
            #
            # add some reco-summary information
            self.addRecSummary(tup, rc_summary)
            #
            tup.write()

        # D0 & D+
        D0Dp = self.select(
            'D0&D+',
            '[ ( Meson -> ( D0 -> K- pi+ ) D+ ) || ( Meson -> ( D0 -> K- pi+ ) D- ) ]CC'
        )
        tup = self.nTuple('D0Dp')
        for pair in D0Dp:
            self.addItem(tup,
                         pair,
                         'D0', 'D0', self.goodD0,
                         'Dp', 'D+', self.goodDp)
            #
            # add some reco-summary information
            self.addRecSummary(tup, rc_summary)
            #
            tup.write()

        # D0 & Ds+
        D0Ds = self.select(
            'D0&Ds+',
            '[ ( Meson -> ( D0 -> K- pi+ ) D_s+ ) || ( Meson -> ( D0 -> K- pi+ ) D_s- ) ]CC'
        )
        tup = self.nTuple('D0Ds')
        for pair in D0Ds:
            self.addItem(tup,
                         pair,
                         'D0', 'D0', self.goodD0,
                         'Ds', 'Ds+', self.goodDs)
            #
            # add some reco-summary information
            self.addRecSummary(tup, rc_summary)
            #
            tup.write()

        # D0 & Lambda_c+
        D0Lc = self.select(
            'D0&Lc+',
            '[ ( Meson -> ( D0 -> K- pi+ ) Lambda_c+ ) || ( Meson -> ( D0 -> K- pi+ ) Lambda_c~- ) ]CC'
        )
        tup = self.nTuple('D0Lc')
        for pair in D0Lc:
            self.addItem(tup,
                         pair,
                         'D0', 'D0', self.goodD0,
                         'Lc', 'Lc+', self.goodLc)
            #
            # add some reco-summary information
            self.addRecSummary(tup, rc_summary)
            #
            tup.write()

        # D+ & D+
        DpDp = self.select(
            'D+&D+',
            '[ ( Meson -> D+ D- ) || ( Meson -> D+ D+ ) ]CC'
        )
        tup = self.nTuple('DpDp')
        for pair in DpDp:
            self.addItem(tup,
                         pair,
                         'Dp1', 'D+', self.goodDp,
                         'Dp2', 'D+', self.goodDp)
            #
            # add some reco-summary information
            self.addRecSummary(tup, rc_summary)
            #
            tup.write()

        # D+ & Ds+
        DpDs = self.select(
            'D+&Ds+',
            '[ ( Meson -> D+ D_s- ) || ( Meson -> D+ D_s+ ) ]CC'
        )
        tup = self.nTuple('DpDs')
        for pair in DpDs:
            self.addItem(tup,
                         pair,
                         'Dp', 'D+', self.goodDp,
                         'Ds', 'Ds+', self.goodDs)
            #
            # add some reco-summary information
            self.addRecSummary(tup, rc_summary)
            #
            tup.write()

        # D+ & Lambda_C+
        DpLc = self.select(
            'D+&Lc+',
            '[ ( Meson -> D+ Lambda_c~- ) || ( Meson -> D+ Lambda_c+ ) ]CC'
        )
        tup = self.nTuple('DpLc')
        for pair in DpLc:
            self.addItem(tup,
                         pair,
                         'Dp', 'D+', self.goodDp,
                         'Lc', 'Lc+', self.goodLc)
            #
            # add some reco-summary information
            self.addRecSummary(tup, rc_summary)
            #
            tup.write()

        # Ds+ & Ds+
        DsDs = self.select(
            'Ds+&Ds+',
            '[ ( Meson -> D_s+ D_s- ) || ( Meson -> D_s+ D_s+ ) ]CC'
        )
        tup = self.nTuple('DsDs')
        for pair in DsDs:
            self.addItem(tup,
                         pair,
                         'Ds1', 'Ds+', self.goodDs,
                         'Ds2', 'Ds+', self.goodDs)
            #
            # add some reco-summary information
            self.addRecSummary(tup, rc_summary)
            #
            tup.write()

        # Ds+ & Lambda_C+
        DsLc = self.select(
            'Ds+&Lc+',
            '[ ( Meson -> D_s+ Lambda_c~- ) || ( Meson -> D_s+ Lambda_c+ ) ]CC'
        )
        tup = self.nTuple('DsLc')
        for pair in DsLc:
            self.addItem(tup,
                         pair,
                         'Ds', 'Ds+', self.goodDs,
                         'Lc', 'Lc+', self.goodLc)

            #
            # add some reco-summary information
            self.addRecSummary(tup, rc_summary)
            #
            tup.write()

        # Lambda_c+ & Lambda_c+
        LcLc = self.select(
            'Lc+&Lc+',
            '[ ( Meson -> Lambda_c+ Lambda_c~- ) || ( Meson -> Lambda_c+ Lambda_c+ ) ]CC'
        )
        tup = self.nTuple('LcLc')
        for pair in LcLc:
            self.addItem(tup,
                         pair,
                         'Lc1', 'Lc+', self.goodLc,
                         'Lc2', 'Lc+', self.goodLc)
            #
            # add some reco-summary information
            self.addRecSummary(tup, rc_summary)
            #
            tup.write()

        return SUCCESS

# =============================================================================
# @class Di2MuAlg
#  Simple template algorithm to study 2xdi-muon
#  @date   2011-05-27
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class Di2MuAlg (DiCharmAlg):

    """
    Simple algorithm to study ``2xdi-muon'' production
    """

    def initialize(self):
        """
        Initialization
        """
        sc = DiCharmAlg.initialize(self)
        if sc.isFailure():
            return sc

        self.delta_m2 = LoKi.PhysKinematics.deltaM2

        dimu1 = LoKi.Child.Selector(1)
        mu11 = LoKi.Child.Selector(1, 1)
        mu12 = LoKi.Child.Selector(1, 2)

        dimu2 = LoKi.Child.Selector(2)
        mu21 = LoKi.Child.Selector(2, 1)
        mu22 = LoKi.Child.Selector(2, 2)

        self.cos_theta1 = COSPOL(mu11, dimu1)
        self.cos_theta2 = COSPOL(mu21, dimu2)
        self.angle_chi = ANGLECHI(mu11, mu12,
                                  mu21, mu22)

        self.triggers['2xDiMuon'] = {}
        self.triggers['J/psi(1S)'] = {}
        self.triggers['psi(2S)'] = {}
        self.triggers['Upsilon(1S)'] = {}
        self.triggers['Upsilon(2S)'] = {}
        self.triggers['Upsilon(3S)'] = {}

        return SUCCESS

    # the only one esential method:
    def analyse(self):
        """
        The only one essential method
        """

        primaries = self.vselect('PVs', ISPRIMARY)
        if primaries.empty():
            return self.Warning('No primary vertices are found', SUCCESS)
        #
        # counters
        # gec_counters = self.get('/Event/Charm/GlobalEventCounters').counters()
        #
        # rec summary
        rc_summary = self.get('/Event/Rec/Summary').summaryData()

        # 2xdimu
        fourmu = self.select(
            '2xdi-mu', ' Meson ->  ( Meson -> mu+ mu-)  ( Meson -> mu+ mu-) ')

        tup = self.nTuple('JJ')

        for pair in fourmu:

            # reset IDs
            p1 = pair.child(1)
            p2 = pair.child(2)
            self.dimuon_ID(p1)
            self.dimuon_ID(p2)

            # collect trigger info
            self.decisions(pair,
                           self.triggers['2xDiMuon'],
                           self.l0tistos,
                           self.  tistos)
            self.decisions(p1,
                           self.triggers[p1.name()],
                           self.l0tistos,
                           self.  tistos)
            self.decisions(p2,
                           self.triggers[p2.name()],
                           self.l0tistos,
                           self.  tistos)

            # fill n-tuple
            self.addItem(tup,
                         pair,
                         'dimu1', 'psi', self.goodDiMu,
                         'dimu2', 'psi', self.goodDiMu)
            # some global event information
            # self.addGecInfo    ( tup , gec_counters )
            # add some reco-summary information
            self.addRecSummary(tup, rc_summary)
            #

            m11 = p1.child(1)
            m12 = p1.child(2)

            m21 = p2.child(1)
            m22 = p2.child(2)

            tup.column('dm2_1', self.delta_m2(m11, m21))  # OK
            tup.column('dm2_2', self.delta_m2(m12, m22))  # OK
            tup.column('dm2_3', self.delta_m2(m11, m22))  # cross-check
            tup.column('dm2_4', self.delta_m2(m12, m21))  # cross-check

            tup.column('cos_t1', self.cos_theta1(pair))
            tup.column('cos_t2', self.cos_theta2(pair))
            tup.column('angchi', self.angle_chi(pair))

            tup.write()

        return SUCCESS

    # finalize it!
    def finalize(self):
        """
        finalize the algorithm
        """
        #
        self.cos_theta1 = None
        self.cos_theta2 = None
        self.angle_chi = None
        #
        return DiCharmAlg.finalize(self)

# =============================================================================
# configure the job


def configure_common(rootInTES,
                     datafiles,
                     catalogs=[],
                     castor=False):
    """
    Job configuration 
    """

    from Configurables import DaVinci  # needed for job configuration
    from Configurables import EventSelector  # needed for job configuration
    from Configurables import FileCatalog  # needed for job configuration
    from Configurables import NTupleSvc  # needed for job configuration

    from PhysConf.Filters import LoKi_Filters

    fltrs = LoKi_Filters(
        STRIP_Code="""
        HLT_PASS_RE ( 'Stripping.*CharmAndDiMuon.*Decision' ) |
        HLT_PASS_RE ( 'Stripping.*DiMuonAndCharm.*Decision' ) |
        HLT_PASS_RE ( 'Stripping.*DiCharm.*Decision'        ) | 
        HLT_PASS_RE ( 'Stripping.*DoubleDiMuon.*Decision'   ) 
        """
    )
    filters = fltrs.filters('Filters')
    # filters.reverse()

    davinci = DaVinci(
        DataType='2011',
        InputType='MDST',
        Simulation=False,
        PrintFreq=1000,
        EventPreFilters=filters,
        EvtMax=-1,
        #
        HistogramFile='2xCharm_Histos.root',
        TupleFile='2xCharm.root',
        #
        Lumi=True,
        #
    )

    from Configurables import CondDB
    CondDB(UseLatestTags=["2011"],
           ## UseOracle     = True
           )

    # from Configurables import Gaudi__IODataManager as IODataManager
    # IODataManager().AgeLimit = 2

    # ------- decoding set-up start ----------
    from Bender.MicroDST import uDstConf
    uDstConf(rootInTES)
    logger.info("Configure RootInTES:  %s " % rootInTES)
    # ------- decoding set-up end  -----------

    # come back to Bender
    setData(datafiles, catalogs, castor)

    return SUCCESS


# =============================================================================
# configure the job
def configure_JD(datafiles, catalogs=[], castor=False):
    """
    Job configuration 
    """
    #
    # common configuration
    #
    rootInTES = '/Event/Charm'
    configure_common(
        rootInTES,
        datafiles,
        catalogs,
        castor)

    #
    # start Gaudi
    #
    gaudi = appMgr()

    #
    # Set properties of the TisTosTools
    #
    for t in (gaudi.tool('DiCharm.L0TriggerTisTos'),
              gaudi.tool('DiCharm.TriggerTisTos')):

        t . UseParticle2LHCbIDsMap = 2
        t . PropertiesPrint = True

    alg = DiCharmAlg(
        'DiCharm',  # Algorithm name ,
        RootInTES=rootInTES,
        Inputs=['Phys/DiMuonAndCharmForPromptCharm/Particles']
    )
    #
    mainSeq = gaudi.algorithm('GaudiSequencer/DaVinciUserSequence', True)
    mainSeq.Members += [alg.name()]

    return SUCCESS

# =============================================================================
# configure the job


def configure_DD(datafiles, catalogs=[], castor=False):
    """
    Job configuration 
    """
    #
    # common configuration
    #
    rootInTES = '/Event/Charm'
    configure_common(
        rootInTES,
        datafiles,
        catalogs,
        castor)

    #
    # start Gaudi
    #
    gaudi = appMgr()

    #
    # Set properties of the TisTosTools
    #
    for t in (gaudi.tool('DiD.L0TriggerTisTos'),
              gaudi.tool('DiD.TriggerTisTos')):

        t . UseParticle2LHCbIDsMap = 2
        t . PropertiesPrint = True

    alg = DiDAlg(
        'DiD',  # Algorithm name ,
        RootInTES=rootInTES,
        Inputs=['Phys/DiCharmForPromptCharm/Particles']
    )
    #
    mainSeq = gaudi.algorithm('GaudiSequencer/DaVinciUserSequence', True)
    mainSeq.Members += [alg.name()]

    return SUCCESS


# =============================================================================
# configure the job
def configure_JJ(datafiles, catalogs=[], castor=False):
    """
    Job configuration 
    """
    #
    # common configuration
    #
    rootInTES = '/Event/Leptonic'
    configure_common(rootInTES,
                     datafiles,
                     catalogs,
                     castor)

    #
    # start Gaudi
    #
    gaudi = appMgr()

    #
    # Set properties of the TisTosTools
    #
    for t in (gaudi.tool('Di2Mu.L0TriggerTisTos'),
              gaudi.tool('Di2Mu.TriggerTisTos')):

        t . UseParticle2LHCbIDsMap = 2
        t . PropertiesPrint = True

    #
    alg = Di2MuAlg(
        'Di2Mu',  # Algorithm name ,
        RootInTES=rootInTES,
        Inputs=['Phys/DoubleDiMuonForCharmAssociative/Particles']
    )
    #
    mainSeq = gaudi.algorithm('GaudiSequencer/DaVinciUserSequence', True)
    mainSeq.Members += [alg.name()]

    return SUCCESS


# =============================================================================
# configure the job
def configure(datafiles, catalogs=[], castor=False):
    """
    Configure the job
    """
    # return configure_JD ( datafiles , catalogs , castor )
    # return configure_DD ( datafiles , catalogs , castor )
    return configure_JJ(datafiles, catalogs, castor)

# =============================================================================
# The actual job steering
if '__main__' == __name__:

    import DiCharm.Samples as Samples
    from DiCharm.Utils import lfn2disk

    data = Samples.DATA_2011_1fb_samples

    # dicharm        = lfn2disk (
    # Samples.DATA_samples [ 'DiCharm/Down'     ] +
    #   Samples.DATA_samples [ 'DiCharm/Up'       ]  )

    # psi_and_charm  = lfn2disk (
    #    Samples.DATA_samples [ 'DiMuon&Charm/Down'] +
    #    Samples.DATA_samples [ 'DiMuon&Charm/Up'  ]
    #    )

    ## lfns = Samples.DATA_samples [ 'DiMuon&Charm/Down']
    ## + Samples.DATA_samples [ 'DiMuon&Charm/Up'  ]

    psi_and_psi = data['2xDiMuon/Up']
    psi_and_charm = data['DiMuon&Charm/Up']
    charm_and_charm = data['DiCharm/Up']

    configure(psi_and_charm, [], True)

    run(100)

    gaudi = appMgr()

# =============================================================================
# The END
# =============================================================================
