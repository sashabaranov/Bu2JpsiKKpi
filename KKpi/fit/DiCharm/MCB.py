#!/usr/bin/env ipython
# =============================================================================
# @file DiCharm/Alg.py
#
#  Get from MC the background from B-decays
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

Get from MC the background from B-decays 

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
#
from GaudiKernel.SystemOfUnits import GeV, MeV, mm
from GaudiKernel.PhysicalConstants import c_light
# =============================================================================
import DiCharm.Efficiency as Eff
# =============================================================================
# @class MCB
#  Simple template algorithm to get background from B-decays
#  @date   2011-05-27
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class MCB(AlgoMC):

    """
    Simple algorithm to study feeddown from B-decays 
    """
    # the only one esential method:

    def initialize(self):

        sc = AlgoMC .     initialize(self)
        if sc.isFailure():
            return sc

        import DiCharm.GoodParticles as GP

        self.prompt = GP.prompt()
        self.dtfchi2 = DTF_CHI2NDOF(True)

        return SUCCESS

    # the only one esential method:
    def analyse(self):
        """
        The only one essential method
        """

        tup = self.nTuple('MCB')

        psi = self.  select(
            'Psi', '[Meson -> J/psi(1S)                   Charm ]CC ')
        D0 = self.  select(
            'D0', '[Meson -> ( D0        -> K- pi+     ) Charm ]CC ')
        Dp = self.  select(
            'Dp', '[Meson -> ( D+        -> K- pi+ pi+ ) Charm ]CC ')
        Ds = self.  select(
            'Ds', '[Meson -> ( D_s+      -> K- K+  pi+ ) Charm ]CC ')
        Lc = self.  select(
            'Lc', '[Meson -> ( Lambda_c+ -> p+ K-  pi+ ) Charm ]CC ')
        charm = self.  select('c2', '[Meson -> Charm Charm ]CC ')

        mcPsi = self.mcselect('mcPsi', '   J/psi(1S)  => mu+ mu-        ')
        mcD0 = self.mcselect('mcD0', ' [ D0         -> K- pi+     ]CC ')
        mcDp = self.mcselect('mcDp', ' [ D+        --> K- pi+ pi+ ]CC ')
        mcDs = self.mcselect('mcDs', ' [ D_s+      --> K- K+  pi+ ]CC ')
        mcLc = self.mcselect('mcLc', ' [ Lambda_c+ --> p+ K-  pi+ ]CC ')

        MCpsi = MCTRUTH(
            self.mcTruth(), mcPsi) if not mcPsi.empty() else PNONE
        MCd0 = MCTRUTH(
            self.mcTruth(), mcD0) if not mcD0 .empty() else PNONE
        MCdp = MCTRUTH(
            self.mcTruth(), mcDp) if not mcDp .empty() else PNONE
        MCds = MCTRUTH(
            self.mcTruth(), mcDs) if not mcDs .empty() else PNONE
        MClc = MCTRUTH(
            self.mcTruth(), mcLc) if not mcLc .empty() else PNONE

        idD0 = 'D0' == ABSID
        idDp = 'D+' == ABSID
        idDs = 'D_s+' == ABSID
        idLc = 'Lambda_c+' == ABSID

        idPsi = 'J/psi(1S)' == ABSID
        idDst = 'D*(2010)+' == ABSID

        tup = self.nTuple('C2')

        for cc in charm:

            c1 = cc.child(1)
            c2 = cc.child(2)

            if idDst(c1) or idDst(c2):
                continue

            tup.column('pid1', ID(c1))
            tup.column('pid2', ID(c2))

            pr = self.prompt(cc)
            tup.column('promt', pr, 0, 1)

            dtf = self.dtfchi2(cc)
            tup.column('dtf', dtf)

            tup.column('pt1', PT(c1) / GeV)
            tup.column('pt2', PT(c2) / GeV)

            tup.column('mc_1_Psi', MCpsi(c1), 0, 1)
            tup.column('mc_1_D0', MCd0(c1), 0, 1)
            tup.column('mc_1_Dp', MCdp(c1), 0, 1)
            tup.column('mc_1_Ds', MCds(c1), 0, 1)
            tup.column('mc_1_Lc', MClc(c1), 0, 1)

            tup.column('mc_2_Psi', MCpsi(c2), 0, 1)
            tup.column('mc_2_D0', MCd0(c2), 0, 1)
            tup.column('mc_2_Dp', MCdp(c2), 0, 1)
            tup.column('mc_2_Ds', MCds(c2), 0, 1)
            tup.column('mc_2_Lc', MClc(c2), 0, 1)

            trg1 = Eff.VE(1.0, 0.0)

            pt1 = PT(c1) / GeV
            y1 = Y(c1)
            pt2 = PT(c2) / GeV
            y2 = Y(c2)

            trg1 = Eff.VE(0.0, 0.0)
            ok1 = False
            if idPsi(c1) and MCpsi(c1):
                trg1 = Eff.tosEff_Jpsi(pt1, y1, 'L0xL1xL2')
                ok1 = 0 <= pt1 <= 12 and 2 <= y1 < 4.5
            elif idD0(c1) and MCd0(c1):
                trg1 = Eff.tosEff_D0(pt1, y1, 'L0xL1xL2')
                ok1 = 2 <= pt1 <= 12 and 2 <= y1 < 4.0
            elif idDp(c1) and MCdp(c1):
                trg1 = Eff.tosEff_Dp(pt1, y1, 'L0xL1xL2')
                ok1 = 2 <= pt1 <= 12 and 2 <= y1 < 4.0
            elif idDs(c1) and MCds(c1):
                trg1 = Eff.tosEff_Ds(pt1, y1, 'L0xL1xL2')
                ok1 = 2 <= pt1 <= 12 and 2 <= y1 < 4.0
            elif idLc(c1) and MClc(c1):
                trg1 = Eff.tosEff_Lc(pt1, y1, 'L0xL1xL2')
                ok1 = 2 <= pt1 <= 12 and 2 <= y1 < 4.0

            trg2 = Eff.VE(0.0, 0.0)
            ok2 = False
            if idPsi(c2) and MCpsi(c2):
                trg2 = Eff.tosEff_Jpsi(pt2, y2, 'L0xL1xL2')
                ok2 = 0 <= pt2 <= 12 and 2 <= y2 < 4.5
            elif idD0(c2) and MCd0(c2):
                trg2 = Eff.tosEff_D0(pt2, y2, 'L0xL1xL2')
                ok2 = 2 <= pt2 <= 12 and 2 <= y2 < 4.0
            elif idDp(c2) and MCdp(c2):
                trg2 = Eff.tosEff_Dp(pt2, y2, 'L0xL1xL2')
                ok2 = 2 <= pt2 <= 12 and 2 <= y2 < 4.0
            elif idDs(c2) and MCds(c2):
                trg2 = Eff.tosEff_Ds(pt2, y2, 'L0xL1xL2')
                ok2 = 2 <= pt2 <= 12 and 2 <= y2 < 4.0
            elif idLc(c2) and MClc(c2):
                trg2 = Eff.tosEff_Lc(pt2, y2, 'L0xL1xL2')
                ok2 = 2 <= pt2 <= 12 and 2 <= y2 < 4.0

            tup.column('trg1', trg1.value())
            tup.column('trg2', trg2.value())

            tup.column('ok1', ok1, 0, 1)
            tup.column('ok2', ok2, 0, 1)

            tup.column('m1', M(c1) / GeV)
            tup.column('m2', M(c2) / GeV)
            tup.column('m2c', M(cc) / GeV)

            print cc.decay(), pr, dtf, PT(c1) / GeV, PT(c2) / GeV, M(cc) / GeV

            tup.write()

        self.setFilterPassed(False)

        return SUCCESS

    def finalize(self):

        self.prompt = None
        self.dtfchi2 = None

        return AlgoMC.finalize(self)

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

    davinci = DaVinci(
        DataType='2011',
        InputType='DST',
        Simulation=True,
        PrintFreq=1000,
        ##EventPreFilters = filters ,
        EvtMax=-1,
        #
        HistogramFile='Dp_2xCharm_fromB_Histos.root',
        TupleFile='Dp_2xCharm_fromB.root',
        #
        Lumi=False
        #
    )

    # come back to Bender
    setData(datafiles, catalogs)

    gaudi = appMgr()

    alg = MCB(
        #
        'MCB',  # Algorithm name
        #
        Inputs=[
            '/Event/CharmAndDiMuonLine/Phys/SelDiMuonAndCharmForPromptCharm/Particles',
            '/Event/DiCharmLine/Phys/SelDiCharmForPromptCharm/Particles'
        ],
        #
        PP2MCs=['Relations/Rec/ProtoP/Charged']
    )
    #
    mainSeq = gaudi.algorithm('GaudiSequencer/DaVinciUserSequence', True)
    mainSeq.Members += [alg.name()]

    return SUCCESS
# =============================================================================
# The actual job steering
if '__main__' == __name__:

    import DiCharm.Samples as Samples
    from DiCharm.Utils import lfn2disk

    psi_from_B = lfn2disk(Samples.MC10_from_B['Jpsi_from_B'])
    D0_from_B = lfn2disk(Samples.MC10_from_B['D0_from_B'])
    Dp_from_B = lfn2disk(Samples.MC10_from_B['Dp_from_B'])
    Ds_from_B = lfn2disk(Samples.MC10_from_B['Ds_from_B'])
    Lc_from_B = lfn2disk(Samples.MC10_from_B['Lc_from_B'])

    configure(Dp_from_B)

    run(100)

    gaudi = appMgr()

# =============================================================================
# The END
# =============================================================================
