#!/bin/env python

#=============================================================================
__author__ = "  "
__date__ = "  "
__version__ = "  "
#=============================================================================

import __builtin__

from math import *

# import everything from BENDER
from Bender.All import *
import BenderTools.Fill
import BenderTools.TisTos

from Gaudi.Configuration import *

# needed for job configuration

# import of the useful units from Gaudi
from GaudiKernel.SystemOfUnits import GeV, MeV, mm, micrometer
from GaudiKernel.PhysicalConstants import c_light
from LoKiTracks.decorators import *  # needed for TrKEY work

import LHCbMath.Types                 ## easy access to various geometry routines

import __builtin__

from LoKiCore.basic import cpp


class fakeK ( object ) :

    def __init__(self, p, pid) :
        self.particle = p
        self.old_pid  = LHCb.ParticleID ( p.particleID() )
        self.new_pid  = pid

    def __enter__  ( self ) :
        self.particle.setParticleID ( self.new_pid )

    def __exit__   ( self , *_ ) :

        self.particle.setParticleID ( self.old_pid )
        self.particle = None



#=============================================================================

AlgoMC.decisions         = BenderTools.TisTos. decisions
AlgoMC.trgDecs           = BenderTools.TisTos. trgDecs
AlgoMC.tisTos            = BenderTools.TisTos. tisTos
AlgoMC.tisTos_initialize = BenderTools.TisTos._tisTosInit
AlgoMC.tisTos_finalize   = BenderTools.TisTos._tisTosFini


class MCAnalysisAlgorithm(AlgoMC):

    def initialize(self):
        triggers = {}
        triggers ['psi'] = {}

        lines = {}
        lines [ "psi"      ] = {}
        lines [ "psi1"     ] = {} ## clean-selection: dimuon
        lines [ "psi2"     ] = {} ## clean-selection: detached
        lines [ "psi3"     ] = {} ## clean-selection: unbiased

        #
        ## J/psi
        #
        lines [ "psi" ][   'L0TOS' ] = 'L0(DiMuon|Muon)Decision'
        lines [ "psi" ][   'L0TIS' ] = 'L0(Hadron|DiMuon|Muon|Electron|Photon)Decision'
        lines [ "psi" ][ 'Hlt1TOS' ] = 'Hlt1(DiMuon|SingleMuon|TrackMuon).*Decision'
        lines [ "psi" ][ 'Hlt1TIS' ] = 'Hlt1(DiMuon|SingleMuon|Track).*Decision'
        lines [ "psi" ][ 'Hlt2TOS' ] = 'Hlt2(DiMuon|ExpressJPsi|SingleMuon).*Decision'
        lines [ "psi" ][ 'Hlt2TIS' ] = 'Hlt2(Charm|Topo|DiMuon|Single).*Decision'
        #
        ## J/psi ("clean-1")
        #
        lines [ "psi1" ][   'L0TOS' ] = 'L0(DiMuon|Muon)Decision'
        lines [ "psi1" ][   'L0TIS' ] = 'L0(Hadron|DiMuon|Muon|Electron|Photon)Decision'
        lines [ "psi1" ][ 'Hlt1TOS' ] = 'Hlt1(DiMuon|TrackMuon).*Decision'
        lines [ "psi1" ][ 'Hlt1TIS' ] = 'Hlt1(DiMuon|SingleMuon|Track).*Decision'
        lines [ "psi1" ][ 'Hlt2TOS' ] = 'Hlt2DiMuon.*Decision'
        lines [ "psi1" ][ 'Hlt2TIS' ] = 'Hlt2(Charm|Topo|DiMuon|Single).*Decision'
        #
        ## J/psi ("clean-2")
        #
        lines [ "psi2" ][   'L0TOS' ] = 'L0(DiMuon|Muon)Decision'
        lines [ "psi2" ][   'L0TIS' ] = 'L0(Hadron|DiMuon|Muon|Electron|Photon)Decision'
        lines [ "psi2" ][ 'Hlt1TOS' ] = 'Hlt1DiMuonHighMass.*Decision'
        lines [ "psi2" ][ 'Hlt1TIS' ] = 'Hlt1(DiMuon|SingleMuon|Track).*Decision'
        lines [ "psi2" ][ 'Hlt2TOS' ] = 'Hlt2DiMuonDetached(Heavy|JPsi)Decision'
        lines [ "psi2" ][ 'Hlt2TIS' ] = 'Hlt2(Charm|Topo|DiMuon|Single).*Decision'
        #
        ## J/psi ("unbiased")
        #
        lines [ "psi3" ][   'L0TOS' ] = 'L0(DiMuon|Muon)Decision'
        lines [ "psi3" ][   'L0TIS' ] = 'L0(Hadron|DiMuon|Muon|Electron|Photon)Decision'
        lines [ "psi3" ][ 'Hlt1TOS' ] = 'Hlt1DiMuonHighMass.*Decision'
        lines [ "psi3" ][ 'Hlt1TIS' ] = 'Hlt1(DiMuon|SingleMuon|Track).*Decision'
        lines [ "psi3" ][ 'Hlt2TOS' ] = 'Hlt2DiMuonJPsiHighPT.*Decision'
        lines [ "psi3" ][ 'Hlt2TIS' ] = 'Hlt2(Charm|Topo|DiMuon|Single).*Decision'

        sc = AlgoMC.initialize(self)
        if sc.isFailure():
            return sc

        sc = self.fill_initialize()
        if sc.isFailure():
            return sc

        sc = self.tisTos_initialize ( triggers , lines )
        if sc.isFailure () :
            return sc

        self._mass = DTF_FUN ( M , True , 'J/psi(1S)' )

        return SUCCESS

    def analyse(self):
        primaries = self.vselect( 'PVs' , ISPRIMARY )
        if primaries.empty() :
            return self.Warning('No primary vertices are found', SUCCESS )

        mcB = self.mcselect(
            'mcB', "[( B+ ==>  ( J/psi(1S) =>  mu+  mu-  )  K+  pi+  pi- )]CC")

        mcB_psi = self.mcselect(
            'mcB_psi', "[( B+ ==> K+ (psi(2S) => ( J/psi(1S) =>  mu+  mu-  ) pi+ pi-))]CC")

        mcB_X = self.mcselect(
            'mcB_X', "[( B+ ==> K+ (X_1(3872) => ( J/psi(1S) =>  mu+  mu-  ) pi+ pi-))]CC")

        mcB_K_Ks = self.mcselect(
            'mcB_K_Ks', "[( B+ ==>  ( J/psi(1S) =>  mu+  mu-  )  (K_1(1270)+  => (K*(892)0 => K+ pi-) pi+) )]CC")

        mcB_K_rho = self.mcselect(
            'mcB_K_rho', "[( B+ ==>  ( J/psi(1S) =>  mu+  mu-  )  (K_1(1270)+  => (rho(770)0 => pi+ pi-) K+) )]CC")

        mcB_K_K0s = self.mcselect(
            'mcB_K_K0s', "[( B+ ==>  ( J/psi(1S) =>  mu+  mu-  )  (K_1(1270)+  => (K*_0(1430)0 => K+ pi-) pi+) )]CC")


        sizes = [
            ("psi", mcB_psi.size()),
            ("X", mcB_X.size()),
            ("K_Ks", mcB_K_Ks.size()),
            ("K_rho", mcB_K_rho.size()),
            ("K_K0s", mcB_K_K0s.size()),
        ]

        nt_sizes = self.nTuple("sizes")

        for name, size in sizes:
            nt_sizes.column('mcB_' + name, size)
        
        nt_sizes.write()


        if mcB.size() != 1 or __builtin__.sum([x[1] for x in sizes]) != 1:
            return self.Warning("Something wrong with MC size " + str(mcB.size()), SUCCESS)


        mcK = self.mcselect(
            "mcK",  "[( B+ ==>  ( J/psi(1S) =>  mu+  mu-  )  ^K+  pi- pi+ )]CC")
        mcPi = self.mcselect(
            "mcPi",  "[( B+ ==>  ( J/psi(1S) =>  mu+  mu-  )  K+  ^pi- ^pi+ )]CC")
        mcMu = self.mcselect(
            "mcMu", "[( B+ ==>  ( J/psi(1S) =>  ^mu+  ^mu-  )  K+  pi- pi+ )]CC")
        mcPsi = self.mcselect(
            "mcPsi", "[( B+ ==>  ^( J/psi(1S) =>  mu+  mu-  )  K+  pi- pi+ )]CC")

        if mcK.empty() or mcMu.empty() or mcPsi.empty() or mcPi.empty():
            return self.Warning('No true MC-decay components are found', SUCCESS )

        match = self.mcTruth()
        trueB = MCTRUTH(match, mcB)

        true_mcB_psi = MCTRUTH(match, mcB_psi)
        true_mcB_X = MCTRUTH(match, mcB_X)
        true_mcB_K_Ks = MCTRUTH(match, mcB_K_Ks)
        true_mcB_K_rho = MCTRUTH(match, mcB_K_rho)
        true_mcB_K_K0s = MCTRUTH(match, mcB_K_K0s)


        trueK = MCTRUTH(match, mcK)
        truePi = MCTRUTH(match, mcPi)
        truePsi = MCTRUTH(match, mcPsi)
        trueMu = MCTRUTH(match, mcMu)


        kaons  = self.select ( 'K'  ,  ('K+'  == ABSID) & trueK )
        if kaons.empty():
            return self.Warning("No reconstructed kaons", SUCCESS)  # RETURN
    


        pions  = self.select ( 'pi'  ,  ('pi+'  == ABSID) & truePi )
        if pions.empty():
            return self.Warning("No reconstructed pions", SUCCESS)  # RETURN
        
        piplus  = self.select ( 'pi+' ,  pions , Q > 0 )
        piminus = self.select ( 'pi-' ,  pions , Q < 0 )
        if piplus.empty() or piminus.empty():
            return self.Warning("Both piplus and piminus are empty", SUCCESS)  # RETURN


        psis = self.select( "Jpsi", ('J/psi(1S)'  == ABSID) & truePsi )
        if psis.empty():
            return self.Warning("No reconstructed psis", SUCCESS)  # RETURN


        # myB = self.select('Bu' , '[( B+ ->  J/psi(1S)  K+  pi+  pi-)]CC' )
        myB = self.loop('Jpsi K pi+ pi-', 'B+')

        cut_pi = ( CLONEDIST   > 5000      ) & \
                ( TRGHOSTPROB < 0.5       ) & \
                ( TRCHI2DOF   < 4         ) & \
                in_range ( 2          , ETA , 5         ) & \
                in_range ( 3.2 * GeV  , P   , 150 * GeV ) & \
                HASRICH                     & \
                ( PROBNNpi     > 0.1      ) & \
                ( MIPCHI2DV()  > 4        )
        
        cut_k = ( PT          > 200 * MeV ) & \
                ( CLONEDIST   > 5000      ) & \
                ( TRGHOSTPROB < 0.5       ) & \
                ( TRCHI2DOF   < 4         ) & \
                in_range ( 2          , ETA , 5         ) & \
                in_range ( 3.2 * GeV  , P   , 150 * GeV ) & \
                HASRICH                     & \
                ( PROBNNk      > 0.1      ) & \
                ( MIPCHI2DV()  > 4        )

        for b in myB:
            if not 0 < VCHI2(b) < 100:
                continue
            k = b(2)
            if Q(k) > 0:
                b.setPID("B+")
            else:
                b.setPID("B-")

            # if not (cut_pi(b(3)) and cut_pi(b(4))):
            #     continue

            # if not cut_k(k):
            #     continue

            b.save("BB")

        bb = self.selected("BB")
        print len(bb)


        # Constrains
        dtffun_ctau = DTF_CTAU(0, True)
        dtffun_chi2 = DTF_CHI2NDOF(True, "J/psi(1S)")
        dtffun_m = DTF_FUN(M, True, "J/psi(1S)")


        nt = self.nTuple("t")

        for myb in bb:
            if not all([myb(i) for i in xrange(0, 5)]):
                continue


            b, jpsi, k, pi1, pi2 = tuple(myb(i) for i in xrange(5))

            # add DTF-applied information
            nt.column('DTFm_b', dtffun_m(myb) / GeV)
            nt.column('DTFctau', dtffun_ctau(myb))
            nt.column('DTFchi2ndof', dtffun_chi2(myb))

            MIPCHI2DVfun = MIPCHI2DV()

            self.treatKine(nt, b, '_b')
            self.treatKine(nt, jpsi, '_jpsi')

            # add the information for Pid efficiency correction
            self.treatPions(nt, b)
            self.treatKaons(nt, b)
            self.treatMuons(nt, b)
            self.treatTracks(nt, b)


            # ==========================================
            # Do fake
            # ==========================================
            nt.column('mass', self._mass ( b )  / GeV )

            ## try with pi1->K
            with fakeK ( pi2, pid = LHCb.ParticleID( int(Q(pi2)) * 321 ) ) :
                nt.column ( 'mass_pi2ask' , self._mass ( b ) / GeV )

            self.fillMasses(nt, myb, "c2", True, "J/psi(1S)")


            nt.column('MIPCHI2DV_k', MIPCHI2DVfun(k))
            nt.column('MIPCHI2DV_pi1', MIPCHI2DVfun(pi1))
            nt.column('MIPCHI2DV_pi2', MIPCHI2DVfun(pi2))

            nt.column ( 'mcTrueB'    , trueB(b)          )
            nt.column ( 'mcTrueB_psi'   , true_mcB_psi(myb))
            nt.column ( 'mcTrueB_X'   , true_mcB_X(myb))
            nt.column ( 'mcTrueB_K_Ks'   , true_mcB_K_Ks(myb))
            nt.column ( 'mcTrueB_K_rho'   , true_mcB_K_rho(myb))
            nt.column ( 'mcTrueB_K_K0s'   , true_mcB_K_K0s(myb))
            nt.column ( 'mcTruePsi' , truePsi(jpsi(0)    ))
            nt.column ( 'mcTrueK'    , trueK(myb(2))     )
            nt.column ( 'mcTruePi1'  , truePi(myb(3))    )
            nt.column ( 'mcTruePi2'  , truePi(myb(4))    )

            nt.column ( 'mcTrueMu1'  , trueMu(jpsi(1))    )
            nt.column ( 'mcTrueMu2'  , trueMu(jpsi(2))    )



            # add the information needed for TisTos
            self.tisTos ( jpsi  , nt  , 'psi_' ,
                          self.lines [ 'psi' ] , self.l0tistos , self.l1tistos , self.l2tistos )

            self.tisTos ( jpsi  , nt  , 'psi1_' ,
                          self.lines [ 'psi1' ] , self.l0tistos , self.l1tistos , self.l2tistos )

            self.tisTos ( jpsi  , nt  , 'psi2_' ,
                          self.lines [ 'psi2' ] , self.l0tistos , self.l1tistos , self.l2tistos )

            self.tisTos ( jpsi  , nt  , 'psi3_' ,
                          self.lines [ 'psi3' ] , self.l0tistos , self.l1tistos , self.l2tistos )

            nt.write()

        return SUCCESS

    # finalize & print histos
    def finalize(self):
        self.fill_finalize()
        self.tisTos_finalize ()

        self._mass = None
        return AlgoMC.finalize(self)

# =============================================================================
# configure the job


def configure(datafiles, catalogs=[], params={}, castor=False):
    """
    Configure the job
    """
    from Configurables           import DaVinci       ## needed for job configuration

    ## get the builder
    from StrippingSelections.StrippingPsiXForBandQ  import PsiX_BQ_Conf as  PSIX

    ## for MC it is better to exclude PID/DLL/PROBNN cuts
    builder_configuration = {
        # 'PionCut'   : """
        # ( PT          > 200 * MeV ) &
        'PionCut'   : """
        ( CLONEDIST   > 5000      ) &
        """
        # ( TRGHOSTPROB < 0.5       ) &
        # ( TRCHI2DOF   < 4         ) &
        # in_range ( 2          , ETA , 5         ) &
        # in_range ( 3.2 * GeV  , P   , 150 * GeV ) &
        # HASRICH                     &
        # ( PROBNNpi     > 0.1      )
        #( MIPCHI2DV()  > 4        )
        ,
        # ( PT          > 200 * MeV ) &
        'KaonCut'   : """
        ( CLONEDIST   > 5000      ) &
        ( TRGHOSTPROB < 0.5       ) &
        ( TRCHI2DOF   < 4         ) &
        in_range ( 2          , ETA , 5         ) &
        in_range ( 3.2 * GeV  , P   , 150 * GeV ) &
        HASRICH                     &
        ( PROBNNk      > 0.1      )
        """
        # ( MIPCHI2DV()  > 4        )
    }


    def _kaons_     ( self ) :
        """
        Kaons for   B -> psi X lines
        """
        from GaudiConfUtils.ConfigurableGenerators import FilterDesktop
        ## from StandardParticles                     import StdAllLooseKaons as inpts
        from StandardParticles                     import StdNoPIDsKaons   as inpts
        ##
        return self.make_selection (
            'Kaon'                 ,
            FilterDesktop          ,
            [ inpts ]              ,
            Code = self['KaonCut'] ,
        )


    def _pions_    ( self ) :
        """
        Pions for   B -> psi X lines
        """
        from GaudiConfUtils.ConfigurableGenerators import FilterDesktop
        ## from StandardParticles                     import StdAllLoosePions as inpts
        from StandardParticles                     import StdNoPIDsPions   as inpts
        ##
        return self.make_selection (
            'Pion'                 ,
            FilterDesktop          ,
            [ inpts ]              ,
            Code = self['PionCut'] ,
        )

    jpsi_name = 'FullDSTDiMuonJpsi2MuMuDetachedLine'
    psi2_name = 'FullDSTDiMuonPsi2MuMuDetachedLine'

    from PhysSelPython.Wrappers import AutomaticData
    jpsi  = AutomaticData ( '/Event/AllStreams/Phys/%s/Particles' % jpsi_name )
    psi2s = AutomaticData ( '/Event/AllStreams/Phys/%s/Particles' % psi2_name )
    #
    ## merged selectoon for J/psi & psi'
    #
    from PhysSelPython.Wrappers import MergedSelection
    psis = MergedSelection (
        'SelDetachedPsisForBandQ' ,
        RequiredSelections = [ jpsi ]
    )

    def _psi_ ( self ) :
        """
        psi(') -> mu+ mu-
        """
        return psis


    PSIX.pions = _pions_
    PSIX.kaons = _kaons_
    PSIX.psi   = _psi_

    ## use builder
    builder = PSIX ( 'PsiX' , builder_configuration  )



    from PhysSelPython.Wrappers import SelectionSequence

    psi3k      = SelectionSequence ( 'Psi3K'       , builder.psi_3K   () )
    psi3kpi    = SelectionSequence ( 'Psi3Kpi'     , builder.psi_3Kpi () )


    davinci = DaVinci(
        InputType     = 'DST'    ,
        Simulation    = True     ,
        PrintFreq     = 1000     ,
        EvtMax        = -1       ,
        Lumi          = True     ,
        DataType = params['Year'],
        DDDBtag = params['DDDB'],
        CondDBtag = params['SIMCOND'],
        # HistogramFile = 'DVHistos.root' ,
        TupleFile     = 'output_loop.root' ,
    )

    from Configurables import GaudiSequencer
    # seq   = GaudiSequencer('SEQ1', Members=[psi3k.sequence()])
    seq   = GaudiSequencer('SEQ2', Members=[psi3kpi.sequence()])


    my_name = "Bplus"


    davinci.UserAlgorithms = [ my_name ]

    setData ( datafiles , catalogs , castor )

    gaudi = appMgr()

    from StandardParticles import StdAllNoPIDsPions, StdAllNoPIDsKaons

    # create local algorithm:
    alg = MCAnalysisAlgorithm(
        my_name,
        Inputs = [
            StdAllNoPIDsPions.outputLocation(),
            StdAllNoPIDsKaons.outputLocation(),
            '/Event/AllStreams/Phys/%s/Particles' % jpsi_name
        ] ,
        PP2MCs = [ 'Relations/Rec/ProtoP/Charged' ],
        ReFitPVs = True
    )

    return SUCCESS





# =============================================================================
# job steering
if __name__ == '__main__':

    # make printout of the own documentations
    print '*' * 120
    print __doc__
    print ' Author  : %s ' % __author__
    print ' Version : %s ' % __version__
    print ' Date    : %s ' % __date__
    print '*' * 120

    inputdata = ['/lhcb/MC/2011/ALLSTREAMS.DST/00030463/0000/00030463_00000001_1.allstreams.dst']
    test_params = {'DDDB': 'Sim08-20130503', 'Mode': 'B+ -> J/psi K+ K+ K-', 'SIMCOND': 'Sim08-20130503-vc-md100', 'Year': '2011'}

    configure(inputdata, params=test_params, castor=True)

    # run the job
    run(3000)

    gaudi = appMgr()

    myalg2 = gaudi.algorithm ( 'Bplus' )
