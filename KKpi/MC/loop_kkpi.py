#!/bin/env python

#=============================================================================
__author__ = "  "
__date__ = "  "
__version__ = "  "
#=============================================================================

# import everything from BENDER
from Bender.All import *
from Gaudi.Configuration import *
import BenderTools.Fill
import BenderTools.TisTos

# needed for job configuration

# import of the useful units from Gaudi
from GaudiKernel.SystemOfUnits import GeV, MeV, mm, micrometer
from GaudiKernel.PhysicalConstants import c_light
from math import *
from LoKiTracks.decorators import *  # needed for TrKEY work


import LHCbMath.Types                 ## easy access to various geometry routines
from   Gaudi.Configuration import *   ## needed for job configuration

from   GaudiKernel.SystemOfUnits     import GeV, MeV, mm
from   GaudiKernel.PhysicalConstants import c_light

import math
from LoKiCore.basic import cpp

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

        return SUCCESS

    def analyse(self):
        prims = self.vselect('prims', ISPRIMARY)
        if prims.empty():
            return self.Warning("No primary vertices", SUCCESS)

        mcB = self.mcselect(
            'mcB', "[( B+ ==>  ( J/psi(1S) =>  mu+  mu-  )  K+  K- pi+ )]CC")

        mcB_nonres = self.mcselect(
            'mcB_nonres', "[( B+ =>  ( J/psi(1S) =>  mu+  mu-  )  K+  K- pi+ )]CC")

        mcB_Ks = self.mcselect(
            'mcB_Ks', "[( B+ ==>  ( J/psi(1S) =>  mu+  mu-  )  (K*(892)~0 => K- pi+) K+ )]CC")

        mcB_f0 = self.mcselect( 
            'mcB_f0', "[( B+ ==>  ( J/psi(1S) =>  mu+  mu-  )  (f_0(980) => K+ K- ) pi+ )]CC")

        mcB_f2 = self.mcselect(
            'mcB_f2', "[( B+ ==>  ( J/psi(1S) =>  mu+  mu-  )  (f_2(1270) => K+ K- ) pi+ )]CC")


        nt_sizes = self.nTuple("sizes")

        nt_sizes.column('mcB', mcB.size())
        nt_sizes.column('mcB_Ks', mcB_Ks.size())
        nt_sizes.column('mcB_f0', mcB_f0.size())
        nt_sizes.column('mcB_f2', mcB_f2.size())

        nt_sizes.write()


        if mcB.size() != 1 or ((mcB_nonres.size() + mcB_Ks.size() + mcB_f0.size() + mcB_f2.size()) != 1):
            return self.Warning("Something wrong with MC size " + str(mcB.size()), SUCCESS)


        mcK = self.mcselect(
            "mcK",  "[( B+ ==>  ( J/psi(1S) =>  mu+  mu-  )  ^K+  ^K- pi+ )]CC")
        mcPi = self.mcselect(
            "mcPi",  "[( B+ ==>  ( J/psi(1S) =>  mu+  mu-  )  K+  K- ^pi+ )]CC")
        mcMu = self.mcselect(
            "mcMu", "[( B+ ==>  ( J/psi(1S) =>  ^mu+  ^mu-  )  K+  K- pi+ )]CC")
        mcPsi = self.mcselect(
            "mcPsi", "[( B+ ==>  ^( J/psi(1S) =>  mu+  mu-  )  K+  K- pi+ )]CC")

        if mcK.empty() or mcMu.empty() or mcPsi.empty() or mcPi.empty():
            return self.Warning('No true MC-decay components are found', SUCCESS )


        match = self.mcTruth()
        trueB = MCTRUTH(match, mcB)
        trueB_NR = MCTRUTH(match, mcB_nonres)
        trueB_Ks = MCTRUTH(match, mcB_Ks)
        trueB_f0 = MCTRUTH(match, mcB_f0)
        trueB_f2 = MCTRUTH(match, mcB_f2)


        trueK = MCTRUTH(match, mcK)
        truePi = MCTRUTH(match, mcPi)
        truePsi = MCTRUTH(match, mcPsi)
        trueMu = MCTRUTH(match, mcMu)



        kaons  = self.select ( 'K'  ,  ('K+'  == ABSID) & trueK )
        kplus  = self.select ( 'K+' ,  kaons , Q > 0 )
        kminus = self.select ( 'K-' ,  kaons , Q < 0 )
        pions  = self.select ( 'pi'  ,  ('pi+'  == ABSID) & truePi )
        psis = self.select( "Jpsi", ('J/psi(1S)'  == ABSID) & truePsi )


        k_counter = self.counter("k_counter")
        pi_counter = self.counter("pi_counter")
        jpsi_counter = self.counter("jpsi_counter")
        b_counter = self.counter("b_counter")

        k_counter += kaons.size()
        pi_counter += pions.size()
        jpsi_counter += psis.size()

        if kaons.empty():
            return self.Warning("No reconstructed kaons", SUCCESS)  # RETURN
    
        if pions.empty():
            return self.Warning("No reconstructed pions", SUCCESS)  # RETURN
        
        if kplus.empty() or kminus.empty():
            return self.Warning("K+/- empty", SUCCESS)  # RETURN

        if psis.empty():
            return self.Warning("No reconstructed psis", SUCCESS)  # RETURN

        # myB = self.select('Bu' , '[( B+ ->  ( J/psi(1S) ->  mu+  mu-  )  K+  K+  K-)]CC' )
        # myB = self.select('myB' , '[( B+ ->  J/psi(1S)  K+  K- pi+)]CC' )
        myB = self.loop('Jpsi pi K+ K-', 'B+')

        cut_pi  = ( PT          > 200 * MeV ) & \
                in_range ( 2          , ETA , 5         ) & \
                in_range ( 3.2 * GeV  , P   , 150 * GeV ) & \
                ( MIPCHI2DV()  > 4        )

                # ( PROBNNpi     > 0.1      ) & \
                # ( CLONEDIST   > 5000      ) & \
                # ( TRGHOSTPROB < 0.5       ) & \
                # ( TRCHI2DOF   < 4         ) & \
                # HASRICH                     #& \
        
        cut_k = ( PT          > 200 * MeV ) & \
                ( CLONEDIST   > 5000      ) & \
                ( TRGHOSTPROB < 0.5       ) & \
                ( TRCHI2DOF   < 4         ) & \
                in_range ( 2          , ETA , 5         ) & \
                in_range ( 3.2 * GeV  , P   , 150 * GeV ) & \
                HASRICH                     & \
                ( PROBNNk      > 0.1      ) & \
                ( MIPCHI2DV()  > 4        )
        cut_b = VFASPF(VCHI2PDOF) < 12 #& \
                #(CTAU > (75 * micrometer))

        for b in myB:
            if not 0 < VCHI2(b) < 100:
                continue
            pi, k1, k2 = b(2), b(3), b(4)
            if Q(pi) > 0:
                b.setPID("B+")
            else:
                b.setPID("B-")
            if not trueB(b):
                continue

            # if not (cut_pi(b(3)) and cut_pi(b(4))):
            #     continue

            if not (cut_k(k1) and cut_k(k1)):
                continue
            if not cut_pi(pi):
                continue
            # if not cut_b(b):
            #     continue

            b.save("BB")

        bb = self.selected("BB")
        b_counter += bb.size()


        # Minimal impact parameter chi2
        mipFun = MIP(prims, self.geo())
        CHI2mipFun = MIPCHI2(prims, self.geo())
        MIPCHI2DVfun= MIPCHI2DV()

        DLLKpi = PIDK - PIDpi
        DLLKp = PIDK - PIDp

        # Constrains
        dtffun_ctau = DTF_CTAU(0, True)
        dtffun_chi2 = DTF_CHI2NDOF(True, "J/psi(1S)")
        dtffun_m = DTF_FUN(M, True, "J/psi(1S)")

        dtffun_m123 = DTF_FUN(MASS(1, 2, 3), True, "J/psi(1S)")
        dtffun_m23 = DTF_FUN(M23, True, "J/psi(1S)")
        dtffun_m34 = DTF_FUN(M34, True, "J/psi(1S)")
        dtffun_m234 = DTF_FUN(MASS(2, 3, 4), True, "J/psi(1S)")


        nt = self.nTuple("t")

        for myb in bb:
            if not all([myb(i) for i in xrange(0, 5)]):
                continue


            b, jpsi, k1, k2, pi = tuple(myb(i) for i in xrange(5))

            if not dtffun_m(myb) / GeV > 5.0:
                continue

            self.treatKine(nt, b, '_b')
            self.treatKine(nt, jpsi, '_jpsi')

            # add the information for Pid efficiency correction
            self.treatPions(nt, b)
            self.treatKaons(nt, b)
            self.treatMuons(nt, b)
            self.treatTracks(nt, b)

            # add DTF-applied information
            nt.column('DTFctau', dtffun_ctau(myb))
            nt.column('DTFchi2ndof', dtffun_chi2(myb))
            nt.column('DTFm_b', dtffun_m(myb) / GeV)

            nt.column('DTFm_jpsikk', dtffun_m123(myb) / GeV)
            nt.column('DTFm_kk', dtffun_m23(myb) / GeV)
            nt.column('DTFm_kpi', dtffun_m34(myb) / GeV)
            nt.column('DTFm_kkpi', dtffun_m234(myb) / GeV)

            nt.column('MIPCHI2DV_k1', MIPCHI2DVfun(k1))
            nt.column('MIPCHI2DV_k2', MIPCHI2DVfun(k2))
            nt.column('MIPCHI2DV_pi', MIPCHI2DVfun(pi))

            # add the information needed for TisTos
            self.tisTos ( jpsi  , nt  , 'psi_' ,
                          self.lines [ 'psi' ] , self.l0tistos , self.l1tistos , self.l2tistos )

            self.tisTos ( jpsi  , nt  , 'psi1_' ,
                          self.lines [ 'psi1' ] , self.l0tistos , self.l1tistos , self.l2tistos )

            self.tisTos ( jpsi  , nt  , 'psi2_' ,
                          self.lines [ 'psi2' ] , self.l0tistos , self.l1tistos , self.l2tistos )

            self.tisTos ( jpsi  , nt  , 'psi3_' ,
                          self.lines [ 'psi3' ] , self.l0tistos , self.l1tistos , self.l2tistos )


            nt.column('mcTrueB', trueB(myb))
            nt.column('mcTrueB_NR', trueB_NR(myb))
            nt.column('mcTrueB_Ks', trueB_Ks(myb))
            nt.column('mcTrueB_f0', trueB_f0(myb))
            nt.column('mcTrueB_f2', trueB_f2(myb))


            nt.column('mcTruePsi', truePsi(jpsi))
            nt.column('mcTrueK1', trueK(k1))
            nt.column('mcTrueK2', trueK(k2))
            nt.column('mcTruePi', truePi(pi))

            nt.column('mcTrueMu1', trueMu(jpsi(1)))
            nt.column('mcTrueMu2', trueMu(jpsi(2)))


            nt.write()

        return SUCCESS

    # finalize & print histos
    def finalize(self):
        self.fill_finalize()
        self.tisTos_finalize ()
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
        'PionCut'   : """
        ( PT          > 200 * MeV ) &
        ( CLONEDIST   > 5000      ) &
        ( TRGHOSTPROB < 0.5       ) &
        ( TRCHI2DOF   < 4         ) &
        in_range ( 2          , ETA , 5         ) &
        in_range ( 3.2 * GeV  , P   , 150 * GeV ) &
        HASRICH                     &
        ( PROBNNpi     > 0.1      ) &
        ( MIPCHI2DV()  > 4        )
        """
        ### ( PROBNNpi     > 0.1      ) &
        ,
        'KaonCut'   : """
        ( PT          > 200 * MeV ) &
        ( CLONEDIST   > 5000      ) &
        ( TRGHOSTPROB < 0.5       ) &
        ( TRCHI2DOF   < 4         ) &
        in_range ( 2          , ETA , 5         ) &
        in_range ( 3.2 * GeV  , P   , 150 * GeV ) &
        HASRICH                     &
        ( PROBNNk      > 0.1      ) &
        ( MIPCHI2DV()  > 4        )
        """
        ### ( PROBNNk      > 0.1      ) &
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

    from PhysConf.Filters import LoKi_Filters
    fltrs   = LoKi_Filters (
        STRIP_Code = """
            HLT_PASS_RE('Stripping.*FullDSTDiMuonJpsi2MuMuDetachedLine.*')
        """
    )
    davinci = DaVinci(
        EventPreFilters = fltrs.filters('WG'),
        InputType     = 'DST'    ,
        Simulation    = True     ,
        PrintFreq     = 1000     ,
        EvtMax        = -1       ,
        Lumi          = True     ,
        DataType = params['Year'],
        DDDBtag = params['DDDB'],
        CondDBtag = params['SIMCOND'],
        # HistogramFile = 'DVHistos.root' ,
        TupleFile     = 'DVNtuples.root' ,
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

    inputdata = ['/lhcb/MC/2012/ALLSTREAMS.DST/00030439/0000/00030439_00000024_1.allstreams.dst']
    test_params = {'DDDB': 'Sim08-20130503-1', 'Mode': 'B+ -> J/psi K+ K- pi+', 'SIMCOND': 'Sim08-20130503-1-vc-md100', 'Year': '2012'}

    configure(inputdata, params=test_params, castor=True)

    # run the job
    run(5000)

    gaudi = appMgr()

    myalg2 = gaudi.algorithm ( 'Bplus' )
