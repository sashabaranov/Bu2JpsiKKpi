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


class MCBu2JpsiKKK(AlgoMC):

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
        primaries = self.vselect( 'PVs' , ISPRIMARY )
        if primaries.empty() :
            return self.Warning('No primary vertices are found', SUCCESS )

        mcB = self.mcselect(
            'mcB', "[( Beauty ==>  ( J/psi(1S) =>  mu+  mu-  )  K+  K+  K- )]CC")

        if 0 == mcB.size():
            return SUCCESS

        mcK = self.mcselect(
            "mcK",  "[( Beauty ==>  ( J/psi(1S) =>  mu+  mu-  )  ^K+  ^K+  ^K- )]CC")
        mcMu = self.mcselect(
            "mcMu", "[( Beauty ==>  ( J/psi(1S) =>  ^mu+  ^mu-  )  K+  K+  K- )]CC")
        mcPsi = self.mcselect(
            "mcPsi", "[( Beauty ==>  ^( J/psi(1S) =>  mu+  mu-  )  K+  K+  K- )]CC")

        if mcK.empty() or mcMu.empty() or mcPsi.empty():
            return self.Warning('No true MC-decay components are found', SUCCESS )

        match = self.mcTruth()
        trueB = MCTRUTH(match, mcB)
        trueK = MCTRUTH(match, mcK)
        truePsi = MCTRUTH(match, mcPsi)
        trueMu = MCTRUTH(match, mcMu)


        # myB = self.select('Bu' , '[( Beauty ->  ( J/psi(1S) ->  mu+  mu-  )  K+  K+  K-)]CC' )
        myB = self.select('Bu' , '[( Beauty ->  J/psi(1S)  K+  K+  K-)]CC' )
        
        # Constrains
        dtffun_ctau = DTF_CTAU(0, True)
        dtffun_chi2 = DTF_CHI2NDOF(True, "J/psi(1S)")
        dtffun_m = DTF_FUN(M, True, "J/psi(1S)")


        nt = self.nTuple("t")

        for myb in myB:
            if not all([myb(i) for i in xrange(0, 5)]):
                continue


            b = myb
            jpsi = myb.child(1)
            
            # add DTF-applied information
            nt.column('DTFm_b', dtffun_m(myb) / GeV)
            nt.column('DTFctau', dtffun_ctau(myb))
            nt.column('DTFchi2ndof', dtffun_chi2(myb))

            self.treatKine(nt, b, '_b')
            self.treatKine(nt, jpsi, '_jpsi')


            # add the information for Pid efficiency correction
            # self.treatPions(nt, b)
            self.treatKaons(nt, b)
            self.treatMuons(nt, b)
            self.treatTracks(nt, b)

            # particles without misid kaon
            particles = [myb(i) for i in xrange(1, 4)]
            # particles with misid kaon
            all_particles = [myb(i) for i in xrange(1, 5)]
            kaon = myb(4)

            pion_mass = 0.1395702 * GeV  # GeV / c^2

            E_wo_misid = reduce(lambda x, y: x + y, [E(p) for p in particles])
            # GeV/c^2 + (GeV / c)^2
            E_misid = sqrt(pion_mass ** 2 + (P(kaon)) ** 2)

            total_PX = reduce(
                lambda x, y: x + y, [PX(p) for p in all_particles])
            total_PY = reduce(
                lambda x, y: x + y, [PY(p) for p in all_particles])
            total_PZ = reduce(
                lambda x, y: x + y, [PZ(p) for p in all_particles])

            total_P_sq = (total_PX) ** 2 + (total_PY) ** 2 + (total_PZ) ** 2

            misid_Bu_M = sqrt((E_wo_misid + E_misid) ** 2 - total_P_sq)


            nt.column("m_b_misid", misid_Bu_M / GeV)

            nt.column ( 'mcTrueB'   , trueB(b)          )
            nt.column ( 'mcTruePsi' , truePsi(jpsi(0)    ))
            nt.column ( 'mcTrueK1'  , trueK(myb(2))     )
            nt.column ( 'mcTrueK2'  , trueK(myb(3))    )
            nt.column ( 'mcTrueK3'  , trueK(myb(4))    )

            nt.column ( 'mcTrueMu1' , trueMu(jpsi.child(1))    )
            nt.column ( 'mcTrueMu2' , trueMu(jpsi.child(1))    )


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
            # myb.save("GoodB")

        # myBd = self.selected("GoodB")
        # if not myBd.empty():
        #     self.setFilterPassed(True)        # IMPORTANT
        # else:
        #     return self.Warning("No reconstructed Bs", SUCCESS)  # RETURN

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
    from Configurables           import EventSelector ## needed for job configuration
    
    from Configurables import MessageSvc
    msg = MessageSvc()
    msg.setError += [ 'HcalDet.Quality'   ,
                      'EcalDet.Quality'   ,
                      'MagneticFieldSvc'  ,
                      'PropertyConfigSvc' ]

    from PhysSelPython.Wrappers import AutomaticData
    jpsi_location = 'FullDSTDiMuonJpsi2MuMuDetachedLine'
    jpsi = AutomaticData(
        Location='/Event/AllStreams/Phys/%s/Particles' % jpsi_location)

    
    # =============================================================================
    from StrippingSelections.StrippingPsiXForBandQ import PsiX_BQ_Conf as PsiX

    ## 1) redefine stripping configurations
    def _psi_(self):
        """
        psi(') -> mu+ mu-
        """
        return jpsi

    PsiX.psi = _psi_

    logger.warning("Redefine PsiX .psi")

    ## 2) unify the pion& kaon  selections
    # =============================================================================
    _PionCut_ = """
    ( CLONEDIST   > 5000   ) &
    ( TRCHI2DOF   < 4      ) &
    ( TRGHOSTPROB < 0.4    ) &
    ( PT          > 200 * MeV               ) &
    in_range ( 2          , ETA , 4.9       ) &
    in_range ( 3.2 * GeV  , P   , 150 * GeV ) &
    HASRICH                  &
    ( PROBNNpi     > 0.15  ) &
    ( MIPCHI2DV()  > 9.    )
    """
    _KaonCut_ = """
    ( CLONEDIST   > 5000   ) &
    ( TRCHI2DOF   < 4      ) &
    ( TRGHOSTPROB < 0.4    ) &
    ( PT          > 200 * MeV               ) &
    in_range ( 2          , ETA , 4.9       ) &
    in_range ( 3.2 * GeV  , P   , 150 * GeV ) &
    HASRICH                  &
    ( PROBNNk      > 0.15  ) &
    ( MIPCHI2DV()  > 9.    )
    """

    from GaudiConfUtils.ConfigurableGenerators import FilterDesktop
    _alg_pi = FilterDesktop(
        ##
        Code=_PionCut_,
        ##
    )

    from PhysSelPython.Wrappers import Selection
    from StandardParticles import StdAllNoPIDsPions as input_pions
    pions = Selection(
        "SelPiForBQ",
        Algorithm=_alg_pi,
        RequiredSelections=[input_pions]
    )

    from GaudiConfUtils.ConfigurableGenerators import FilterDesktop
    _alg_k = FilterDesktop(
        ##
        Code=_KaonCut_,
        ##
    )

    from PhysSelPython.Wrappers import Selection
    from StandardParticles import StdAllNoPIDsKaons as input_kaons
    kaons = Selection(
        "SelKForBQ",
        Algorithm=_alg_k,
        RequiredSelections=[input_kaons]
    )

    def _kaons_(self):
        return kaons

    def _pions_(self):
        return pions

    #
    ## get the selections
    #

    for s in [PsiX]:
        s.pions = _pions_
        s.kaons = _kaons_

    logger.warning("Redefine PsiX.kaons          ")
    logger.warning("Redefine PsiX.kaons          ")

    psix = PsiX('PsiX', {})

    for s in [psix.psi_3K()]:
        a = s.algorithm()
        a.ParticleCombiners = {'': 'LoKi::VertexFitter:PUBLIC'}

    from PhysSelPython.Wrappers import SelectionSequence
    sel_seq = SelectionSequence('B2Psi3K', psix . psi_3K())


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
        TupleFile     = 'DVNtuples.root' ,
    )
    

    my_name = "Bplus"
    from Configurables import GaudiSequencer
    
    davinci.UserAlgorithms = [ sel_seq.sequence() , my_name ] 
    
    setData ( datafiles , catalogs , castor )
    
    gaudi = appMgr()

    print 'seq.outputLocation()= ', sel_seq.outputLocation()    # Phys/SelPsi3KPiForPsiX/Particles
    # create local algorithm:
    alg = MCBu2JpsiKKK(
        my_name,
        Inputs = [
        sel_seq.outputLocation()
        ] ,
        PP2MCs = [ 'Relations/Rec/ProtoP/Charged' ]
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

    inputdata = ['/lhcb/MC/2011/ALLSTREAMS.DST/00030465/0000/00030465_00000015_1.allstreams.dst']
    test_params = {'DDDB': 'Sim08-20130503', 'Mode': 'B+ -> J/psi K+ pi+ pi-', 'SIMCOND': 'Sim08-20130503-vc-md100', 'Year': '2011'}

    configure(inputdata, params=test_params, castor=True)

    # run the job
    run(1000)

    gaudi = appMgr()
    
    myalg2 = gaudi.algorithm ( 'Bplus' )
