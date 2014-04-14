#!/usr/bin/env python
# ==========================================================================================
## @file RealConv.py
#
#  reconstraction for gamma -> e+ e- decay
# ===========================================================================================
## import ROOT                           ## needed to produce/visualize the histograms

from   Bender.Main         import *   ## import all bender goodies
from   Bender.MainMC       import *   ## import all bender goodies/
### Dasha
###
import LHCbMath.Types                 ## easy access to various geometry routines 
from   Gaudi.Configuration import *   ## needed for job configuration

from   GaudiKernel.SystemOfUnits     import GeV, MeV, mm
from   GaudiKernel.PhysicalConstants import c_light

import math
from LoKiCore.basic import cpp

import BenderTools.TisTos

AlgoMC.decisions         = BenderTools.TisTos. decisions
AlgoMC.trgDecs           = BenderTools.TisTos. trgDecs
AlgoMC.tisTos            = BenderTools.TisTos. tisTos
AlgoMC.tisTos_initialize = BenderTools.TisTos._tisTosInit
AlgoMC.tisTos_finalize   = BenderTools.TisTos._tisTosFini

# ===========================================================================================
#from LoKiHlt.decorators import L0_DATA    
# ============================================================================================
## @class Conv
#  Simple algorithm to study XXX -> Jpsi Jpsi

def track ( part ) :

    if hasattr ( part , 'track'    ) : return         part.track()
    if hasattr ( part , 'proto'    ) : return track ( part.proto()    )
    if hasattr ( part , 'particle' ) : return track ( part.particle() )

    return None

class Jpsi_mu(AlgoMC) :
    """
    Simple algorithm to study XXX -> Jpsi Jpsi
    """
    
    def initialize ( self ) :
        """
        Initialization
        """
        sc = AlgoMC.initialize ( self )
        if sc.isFailure() : return sc
        
        self._fun_ctau = BPVLTIME ( 9 ) * c_light

        triggers = {}
        triggers [ 'J/psi'  ] = {}
        triggers [ 'B'      ] = {}

        lines            = {}
        lines ['J/psi' ] = {}
        lines ['B'     ] = {}

        lines['J/psi'][  'L0TOS'] = 'L0(Muon|DiMuon)Decision'
        lines['J/psi'][  'L0TIS'] = 'L0(Muon|DiMuon|Electron|Hadron|Photon)Decision'

        lines['J/psi']['Hlt1TOS'] = 'Hlt1(DiMuon|SingleMuon|MuTrack).*Decision'
        lines['J/psi']['Hlt1TIS'] = 'Hlt1(DiMuon|SingleMuon|MuTrack|DiHadron|TrackAllL0).*Decision'

        lines['J/psi']['Hlt2TOS'] = 'Hlt2(DiMuon|SingleMuon|MuTrack).*Decision'
        lines['J/psi']['Hlt2TIS'] =  lines['J/psi']['Hlt2TOS']

        lines ['B'] = lines[ 'J/psi' ]

        sc = self.tisTos_initialize ( triggers , lines )
        
        return SUCCESS
    
    ## the only one esential method: 
    def analyse  (self ) :
        """
        The only one essential method
        """
        minDLLmu = MINTREE( 'mu+' == ABSID , PIDmu )
        minDLLk  = MINTREE( 'K+'  == ABSID , PIDK  - PIDpi )
        minDLLpi = MINTREE( 'pi+' == ABSID , PIDpi - PIDK  )
        minProbNNpion = MINTREE ('pi+' == ABSID , PROBNNpi)
        minProbNNmuon = MINTREE ('mu+' == ABSID , PROBNNmu)
        
        minCloneDist = MINTREE ( ISBASIC & HASTRACK , CLONEDIST   )
        maxGhostProb = MAXTREE ( ISBASIC & HASTRACK , TRGHOSTPROB )        
        maxTrChi2    = MAXTREE( ISBASIC  , TRCHI2DOF)

        primaries = self.vselect( 'PVs' , ISPRIMARY )
        if primaries.empty() :
            return self.Warning('No primary vertices are found', SUCCESS )

        mips = MIPCHI2( primaries, self.geo() )
        maxMips  = MINTREE( ISBASIC  , mips )
        
        dtf_chi2  = DTF_CHI2NDOF (        True , 'J/psi(1S)' )
        dtf_M     = DTF_FUN      ( M    , True , 'J/psi(1S)' )
        dtf_M1    = DTF_FUN      ( M1   , True , 'J/psi(1S)' )
        
        mass_134 = MASS(1,3,4)
        mass_234 = MASS(2,3,4)
        mass_123 = MASS(1,2,3)
        mass_23 = MASS(2,3)
        mass_24 = MASS(2,4)
        mass_34 = MASS(3,4)

        tup     = self.nTuple ( 'Bp' )
        tup2 = self.nTuple ( 'mcB' )

        bmesons = self.mcselect ( 'bc'  , ' [ Beauty ==>  J/psi(1S) pi+ ]CC' )
        
        if bmesons.empty() :
            return self.Warning('NO MC-decays', SUCCESS )

        for b in bmesons:
            tup2.column("size", bmesons.size() )
            tup2.column("pt" , MCPT  ( b ) )
            tup2.column("eta", MCETA ( b ) )
            tup2.column("y"  , MCY   ( b ) )
            tup2.column("phi", MCPHI ( b ) )
            tup2.write()

        if 1 != len( bmesons ) :
            return self.Warning('illegal number of MC-decauys', SUCCESS )

        mcBc      = self.mcselect ( 'mcBc'  , ' [ Beauty ==>  J/psi(1S) pi+ ]CC' )
        mcMu      = self.mcselect ( 'mcMu'  , ' [ Beauty ==> (J/psi(1S) => ^mu+ ^mu- )   pi+ ]CC' )
        mcPsi     = self.mcselect ( 'mcPsi' , ' [ Beauty ==> ^(J/psi(1S) =>  mu+  mu- )  pi+ ]CC' )
        mcPi      = self.mcselect ( 'mcPi'  , ' [ Beauty ==>  (J/psi(1S) =>  mu+  mu- ) ^pi+ ]CC' )

        if mcMu.empty() or mcPsi.empty() or mcPi.empty() or mcBc.empty():
            return self.Warning('No true MC-decay components are found', SUCCESS )

        trueB   = MCTRUTH ( self.mcTruth() , mcBc  )
        trueMu  = MCTRUTH ( self.mcTruth() , mcMu  )
        truePsi = MCTRUTH ( self.mcTruth() , mcPsi )
        truePi  = MCTRUTH ( self.mcTruth() , mcPi  )

        bmes   = self.select ( 'Bs' , '[Beauty -> J/psi(1S) pi+]CC' )

        for B in bmes:

            mass = M ( B )
            if not 5. * GeV < mass < 7. * GeV : continue

            psi  = B.child ( 1 )
            mpsi = M ( psi )
            if not 3.0 * GeV < mpsi < 3.2 * GeV : continue

            if minCloneDist ( B ) < 5000 : continue
            if maxGhostProb ( B ) > 0.5  : continue
            if maxTrChi2    ( B ) > 4    : continue

            vchi2 = VCHI2 ( B )
            if not 0 <= vchi2 < 500 : continue

            ctau = self._fun_ctau ( B )
            if not 50 * micrometer <= ctau < 2 * millimeter : continue

            pion_1 = B.child(2)
            muon_1 = psi.child(1)
            muon_2 = psi.child(2)

            pt_b    = PT ( B ) / GeV
            pt_psi  = PT ( psi ) / GeV
            pt_m1   = PT ( muon_1 ) / GeV
            pt_m2   = PT ( muon_2 ) / GeV
            pt_pi   = PT ( pion_1 ) / GeV

            tup.column ( 'vchi2'    ,     vchi2            )

            tup.column ( 'mjkk'     ,     M    ( B ) / GeV )
            tup.column ( 'mmjp'     ,     M1   ( B ) / GeV )
            tup.column ( 'm1jkk'    , dtf_M    ( B ) / GeV )
            tup.column ( 'mcjp'     , dtf_M1   ( B ) / GeV )
            tup.column ( 'chi2dtf2' , dtf_chi2 ( B )       )
            tup.column ( 'ctau'     , ctau                 )

            tup.column ( 'ptb'  , pt_b     )
            tup.column ( 'ptpsi', pt_psi   )
            tup.column ( 'ptm1' , pt_m1    )
            tup.column ( 'ptm2' , pt_m2    )
            tup.column ( 'ptpi' , pt_pi    )

            tup.column ( "dllm"    , minDLLmu( B ) )
            tup.column ( "dllp"    , minDLLpi( B ) )
            tup.column ( 'nnpi'    , minProbNNpion( B ) )
            tup.column ( 'nnmu'    , minProbNNmuon( B ) )
            tup.column ( "maxMips" , maxMips ( B ) )

            tup.column ( 'mcTrueB'     , trueB(B)          )
            tup.column ( 'mcTrueJ'     , truePsi(psi)      )
            tup.column ( 'mcTruePi'    , truePi(pion_1)    )
            tup.column ( 'mcTrueMu1'   , trueMu(muon_1)    )
            tup.column ( 'mcTrueMu2'   , trueMu(muon_2)    )

            tup.column ( "eta" , ETA ( B ) )
            tup.column ( "phi" , PHI ( B ) )
            tup.column ( "y"   , Y   ( B ) )

            self.tisTos ( psi , tup  , 'psi_' , self.lines['B' ] )

            tup.write  ()

        return SUCCESS
        
    def finalize ( self ) :

        self._fun_ctau = None 
        return AlgoMC.finalize ( self )

        
# =============================================================================
## configure the job 

## configure the job
def configure ( datafiles , catalogs = [] , castor = False ) :
    """
    Job configuration
    """

    from Configurables           import DaVinci       ## needed for job configuration
    from Configurables           import EventSelector ## needed for job configuration
    
    from Configurables import MessageSvc
    msg = MessageSvc()
    msg.setError += [ 'HcalDet.Quality'   ,
                      'EcalDet.Quality'   ,
                      'MagneticFieldSvc'  ,
                      'PropertyConfigSvc' ]


##     # =========================================================================
##     ## 0) Rerun stripping if MC is nt MC/2011 or MC/2012 
##     # =========================================================================
##     if   '2012' == the_year : 
##         import StrippingArchive.Stripping20.StrippingDiMuonNew              as DiMuon 
##         import StrippingSettings.Stripping20.LineConfigDictionaries_BandQ   as LineSettings
##     elif '2011' == the_year :
##         import StrippingArchive.Stripping20r1.StrippingDiMuonNew            as DiMuon 
##         import StrippingSettings.Stripping20r1.LineConfigDictionaries_BandQ as LineSettings

##     config  = LineSettings.FullDSTDiMuon['CONFIG']
##     name    = 'FullDST'
##     builder = DiMuon.DiMuonConf ( name , config )

##     ## selection
##     jpsi  = builder.SelJpsi2MuMuDetached


    # =========================================================================
    ## 0) Otherwise use existing stripping ilne 
    # =========================================================================
    
    from PhysSelPython.Wrappers import AutomaticData
    jpsi_location = 'FullDSTDiMuonJpsi2MuMuDetachedLine'
    jpsi = AutomaticData ( Location = '/Event/AllStreams/Phys/%s/Particles' % jpsi_location )

    
    # =============================================================================
    from StrippingSelections.StrippingPsiXForBandQ import PsiX_BQ_Conf    as PsiX
    
    # =============================================================================
    ## 1) redefine stripping configurations 
    # ============================================================================= 
    
    #
    ## redefine psi(') -> mu+ mu-
    # 
    def _psi_ ( self ) :
        """
        psi(') -> mu+ mu- 
        """
        return jpsi
    
    PsiX  . psi = _psi_
    
    logger.warning ( "Redefine PsiX .psi" )
    
    # =============================================================================
    ## 2) unify the pion& kaon  selections 
    # =============================================================================
    _PionCut_  = """
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
    _KaonCut_  = """
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
    _alg_pi  = FilterDesktop (
        ##
        Code = _PionCut_ ,
        ##
        )
    
    from PhysSelPython.Wrappers  import Selection
    from StandardParticles       import StdAllNoPIDsPions as input_pions 
    pions  = Selection (
        "SelPiForBQ"                       ,
        Algorithm          =  _alg_pi      ,
        RequiredSelections = [ input_pions ]  
        )
    
    from GaudiConfUtils.ConfigurableGenerators import FilterDesktop
    _alg_k  = FilterDesktop (
        ##
        Code = _KaonCut_ ,
        ##
        )
    
    from PhysSelPython.Wrappers  import Selection
    from StandardParticles       import StdAllNoPIDsKaons as input_kaons 
    kaons  = Selection (
        "SelKForBQ"      ,
        Algorithm          =   _alg_k      ,
        RequiredSelections = [ input_kaons ]  
        )
    
    
    def _kaons_   ( self ) : return kaons 
    def _pions_   ( self ) : return pions 
    
    #
    ## get the selections 
    #
    
    for s in [ PsiX ]  :
        s.pions   = _pions_
        s.kaons   = _kaons_
        
    logger.warning ( "Redefine PsiX.kaons          " )
    logger.warning ( "Redefine PsiX.kaons          " )

        
    psix   = PsiX   ( 'PsiX'  , {} )


    for s in [ psix.psi_pi() ] :
        
        a = s.algorithm ()
        a.ParticleCombiners = { '' : 'LoKi::VertexFitter:PUBLIC' } 
        
        
    from PhysSelPython.Wrappers import      SelectionSequence    
    sel_seq = SelectionSequence ( 'B2PsiPi'   , psix . psi_pi   () )
    

    the_year = '2012'
    davinci = DaVinci (
        DataType      = the_year ,
        InputType     = 'DST'    ,
        Simulation    = True     ,
        PrintFreq     = 1000     ,
        EvtMax        = -1       , 
        HistogramFile = 'DVHistos.root' ,
        TupleFile     = 'DVNtuples.root' ,
        Lumi          = True ,
        ##
        # MC : 
        ## SIMCOND : 'Sim08-20130503-1', 'Sim08-20130503-1-vc-md100'
        #
        DDDBtag   = "Sim08-20130503-1"         ,
        CondDBtag = "Sim08-20130503-1-vc-md100"    
        )
    
    
    my_name = "Bplus"
    from Configurables import GaudiSequencer
    
    davinci.UserAlgorithms = [ sel_seq.sequence() , my_name ] 
    
    setData ( datafiles , catalogs , castor )
    
    gaudi = appMgr()
    
    print 'seq.outputLocation()= ', sel_seq.outputLocation()    # Phys/SelPsi3KPiForPsiX/Particles
    alg = Jpsi_mu(
        my_name               ,   ## Algorithm name
        Inputs = [
        sel_seq.outputLocation()
        ] ,
        PP2MCs = [ 'Relations/Rec/ProtoP/Charged' ]
        )

    return SUCCESS 

# =============================================================================
# The actual job steering
if '__main__' == __name__ :
    
    data1_3 = [
        '/lhcb/MC/2012/ALLSTREAMS.DST/00022044/0000/00022044_00000002_1.allstreams.dst'
    ]        
    
    files = data1_3
    
    configure ( files  , castor = True )

    run ( 1000 )
    
    gaudi = appMgr()
    
    myalg2 = gaudi.algorithm ( 'Bplus' )

# =============================================================================
# The END 
# =============================================================================
