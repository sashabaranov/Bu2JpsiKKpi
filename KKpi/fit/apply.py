import ROOT,glob 
from   Ostap.PyRoUts   import * 
from   Ostap.GetLumi   import getLumi 
import Ostap.ZipShelve as     ZipShelve 
from   AnalysisPython.Logger    import getLogger
# =============================================================================
logger = getLogger( __name__ ) 


from cuts import cuts_trg
from variables import m_Bu


## Variables in dataset
variables = [
    ## B-mass 
    ( m_Bu       , lambda s : float(s.DTFm_b) ) ,
    ## 
    # ( 'm23'      , 'm_{#pi^{+} #pi^{0}}'   ,  0.200 , 2.00  , lambda s : s.m23                 ) ,
    # ( 'mpi0'     , 'm_{#gamma #gamma}'     ,  0.080 , 0.200 , lambda s : s.m_pi0               ) ,
    # ( 'cl1'      , 'CL_{#gamma1#}'         , -0.001 , 1.001 , lambda s : s.CL_photon[0]        ) ,
    # ( 'cl2'      , 'CL_{#gamma2#}'         , -0.001 , 1.001 , lambda s : s.CL_photon[1]        ) ,
    # ( 'ann_pion' , 'P_{#pi}'               , -0.001 , 1.001 , lambda s : s.ann_pion[0]         ) ,
    # ( 'pt_pi'   , 'p_{T}(#pi^{0})'        , -0.001 , 50    , lambda s :pt_pion[0]  ) ,
    ( 'pt_psi'   , 'p_{T}(J/#psi)'         , -0.001 , 50    , lambda s : min( s.pt_jpsi , 49.99) ) ,
]

#
## prepare
# 
from Ostap.PyTMVA import Reader
mlp  = Reader (
    'Bc2Rho' ,
    ## the list of variables 
    variables = [
    ( 'DTFchi2ndof'  , lambda s: s.DTFchi2ndof ) ,
    ( 'DTFctau'      , lambda s: s.DTFctau ) ,
    #
    ( 'pt_b'         , lambda s: s.pt_b ) ,
    ( 'y_b'          , lambda s: s.y_b ) ,
    ( 'c2ip_b'       , lambda s: s.c2ip_b ) ,
    ( 'vchi2_b'      , lambda s: s.vchi2_b ) ,
    #
    ( 'pt_pion'        , lambda s: s.pt_pion[0] ) ,
    ( 'eta_pion'       , lambda s: s.eta_pion[0] ) ,
    #
    ( 'pt_kaon[0]'        , lambda s: s.pt_kaon[0] ) ,
    ( 'eta_kaon[0]'       , lambda s: s.eta_kaon[0]) ,
    #
    ( 'pt_kaon[1]'        , lambda s: s.pt_kaon[1]) ,
    ( 'eta_kaon[1]'       , lambda s: s.eta_kaon[1]) ,

    ] , 
    weights_file = 'tmva_MLP_weights.xml'
)


variables += [ ('mlp' , 'mlp' , -0.05 , 1.05 , lambda s : float(mlp(s)) ) ]

local_cuts  =  cuts_trg + " && %f<=DTFm_b  && DTFm_b<=%f " % ( m_Bu.getMin() , m_Bu.getMax() ) 

from Ostap.Selectors import  SelectorWithVars
sel_Bu = SelectorWithVars(variables, local_cuts)
sel_MC = SelectorWithVars(variables, local_cuts)

from data import selection7, mc_Pythia8

selection7.chain.process ( sel_Bu )
# mc_Pythia8.chain.process ( sel_MC )

# ## sel_MC2 = SWV (
# ##     variables ,
# ##    local_cuts
# ##    )
# ## from data_Bc   import tMC2 
# ## tMC2 .process ( sel_MC2 )
