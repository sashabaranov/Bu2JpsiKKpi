import ROOT
import Ostap.ZipShelve as     ZipShelve 
from   Ostap.PyRoUts   import *
from   Ostap.Utils     import timing
from   Ostap.Utils     import rooSilent as mute 
# =============================================================================
## logging
# =============================================================================
from Bender.Logger import getLogger 
logger = getLogger( __name__ )
# =============================================================================
#
## train TMVA
#
from data import mc_Pythia6, mc_Pythia8, selection7

#
m_ins = " && (DTFm_b > 5.24 && DTFm_b < 5.36)"
m_out = " && (DTFm_b > 5.6 || DTFm_b < 5.0) && DTFm_b < 6.0"
#
from Ostap.PyTMVA import Trainer

t = Trainer( [
    ##  type                 name    configuration 
    ## ( ROOT.TMVA.Types.kMLP , "MLP", "H:!V:EstimatorType=CE:VarTransform=N:NCycles=600:HiddenLayers=N+7:TestRate=5:!UseRegulator" ) 
    ( ROOT.TMVA.Types.kMLP , "MLP", "H:!V:EstimatorType=CE:VarTransform=Norm:NCycles=600:HiddenLayers=N+7:TestRate=5:!UseRegulator" ) 
    ] )

var_list = [
    #
    ( 'DTFchi2ndof'  , 'F' ) ,
    ( 'DTFctau'      , 'F' ) ,
    #
    ( 'pt_b'         , 'F' ) ,
    ( 'y_b'          , 'F' ) ,
    ( 'c2ip_b'       , 'F' ) ,
    ( 'vchi2_b'      , 'F' ) ,
    #
    ( 'pt_pion'        , 'F' ) ,
    ( 'eta_pion'       , 'F' ) ,
    #
    ( 'pt_kaon[0]'        , 'F' ) ,
    ( 'eta_kaon[0]'       , 'F' ) ,
    #
    ( 'pt_kaon[1]'        , 'F' ) ,
    ( 'eta_kaon[1]'       , 'F' ) ,
    #
    #( 'CL_photon[0]' , 'F' ) , 
    #( 'CL_photon[1]' , 'F' ) , 
    #
    ]

#
## cuts:
# 
from cuts import cuts_trg 
local_cuts = cuts_trg

t.train ( var_list,
          signal          = mc_Pythia8.chain,
          background      = selection7.chain,
          outputfile      = 'tmva.root'    ,
          signal_cuts     = local_cuts + m_ins,
          background_cuts = local_cuts + m_out ,
          spectators      = []              ) 

from Ostap.PyTMVA import tmvaGUI
tmvaGUI ( 'tmva.root' )
