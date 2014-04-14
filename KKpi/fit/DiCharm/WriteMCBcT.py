#!/usr/bin/env ipython
# ========================================================================
# @file DiCharm/BcTMC.py
#
#  Simple algorithm for MC-Bc
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

Simple algorithm for Bc

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
# =============================================================================
from GaudiKernel.SystemOfUnits import GeV, MeV, mm, meter
from GaudiKernel.PhysicalConstants import c_light
# =============================================================================
# logging
# =============================================================================
from Bender.Logger import getLogger
logger = getLogger(__name__)
# =============================================================================

# configure the job


def configure(datafiles,
              catalogs=[],
              castor=False,
              params={}):
    """
    Job configuration 
    """

    logger.info("start: params: %s " % params)
    from BenderTools.Parser import hasInFile

    the_year = "2011"

    if params:
        the_year = params['Year']
        logger.info('Year is set from params to be %s ' % the_year)
    else:
        if hasInFile(datafiles, 'Collision11'):
            the_year = '2011'
        elif hasInFile(datafiles, 'Collision12'):
            the_year = '2012'
        elif hasInFile(datafiles, 'Stripping17'):
            the_year = '2011'
        elif hasInFile(datafiles, 'Stripping13'):
            the_year = '2011'
        elif hasInFile(datafiles, 'Stripping15'):
            the_year = '2011'
        elif hasInFile(datafiles, 'Stripping17'):
            the_year = '2011'
        elif hasInFile(datafiles, 'Stripping19'):
            the_year = '2012'
        logger.info('Year is set from files  to be %s ' % the_year)

    #
    # check
    #
    if '2011' == the_year and hasInFile(datafiles, 'Collision12'):
        raise AttributeError, 'Invalid Year %s ' % the_year
    if '2012' == the_year and hasInFile(datafiles, 'Collision11'):
        raise AttributeError, 'Invalid Year %s ' % the_year

    logger.info('Use the Year = %s ' % the_year)

    #
    # dimuon locations in DIMUON.DST
    #
    Jpsi_det_location = '/Event/AllStreams/Phys/FullDSTDiMuonJpsi2MuMuDetachedLine/Particles'
    Jpsi_unb_location = '/Event/AllStreams/Phys/FullDSTDiMuonJpsi2MuMuTOSLine/Particles'

    if hasInFile(datafiles, '/MC/MC11a/'):
        Jpsi_unb_location = '/Event/AllStreams/Phys/FullDSTDiMuonJpsi2MuMuLine/Particles'
        logger.warning('Rename ``unbised line'' to be %s ' %
                       Jpsi_unb_location)

    from PhysSelPython.Wrappers import AutomaticData
    jpsi_det = AutomaticData(Location=Jpsi_det_location)
    jpsi_unb = AutomaticData(Location=Jpsi_unb_location)

    # pions :
    from GaudiConfUtils.ConfigurableGenerators import FilterDesktop
    alg_pions = FilterDesktop(
        #
        Preambulo=["from LoKiPhysMC.decorators import mcMatch"],
        Code="""
        ( PT > 500  * MeV      ) & 
        ( CLONEDIST   > 5000   ) & 
        ( TRGHOSTPROB < 0.5    ) &
        ( TRCHI2DOF   < 4      ) & 
        in_range ( 2 , ETA , 5 ) &
        HASRICH                  
        & mcMatch( '[B_c+  =>  ( J/psi(1S) => mu+ mu- ) ^pi+ ]CC' , 1 )
        """ ,
    )

    ## ( PROBNNpi > 0.10      )

    from PhysSelPython.Wrappers import Selection
    from StandardParticles import StdAllLoosePions
    pions = Selection(
        "ThePions",
        Algorithm=alg_pions,
        RequiredSelections=[StdAllLoosePions]
    )

    alg_kaons = FilterDesktop(
        #
        Preambulo=["from LoKiPhysMC.decorators import mcMatch"],
        Code="""
        ( PT > 500  * MeV      ) & 
        ( CLONEDIST   > 5000   ) & 
        ( TRCHI2DOF   < 4      ) & 
        ( TRGHOSTPROB < 0.5    ) & 
        in_range ( 2 , ETA , 5 ) &
        HASRICH                  
        & mcMatch( '[B+  =>  ( J/psi(1S) => mu+ mu- ) ^K+ ]CC' , 1 )
        """
    )

    ## ( PROBNNk > 0.10       )

    from StandardParticles import StdAllLooseKaons
    kaons = Selection(
        "TheKaons",
        Algorithm=alg_kaons,
        RequiredSelections=[StdAllLooseKaons]
    )

    #
    Preambulo = [
        # shortcut for chi2 of vertex fit
        "chi2vx   = VFASPF(VCHI2)",
        # shortcut for the c*tau
        "from GaudiKernel.PhysicalConstants import c_light",
        # use the embedded cut for chi2(LifetimeFit)<16
        "ctau_25   = BPVLTIME ( 25 ) * c_light ",
        "mbc_acut  = in_range ( 6.050 * GeV , AM , 6.550 * GeV ) ",
        "mbp_acut  = in_range ( 5.100 * GeV , AM , 5.550 * GeV ) ",
        # mass-cut for beauty particles
        "mbc_cut   = in_range ( 6.100 * GeV ,  M , 6.500 * GeV ) ",
        "mbp_cut   = in_range ( 5.150 * GeV ,  M , 5.500 * GeV ) ",
        #
        "from LoKiPhysMC.decorators import mcMatch"
    ]

    # Bc :
    from GaudiConfUtils.ConfigurableGenerators import CombineParticles
    alg_bc = CombineParticles(
        DecayDescriptor="[B_c+ -> J/psi(1S) pi+ ]cc",
        DaughtersCuts={
            "J/psi(1S)": "mcMatch( '[ B_c+  =>  ^( J/psi(1S) => mu+ mu- ) pi+ ]CC' , 1 )"
        },
        Preambulo=Preambulo,
        CombinationCut=" mbc_acut ",
        MotherCut="""
        ( chi2vx < 25  ) & mbc_cut &  ( ctau_25 > %s )
        & mcMatch( '[ B_c+  => ( J/psi(1S) => mu+ mu- ) pi+ ]CC' , 1 )
        """ %  ( 40 * micrometer )
    )

    # B+ :
    alg_bu = CombineParticles(
        DecayDescriptor="[B+ -> J/psi(1S) K+ ]cc",
        DaughtersCuts={
            "J/psi(1S)": "mcMatch( '[ B+  => ^( J/psi(1S) => mu+ mu- ) K+ ]CC' , 1 )"
        },
        Preambulo=Preambulo,
        CombinationCut=" mbp_acut ",
        MotherCut="""
        ( chi2vx < 25  ) & mbp_cut &  ( ctau_25 > %s )
        & mcMatch( '[ B+  =>  ( J/psi(1S) => mu+ mu- ) K+ ]CC' , 1 )
        """ %  ( 40 * micrometer )
    )

    sel_bc_det = Selection(
        "TheBc_det",
        Algorithm=alg_bc,
        RequiredSelections=[jpsi_det, pions]
    )

    sel_bc_unb = Selection(
        "TheBc_unb",
        Algorithm=alg_bc,
        RequiredSelections=[jpsi_unb, pions]
    )

    sel_bu_det = Selection(
        "TheB_det",
        Algorithm=alg_bu,
        RequiredSelections=[jpsi_det, kaons]
    )

    sel_bu_unb = Selection(
        "TheB_unb",
        Algorithm=alg_bu,
        RequiredSelections=[jpsi_unb, kaons]
    )

    from PhysSelPython.Wrappers import SelectionSequence
    Bc_det_seq = SelectionSequence("Bc_det", TopSelection=sel_bc_det)
    Bu_det_seq = SelectionSequence("Bu_det", TopSelection=sel_bu_det)
    Bc_unb_seq = SelectionSequence("Bc_unb", TopSelection=sel_bc_unb)
    Bu_unb_seq = SelectionSequence("Bu_unb", TopSelection=sel_bu_unb)

    #
    # selection sequence
    #

    mode = params['Mode']
    mode = mode.upper()
    if 0 <= mode.find('BC') or 0 <= mode.find('PI'):
        seqs = [Bc_det_seq, Bc_unb_seq]
        # seqs = [ Bc_det_seq ]
    elif 0 <= mode.find('B+') or 0 <= mode.find('K+'):
        seqs = [Bu_det_seq, Bu_unb_seq]
        # seqs = [ Bu_det_seq ]
    else:
        raise

    from PhysSelPython.Wrappers import MultiSelectionSequence
    B_SEQ = MultiSelectionSequence(
        "B2PSI",
        Sequences=seqs,
    )

    from DSTWriters.Configuration import (SelDSTWriter,
                                          stripMicroDSTStreamConf,
                                          stripMicroDSTElements)

    # Configuration of SelDSTWriter
    SelDSTWriterConf = {
        'default': stripMicroDSTStreamConf(pack=False)}
    SelDSTWriterElements = {
        'default': stripMicroDSTElements(pack=False)}

    udstWriter = SelDSTWriter(
        "MyMicroDSTWriter",
        StreamConf=SelDSTWriterConf,
        MicroDSTElements=SelDSTWriterElements,
        OutputFileSuffix='MCBSEQ',
        SelectionSequences=[B_SEQ]
    )

    # Read only fired events to speed up
    from PhysConf.Filters import LoKi_Filters
    fltrs = LoKi_Filters(
        STRIP_Code="HLT_PASS_RE('Stripping.*DiMuonJpsi2MuMu.*')",
        VOID_Code="""
        ( 0.5 < CONTAINS('%s') ) | 
        ( 0.5 < CONTAINS('%s') ) 
        """ %  ( Jpsi_det_location , Jpsi_unb_location )
    )
    #
    fltrs_0 = LoKi_Filters(
        VOID_Code="""
        ( EXISTS ( '/Event/DAQ/RawEvent') | EXISTS('/Event/Trigger/RawEvent' ) ) 
        & EXISTS ( '/Event/Strip/Phys/DecReports')
        & ( RECSUMMARY (  0 , -1 ) >  0.5 )
        """
    )

    #
    # finally: DaVinci
    #
    ## needed for job configuration
    from Configurables import DaVinci
    davinci = DaVinci(
        EventPreFilters=fltrs_0.filters(
            'Filters0') + fltrs.filters('Filters1'),
        DataType=the_year,
        InputType='MDST',
        Simulation=True,
        PrintFreq=1000,
        EvtMax=-1,
        #
        HistogramFile='MCBc_Histos.root',
        TupleFile='MCBc.root',
        #
        DDDBtag=params['DDDB'],
        CondDBtag=params['SIMCOND'],
        #
        Lumi=False  # True ,
        #
    )

    davinci.appendToMainSequence([udstWriter.sequence()])

    #
    # come back to Bender
    #
    setData(datafiles, catalogs, castor)

    from BenderTools.Utils import silence
    silence()

    #
    gaudi = appMgr()
    #

    logger.info("end: params: %s " % params)

    return SUCCESS

# =============================================================================
# The actual job steering
if '__main__' == __name__:

    #
    # Bc -> J/psi pi
    #
    # path = '/MC/MC11a/Beam3500GeV-2011-MagDown-Nu2-EmNoCuts/Sim05/Trig0x40760037Flagged/Reco12a/Stripping17NoPrescalingFlagged/42122001/ALLSTREAMS.DST'
    inputs_Bc = [
        '/lhcb/MC/MC11a/ALLSTREAMS.DST/00014472/0000/00014472_00000006_1.allstreams.dst',
        '/lhcb/MC/MC11a/ALLSTREAMS.DST/00014472/0000/00014472_00000015_1.allstreams.dst',
        '/lhcb/MC/MC11a/ALLSTREAMS.DST/00014472/0000/00014472_00000053_1.allstreams.dst',
        '/lhcb/MC/MC11a/ALLSTREAMS.DST/00014472/0000/00014472_00000011_1.allstreams.dst'
    ]
    params_Bc = {
        'Mode': 'Bc -> J/psi pi',
        'Year': '2011',
        'DDDB': 'MC11-20111102',
        'SIMCOND': 'sim-20111111-vc-md100'
    }

    #
    # B+ -> J/psi K
    #
    # path = '/MC/MC11a/Beam3500GeV-2011-MagDown-Nu2-EmNoCuts/Sim05/Trig0x40760037Flagged/Reco12/Stripping17NoPrescalingFlagged/12143001/ALLSTREAMS.DST'
    inputs_Bu = [
        '/lhcb/MC/MC11a/ALLSTREAMS.DST/00013192/0000/00013192_00000001_1.allstreams.dst',
        '/lhcb/MC/MC11a/ALLSTREAMS.DST/00013416/0000/00013416_00000022_1.allstreams.dst',
        '/lhcb/MC/MC11a/ALLSTREAMS.DST/00013192/0000/00013192_00000010_1.allstreams.dst',
        '/lhcb/MC/MC11a/ALLSTREAMS.DST/00013416/0000/00013416_00000017_1.allstreams.dst',
        '/lhcb/MC/MC11a/ALLSTREAMS.DST/00013416/0000/00013416_00000031_1.allstreams.dst'
    ]
    params_Bu = {
        'Mode': 'B+ -> J/psi K+',
        'Year': '2011',
        'DDDB': 'MC11-20111102',
        'SIMCOND': 'sim-20111111-vc-md100'
    }
    #

    configure(
        inputs_Bu,
        castor=True,
        params=params_Bu
    )

    run(1000)

    gaudi = appMgr()


# =============================================================================
# The END
# =============================================================================
