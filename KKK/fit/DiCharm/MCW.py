#!/usr/bin/env python
# ========================================================================
# @file DiCharm/MCCW.py
#
#  Helper script to get MC for W+/-
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

Helper script to get MC for W+/-

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
__date__ = "2013-05-07"
__version__ = '$Revision$'
# =============================================================================
## needed to produce/visualize the histograms
import ROOT
from Bender.MainMC import *  # import all bender goodies
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
import BenderTools.TisTosMC  # add methods for TisTos
import BenderTools.Fill  # add methods to fill n-tuple info
# =============================================================================

# =============================================================================
# @class MCW
#  Helper script to get fakes for open charm + W analysis
#  @date   2011-05-27
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class MCW (AlgoMC):

    """
    Simple algorithm 
    """

    def initialize(self):

        sc = AlgoMC.       initialize(self)
        if sc.isFailure():
            return sc

        triggers = {}
        triggers['Z0'] = {}
        triggers['W+'] = {}
        triggers['W-'] = {}
        from DiCharm.TisTos import lines

        sc = self.tisTos_initialize(triggers, lines)
        if sc.isFailure():
            return sc

        sc = self.  fill_initialize()
        if sc.isFailure():
            return sc

        #

        self.mfit = DTF_FUN(M, True)
        self.c2dtf = DTF_CHI2NDOF(True)
        self.ip = BPVIP()
        self.perr2 = PERR2 / P2

        ## self.ptCone_ = SUMCONE (   0.25 , PT , '/Event/Phys/StdAllLoosePions/Particles'   )
        ## self.etCone_ = SUMCONE (   0.25 , PT , '/Event/Phys/StdLooseAllPhotons/Particles' )

        self.ptCone_ = PINFO(55001, -100 * GeV)
        self.etCone_ = PINFO(55002, -100 * GeV)
        self.ptCone_2 = PINFO(55003, -100 * GeV)

        # print 'mu1/PT:' , self.ptCone_ ( mu1 ) , PINFO   ( 55001 , -100 * GeV )  ( mu1 )
        # print 'mu1/ET:' , self.etCone_ ( mu1 ) , PINFO   ( 55002 , -100 * GeV
        # )  ( mu1 )

        return SUCCESS

    # finalize the algorithm
    def finalize(self):
        """
        Finalize the action 
        """
        self.mfit = None
        self.c2dtf = None
        self.ip = None
        self.perr2 = None

        self.ptCone_ = None
        self.ptCone_2 = None
        self.etCone_ = None

        self.  fill_finalize()
        self.tisTos_finalize()

        return Algo.finalize(self)

    # the only one esential method:
    def analyse(self):
        """
        The only one essential method
        """

        mcmu = self.mcselect("mcmu", "[W+ => ^mu+ Nu]CC")
        if not mcmu:
            return self.Warning('No muons from W are found', SUCCESS)

        mcmufast = self.mcselect('goodmcmu', mcmu, MCPT > 10 * GeV)
        if not mcmu:
            return self.Warning('No fast muons from W are found', SUCCESS)
        if 1 < len(mcmu):
            self.Warning('too many fast muons are found', SUCCESS)

        primaries = self.vselect('PVs', ISPRIMARY)
        if primaries.empty():
            return self.Warning('No primary vertices are found', SUCCESS)

        #
        # rec summary
        rc_summary = self.get_(
            '/Event/Rec/Summary').summaryData()
        ew_counter = self.get_('/Event/Counters/CharmEW', False)
        #

        W = self.select('mu', 'mu+' == ABSID)

        mcMu = MCTRUTH(self.mcTruth(), mcmu)

        tupW = self.nTuple('W')

        for w in W:

            mu1 = w

            self.treatKine(tupW, mu1, '_mu1')

            tupW.column_float('c2dtf_mu1', self.c2dtf(mu1))
            tupW.column_float('ip_mu1', self.ip(mu1))
            tupW.column_bool('ismuon_mu1', ISMUON(mu1))
            tupW.column_bool('inmuon_mu1', INMUON(mu1))
            tupW.column_bool('isloose_mu1', ISMUONLOOSE(mu1))

            # mc-matched ?
            tupW.column_bool('mc', mcMu(mu1))

            pv = self.bestVertex(w)
            if not pv:
                self.Warning('Illegal primary vertex for W!')
                continue

            ip2 = IP(self.geo(), pv)

            tupW.column_float('ip2_mu1', ip2(mu1))
            tupW.column_float('perr2_mu1', self.perr2(mu1))
            tupW.column_float('ptCone1_mu1', self.ptCone_(mu1))
            tupW.column_float('ptCone2_mu1', self.ptCone_2(mu1))
            tupW.column_float('etCone_mu1', self.etCone_(mu1))

            # add the information needed for TisTos
            self.tisTos(mu1, tupW, 'W1_',
                        self.lines['W+'],
                        self.l0tistos,
                        self.tistos)

            # add the information for Pid efficiency correction
            self.treatMuons(tupW,  w)
            # add the infomation needed for track efficiency correction
            self.treatTracks(tupW,  w)

            # add some reco-summary information
            self.addRecSummary(tupW, rc_summary)
            # more counters
            tupW.column_aux(ew_counter)

            #
            tupW.write()

        return SUCCESS

# =============================================================================
# configure the job


def configure(datafiles,
              catalogs=[],
              castor=False,
              params={}):
    """
    Job configuration 
    """

    ## needed for job configuration
    from Configurables import DaVinci

    the_year = "2011"

    from BenderTools.Parser import hasInFile

    if params:
        the_year = params['Year']
        logger.info('Year is set from params to be %s ' % the_year)
    else:
        if hasInFile(datafiles, 'Collision11'):
            the_year = '2011'
        elif hasInFile(datafiles, 'Collision12'):
            the_year = '2012'
        elif hasInFile(datafiles, 'Collision13'):
            the_year = '2013'
        elif hasInFile(datafiles, 'Stripping17'):
            the_year = '2011'
        elif hasInFile(datafiles, 'Stripping13'):
            the_year = '2011'
        elif hasInFile(datafiles, 'Stripping15'):
            the_year = '2011'
        elif hasInFile(datafiles, 'Stripping19'):
            the_year = '2012'
        elif hasInFile(datafiles, 'Stripping20r1'):
            the_year = '2011'
        elif hasInFile(datafiles, 'Stripping20r1p1'):
            the_year = '2011'
        elif hasInFile(datafiles, 'Stripping20r0p1'):
            the_year = '2012'
        elif hasInFile(datafiles, 'MC11'):
            the_year = '2011'
        logger.info('Year is set from files  to be %s ' % the_year)

    #
    # check
    #
    if '2011' == the_year and hasInFile(datafiles, 'Collision12'):
        raise AttributeError, 'Invalid Year %s ' % the_year
    if '2012' == the_year and hasInFile(datafiles, 'Collision11'):
        raise AttributeError, 'Invalid Year %s ' % the_year

    logger.info('Use the Year = %s ' % the_year)

    W_Location = '/Event/AllStreams/Phys/WMuLine/Particles'
    from PhysSelPython.Wrappers import AutomaticData
    W_Strip = AutomaticData(Location=W_Location)

    EW_preambulo = [
        "pion_cuts  = in_range ( 300 * MeV , PT , 10 * GeV ) & ( CLONEDIST > 5000 ) & ( TRCHI2DOF < 5 ) & ( TRGHOSTPROB < 0.5 ) & ( PERR2/P2 < 0.05**2 ) ",
        "ptCone_    =  SUMCONE (   0.25 , PT , '/Event/Phys/StdAllLoosePions/Particles'               )",
        "ptCone_2   =  SUMCONE (   0.25 , PT , '/Event/Phys/StdAllLoosePions/Particles'   , pion_cuts )",
        "etCone_    =  SUMCONE (   0.25 , PT , '/Event/Phys/StdLooseAllPhotons/Particles'             )",
        "ptCone     =    SINFO (  55001 , ptCone_  , True ) ",
        "ptCone2    =    SINFO (  55003 , ptCone_2 , True ) ",
        "etCone     =    SINFO (  55002 , etCone_  , True ) ",
    ]

    # ========================================================================
    # good W
    # ========================================================================
    from GaudiConfUtils.ConfigurableGenerators import FilterDesktop
    gW = FilterDesktop(
        Preambulo=EW_preambulo,
        Code="""
        in_range ( 15 * GeV , PT , 100 * GeV ) &
        ( -1e+10 * GeV < ptCone  ) &
        ( -1e+10 * GeV < ptCone2 ) &
        ( -1e+10 * GeV < etCone  ) 
        """
    )
    from PhysSelPython.Wrappers import Selection
    W_Data = Selection(
        'W',
        Algorithm=gW,
        RequiredSelections=[W_Strip]
    )

    from PhysSelPython.Wrappers import SelectionSequence
    seq = SelectionSequence("Wseq", TopSelection=W_Data)

    # counters
    from Configurables import LoKi__CounterAlg as CounterAlg
    cnt = CounterAlg(
        'CharmEWCounters',
        Location="Counters/CharmEW",
        Preambulo=[
            "from LoKiPhys.decorators import *",
            "from LoKiCore.functions  import *",
            "pion_cuts  = in_range ( 300 * MeV , PT , 120 * GeV ) & ( CLONEDIST > 5000 ) & ( TRCHI2DOF < 5 ) ",
            "gamma_cuts = in_range ( 300 * MeV , PT ,  10 * GeV )  ",
            "pions      = SOURCE ( '/Event/Phys/StdAllNoPIDsPions/Particles'  ,  pion_cuts ) ",
            "gammas     = SOURCE ( '/Event/Phys/StdLooseAllPhotons/Particles' , gamma_cuts ) ",
        ],
        Variables={
            "px_c": " pions  >> sum ( PX ) ",
            "py_c": " pions  >> sum ( PY ) ",
            "px_g": " gammas >> sum ( PX ) ",
            "py_g": " gammas >> sum ( PY ) ",
            "n_c": " pions  >> SIZE       ",
            "g_c": " gammas >> SIZE       ",
        }
    )
    from Configurables import DataOnDemandSvc
    dod = DataOnDemandSvc()
    dod.AlgMap['/Event/Counters/CharmEW'] = cnt

    # ========================================================================
    # prefilters for drastical speedup in the reading of input data
    # ========================================================================
    from PhysConf.Filters import LoKi_Filters
    fltrs = LoKi_Filters(
        STRIP_Code=" HLT_PASS_RE ( 'Stripping.*WMuLine.*Decision' ) "
    )

    davinci = DaVinci(
        EventPreFilters=fltrs.filters('Filters'),  # PREFILTERS
        DataType=the_year,
        InputType='DST',
        Simulation=True,
        PrintFreq=10000,
        EvtMax=-1,
        #
        HistogramFile='MCW_Histos.root',
        TupleFile='MCW.root',
        #
    )

    # connect to DaVinci
    from Configurables import GaudiSequencer
    davinci.UserAlgorithms = [
        GaudiSequencer('MySeq', Members=[seq.sequence(), 'MCW'])
    ]

    #
    # take care abotu DB-tags:
    #
    # try to get the tags from Rec/Header
    from BenderTools.GetDBtags import getDBTags
    tags = getDBTags(
        datafiles[0],
        castor
    )
    logger.info('Extract tags from DATA : %s' % tags)
    if tags.has_key('DDDB') and tags['DDDB']:
        davinci.DDDBtag = tags['DDDB']
        logger.info('Set DDDB    %s ' % davinci.DDDBtag)
    if tags.has_key('CONDDB') and tags['CONDDB']:
        davinci.CondDBtag = tags['CONDDB']
        logger.info('Set CONDDB  %s ' % davinci.CondDBtag)
    if tags.has_key('SIMCOND') and tags['SIMCOND']:
        davinci.CondDBtag = tags['SIMCOND']
        logger.info('Set SIMCOND %s ' % davinci.CondDBtag)

    #
    # remove excessive printout
    #
    from Configurables import MessageSvc
    msg = MessageSvc()
    msg.setError += ['HcalDet.Quality',
                     'EcalDet.Quality',
                     'MagneticFieldSvc',
                     'PropertyConfigSvc',
                     'ToolSvc.L0DUConfig',
                     'ToolSvc.L0CondDBProvider',
                     'L0MuonFromRaw',
                     'IntegrateBeamCrossing']

    #
    # come back to Bender
    #
    setData(datafiles, catalogs, castor)

    #
    # start Gaudi
    #
    gaudi = appMgr()

    #
    # more silence
    #
    _a = gaudi.tool('ToolSvc.L0DUConfig')
    _a.OutputLevel = 4

    alg = MCW(
        'MCW',
        Inputs=[seq.outputLocation()],
        PP2MCs=['Relations/Rec/ProtoP/Charged']
    )

    return SUCCESS

# =============================================================================
# The actual job steering
if '__main__' == __name__:

    input = [
        '/lhcb/MC/MC11a/ALLSTREAMS.DST/00015831/0000/00015831_00000163_1.allstreams.dst'
    ]

    configure(input, castor=True)

    run(500)

    gaudi = appMgr()

# =============================================================================
# The END
# =============================================================================
