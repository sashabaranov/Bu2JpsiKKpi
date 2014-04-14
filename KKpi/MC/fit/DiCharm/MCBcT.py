#!/usr/bin/env ipython
# ========================================================================
# @file DiCharm/MCBc.py
#
#  Simple template algorithm for Bc
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

Simple template algorithm for Bc

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
from GaudiKernel.SystemOfUnits import GeV, MeV, mm
from GaudiKernel.PhysicalConstants import c_light
# =============================================================================
# logging
# =============================================================================
from Bender.Logger import getLogger
logger = getLogger(__name__)
# =============================================================================
import BenderTools.TisTosMC  # add methods for TisTos
import BenderTools.Fill  # add methods for Fill
# =============================================================================

# =============================================================================
# @class AlgB
#  Simple template algorithm to study Bc
#  @date   2011-05-27
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class AlgB (Algo):

    """
    Simple algorithm 
    """
    # standard constructor

    def __init__(self, name, **kwargs):
        """
        Standard constructor
        """
        Algo.__init__(self, name, **kwargs)

    def initialize(self):
        """
        Initialization
        """
        sc = Algo.    initialize(self)
        if sc.isFailure():
            return sc

        triggers = {}
        triggers['J/psi(1S)'] = {}

        lines = {}
        lines['psi'] = {}
        lines['psi_det'] = {}
        lines['psi_unb'] = {}

        lines["psi"]['L0TOS'] = 'L0(DiMuon|Muon)Decision'
        lines["psi"][
            'L0TIS'] = 'L0(Hadron|DiMuon|Muon|Electron|Photon)Decision'
        lines["psi"][
            'Hlt1TOS'] = 'Hlt1(DiMuon|SingleMuon|TrackMuon).*Decision'
        lines["psi"][
            'Hlt1TIS'] = 'Hlt1(DiMuon|SingleMuon|Track).*Decision'
        lines["psi"][
            'Hlt2TOS'] = 'Hlt2(DiMuon|ExpressJPsi|SingleMuon).*Decision'
        lines["psi"][
            'Hlt2TIS'] = 'Hlt2(Charm|Topo|DiMuon|Single).*Decision'

        lines["psi_det"]['L0TOS'] = 'L0(DiMuon|Muon)Decision'
        lines["psi_det"][
            'L0TIS'] = 'L0(Hadron|DiMuon|Muon|Electron|Photon)Decision'
        lines["psi_det"]['Hlt1TOS'] = 'Hlt1DiMuonHighMass.*Decision'
        lines["psi_det"][
            'Hlt1TIS'] = 'Hlt1(DiMuon|SingleMuon|Track).*Decision'
        lines["psi_det"][
            'Hlt2TOS'] = 'Hlt2DiMuonDetached(Heavy|JPsi)Decision'
        lines["psi_det"][
            'Hlt2TIS'] = 'Hlt2(Charm|Topo|DiMuon|Single).*Decision'

        lines["psi_unb"]['L0TOS'] = 'L0(DiMuon|Muon)Decision'
        lines["psi_unb"][
            'L0TIS'] = 'L0(Hadron|DiMuon|Muon|Electron|Photon)Decision'
        lines["psi_unb"]['Hlt1TOS'] = 'Hlt1DiMuonHighMass.*Decision'
        lines["psi_unb"][
            'Hlt1TIS'] = 'Hlt1(DiMuon|SingleMuon|Track).*Decision'
        lines["psi_unb"]['Hlt2TOS'] = 'Hlt2DiMuonJPsiHighPT.*Decision'
        lines["psi_unb"][
            'Hlt2TIS'] = 'Hlt2(Charm|Topo|DiMuon|Single).*Decision'

        sc = self.tisTos_initialize(triggers, lines)
        if sc.isFailure():
            return sc

        sc = self.  fill_initialize()
        if sc.isFailure():
            return sc

        return SUCCESS

    # finalize it!
    def finalize(self):
        """
        finalize the algorithm
        """
        #
        self . tisTos_finalize()
        self .   fill_finalize()
        #
        self . dumpHistos()
        #
        return Algo.finalize(self)

    # the only one esential method:
    def analyse(self):
        """
        The only one essential method
        """
        primaries = self.vselect('PVs', ISPRIMARY)
        if primaries.empty():
            self.Warning('No primary vertices are found', SUCCESS)

        #
        # rec summary &ODIN
        #
        rc_summary = self.get('/Event/Rec/Summary').summaryData()
        odin = self.get_('/Event/DAQ/ODIN', True)
        #

        bc = self.  select('bc',  self.decay)  # NB: self.decay!

        if bc.empty():
            return SUCCESS

        c2dtf_fit = DTF_CHI2NDOF(True, strings(['J/psi(1S)']))
        mfit_fit = DTF_FUN(M, True, strings(['J/psi(1S)']))
        ctau_fit = DTF_CTAU(0, True, strings(['J/psi(1S)']))

        c2dtf = DTF_CHI2NDOF(True)
        mfit = DTF_FUN(M, True)
        mpsi = DTF_FUN(M1, True)
        ctau = DTF_CTAU(0, True)

        tup = self.nTuple('B')
        for b in bc:
            #
            m = M(b) / GeV
            if not m < 7.1:
                continue
            #
            self.plot(m, 'B-mass', 5, 7)
            #

            psi = b(1)
            had = b(2)

            #
            # with J/psi constraint:
            #
            tup.column_float('mass', mfit_fit(b) / GeV)
            tup.column_float('c2dtf', c2dtf_fit(b))
            tup.column_float('ctau', ctau_fit(b))

            #
            # without J/psi mass constraint
            #
            tup.column_float('mass_0', mfit(b) / GeV)
            tup.column_float('m1_0', mpsi(b) / GeV)
            tup.column_float('ctau_0', ctau(b))
            tup.column_float('c2dtf_0', c2dtf(b))

            #
            # general kinematic
            #
            self.treatKine(tup, b, "_b")
            self.treatKine(tup, psi, "_psi")
            self.treatKine(tup, had, "_h")

            # add the information needed for TisTos
            self.tisTos(psi, tup, 'psi0_',
                        self.lines['psi'],
                        self.l0tistos,
                        self.tistos)

            self.tisTos(psi, tup, 'psi1_',
                        self.lines['psi_det'],
                        self.l0tistos,
                        self.tistos)

            self.tisTos(psi, tup, 'psi2_',
                        self.lines['psi_unb'],
                        self.l0tistos,
                        self.tistos)

            # add the information for Pid efficiency correction
            self.treatPions(tup, b)
            self.treatKaons(tup, b)
            self.treatMuons(tup, b)
            self.treatTracks(tup, b)

            # add some reco-summary information
            self.addRecSummary(tup, rc_summary)
            # add ODIN information
            tup.column_aux(odin)

            tup.write()

        return SUCCESS

# =============================================================================
# @class AlgBc
#  simple class to fill n-tupel for Bc-decays


class AlgBc (AlgB):

    """
    Simple class to fill n-tuple for Bc-decays 
    """

    def __init__(self, name, **kwargs):

        AlgB.__init__(self, name, **kwargs)

        self.decay = "[ B_c+ -> ( J/psi(1S) -> mu+ mu- ) pi+]CC "

# =============================================================================
# @class AlgBu
#  simple class to fill n-tuple for Bu-decays


class AlgBu (AlgB):

    """
    Simple class to fill n-tuple for B+-decays 
    """

    def __init__(self, name, **kwargs):

        AlgB.__init__(self, name, **kwargs)

        self.decay = "[ B+ -> ( J/psi(1S) -> mu+ mu- ) K+ ]CC "

# =============================================================================
# configure the job


def configure_COMMON(datafiles,
                     catalogs,
                     castor,
                     params):
    """
    The actual job configuration 
    """
    # =====================================================================
    from Configurables import DaVinci  # needed for job configuration

    the_year = params['Year']

    rootInTES = '/Event/B2PSI'

    davinci = DaVinci(
        #
        DataType=params['Year'],
        RootInTES=rootInTES,
        InputType='MDST',
        #
        PrintFreq=1000,
        EvtMax=-1,
        #
    )

    #
    # define input files
    #
    setData(datafiles, catalogs, castor)

    #
    # suppress some extra printput
    #
    from BenderTools.Utils import silence
    silence()

    mode = params.get('Mode', '')
    logger.info("Mode is '%s' " % mode)

    mode = mode.upper()

    if not mode:

        davinci.UserAlgorithms += ['Bc_det', 'Bc_unb']
        davinci.UserAlgorithms += ['Bu_det', 'Bu_unb']

    elif 0 <= mode.find('BC') or 0 <= mode.find('PI'):

        davinci.UserAlgorithms += ['Bc_det', 'Bc_unb']

    elif 0 <= mode.find('B+') or 0 <= mode.find('BU') or 0 <= mode.find('K+'):

        davinci.UserAlgorithms += ['Bu_det', 'Bu_unb']

    else:

        raise AttibuteError, 'Mode is not specified!'

    gaudi = appMgr()

    bc_det = AlgBc('Bc_det',
                   RootInTES='/Event/B2PSI',
                   Inputs=['Phys/TheBc_det/Particles'],
                   ReFitPVs=True)

    bc_unb = AlgBc('Bc_unb',
                   RootInTES='/Event/B2PSI',
                   Inputs=['Phys/TheBc_unb/Particles'],
                   ReFitPVs=True)

    bu_det = AlgBu('Bu_det',
                   RootInTES='/Event/B2PSI',
                   Inputs=['Phys/TheB_det/Particles'],
                   ReFitPVs=True)

    bu_unb = AlgBu('Bu_unb',
                   RootInTES='/Event/B2PSI',
                   Inputs=['Phys/TheB_unb/Particles'],
                   ReFitPVs=True)

    return SUCCESS

# =============================================================================
# configure the job


def configure_DATA(datafiles,
                   catalogs,
                   castor,
                   params):

    logger.info('Configure DATA')

    the_year = params['Year']

    from Configurables import CondDB
    CondDB(LatestGlobalTagByDataType=the_year)

    from Configurables import DaVinci
    dv = DaVinci(Lumi=True,
                 #
                 HistogramFile='B_Histos.root',
                 TupleFile='B.root'
                 )

    from Configurables import TrackScaleState
    state_scale = TrackScaleState('StateScale')

    dv.UserAlgorithms = [state_scale]
    logger.info('Momentum scaling is applied!')

    return configure_COMMON(datafiles, catalogs, castor, params)


# =============================================================================
# configure the job
def configure_MC(datafiles,
                 catalogs,
                 castor,
                 params):

    logger.info('Configure MC')
    the_year = params['Year']

    from Configurables import DaVinci
    dv = DaVinci(Lumi=False,
                 Simulation=True,
                 DDDBtag=params['DDDB'],
                 CondDBtag=params['SIMCOND'],
                 #
                 HistogramFile='MCB_Histos.root',
                 TupleFile='MCB.root')

    logger.info('DDDB   tag is %s ' % params['DDDB'])
    logger.info('CondDB tag is %s ' % params['SIMCOND'])

    from PhysConf.MicroDST import uDstConf
    uDstConf('/Event/B2PSI',
             logger=logger)  # use it!

    return configure_COMMON(datafiles, catalogs, castor, params)

# =============================================================================


def configure(datafiles,
              catalogs=[],
              castor=False,
              params={}):

    data = params.has_key('DATA') and params['DATA']
    data = data or (params.has_key('data') and params['data'])
    data = data or (params.has_key('MC') and not params['MC'])
    data = data or (params.has_key('mc') and not params['mc'])

    mc = params.has_key('DATA') and not params['DATA']
    mc = mc or (params.has_key('data') and not params['data'])
    mc = mc or (params.has_key('MC') and params['MC'])
    mc = mc or (params.has_key('mc') and params['mc'])

    if data:
        return configure_DATA(datafiles, catalogs, castor, params)
    elif mc:
        return configure_MC(datafiles, catalogs, castor, params)
    else:
        raise AttibuteError, 'DATA/MC is not specified!'

# =============================================================================
# The actual job steering
if '__main__' == __name__:

    import shelve
    db = shelve.open('/afs/cern.ch/user/i/ibelyaev/public/LFNs/lfns.db', 'r')
    keys = db.keys()
    keys.sort()
    for k in keys:
        if 0 != k.find('Bc-ilfetime'):
            continue
        print 'Valid key:', k
    #
    # DATA:
    #
    inputs, params = db['Bc-lifetime:2k+11,MagDown']
    inputs, params = db['Bc-lifetime:2k+11,MagUp']
    inputs, params = db['Bc-lifetime:2k+12,MagDown']
    inputs, params = db['Bc-lifetime:2k+12,MagUp']

#
# MC :
# B+ -> J/psiK+
# Bu_keys = [
##         'Bc-lifetime:B+->J/psiK+;MC/2011;2011;MagDown;Pythia6;Reco14a' ,
##         'Bc-lifetime:B+->J/psiK+;MC/2011;2011;MagDown;Pythia8;Reco14a' ,
##         'Bc-lifetime:B+->J/psiK+;MC/2011;2011;MagUp;Pythia6;Reco14a'   ,
##         'Bc-lifetime:B+->J/psiK+;MC/2011;2011;MagUp;Pythia8;Reco14a'   ,
##         'Bc-lifetime:B+->J/psiK+;MC/2012;2012;MagDown;Pythia6;Reco14a' ,
##         'Bc-lifetime:B+->J/psiK+;MC/2012;2012;MagDown;Pythia8;Reco14a' ,
##         'Bc-lifetime:B+->J/psiK+;MC/2012;2012;MagUp;Pythia6;Reco14a'   ,
##         'Bc-lifetime:B+->J/psiK+;MC/2012;2012;MagUp;Pythia8;Reco14a'   ,
##         'Bc-lifetime:B+->J/psiK+;MC/2012;2012JulSep2012;MagDown;Pythia6;Reco14'     ,
##         'Bc-lifetime:B+->J/psiK+;MC/2012;2012JulSep2012;MagUp;Pythia6;Reco14'       ,
##         'Bc-lifetime:B+->J/psiK+;MC/2012;2012MayJune2012;MagDown;Pythia6;Reco13a'   ,
##         'Bc-lifetime:B+->J/psiK+;MC/2012;2012MayJune2012;MagDown;Pythia6;Reco13a_1' ,
##         'Bc-lifetime:B+->J/psiK+;MC/2012;2012MayJune2012;MagDown;Pythia6;Reco14' ,
##         'Bc-lifetime:B+->J/psiK+;MC/2012;2012MayJune2012;MagUp;Pythia6;Reco13a'  ,
##         'Bc-lifetime:B+->J/psiK+;MC/2012;2012MayJune2012;MagUp;Pythia6;Reco14'   ,
##         'Bc-lifetime:B+->J/psiK+;MC/MC11a;2011;MagDown;Pythia6;Reco12'    ,
##         'Bc-lifetime:B+->J/psiK+;MC/MC11a;2011;MagDown;Pythia6;Reco12_1'  ,
##         'Bc-lifetime:B+->J/psiK+;MC/MC11a;2011;MagDown;Pythia6;Reco12a'   ,
##         'Bc-lifetime:B+->J/psiK+;MC/MC11a;2011;MagDown;Pythia6;Reco12a_1' ,
##         'Bc-lifetime:B+->J/psiK+;MC/MC11a;2011;MagUp;Pythia6;Reco12'    ,
##         'Bc-lifetime:B+->J/psiK+;MC/MC11a;2011;MagUp;Pythia6;Reco12_1'  ,
##         'Bc-lifetime:B+->J/psiK+;MC/MC11a;2011;MagUp;Pythia6;Reco12a'   ,
# 'Bc-lifetime:B+->J/psiK+;MC/MC11a;2011;MagUp;Pythia6;Reco12a_1'
# ]
##     inputs,params = db[ Bu_keys[0] ]

# Bc -> J/psi pi
# Bc_keys = [
##         'Bc-lifetime:Bc->J/psiPi+;MC/2012;2012;MagDown;Pythia6;Reco14a'           ,
##         'Bc-lifetime:Bc->J/psiPi+;MC/2012;2012;MagUp;Pythia6;Reco14a'             ,
##         'Bc-lifetime:Bc->J/psiPi+;MC/2012;2012JulSep2012;MagDown;Pythia6;Reco14'  ,
##         'Bc-lifetime:Bc->J/psiPi+;MC/2012;2012JulSep2012;MagUp;Pythia6;Reco14'    ,
##         'Bc-lifetime:Bc->J/psiPi+;MC/2012;2012MayJune2012;MagDown;Pythia6;Reco14' ,
##         'Bc-lifetime:Bc->J/psiPi+;MC/2012;2012MayJune2012;MagUp;Pythia6;Reco14'   ,
##         'Bc-lifetime:Bc->J/psiPi+;MC/MC11a;2011;MagDown;Pythia6;Reco12a'          ,
# 'Bc-lifetime:Bc->J/psiPi+;MC/MC11a;2011;MagUp;Pythia6;Reco12a'
# ]
##     inputs,params = db[ Bc_keys[1] ]

    db.close()

    configure(inputs,
              catalogs=[],
              castor=True, params=params)

    run(10000)
# =============================================================================
# The END
# =============================================================================
