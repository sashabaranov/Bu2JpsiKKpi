#!/usr/bin/env ipython
# ========================================================================
# @file DiCharm/Bc.py
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
from Bender.Main import *  # import all bender goodies
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
import BenderTools.TisTos  # add methods for TisTos
import DiCharm.Pids  # add methods for Pid/Track information
# =============================================================================
from DiCharm.Alg import DiCharmAlg
# =============================================================================


# =============================================================================
# @class B2Ds
#  Simple template algorithm to study Charm & Y
#  @date   2011-05-27
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
class B2Ds (DiCharmAlg):

    """
    Simple algorithm 
    """

    def initialize(self):

        sc = DiCharmAlg.initialize(self)

        mDs = LoKi.Particles.massFromName('D_s+')

        MASSES = LoKi.AuxDTFBase.MASSES

        self._masses1 = MASSES()
        self._masses1['D_s+'] = mDs - 0.16 * MeV

        self._masses2 = MASSES()
        self._masses2['D_s+'] = mDs + 0.16 * MeV

        self._masses3 = MASSES()
        self._masses3['D_s+'] = mDs - 0.33 * MeV

        self._masses4 = MASSES()
        self._masses4['D_s+'] = mDs + 0.33 * MeV

        return sc

    # the only one esential method:
    def analyse(self):
        """
        The only one essential method
        """
        primaries = self.vselect('PVs', ISPRIMARY)
        if primaries.empty():
            self.Warning('No primary vertices are found', SUCCESS)

        #
        # rec summary
        rc_summary = self.get('/Event/Rec/Summary').summaryData()
        #

        bc = self.select(
            'bc', '[ Beauty -> ( J/psi(1S) -> mu+ mu- ) ( D_s+ --> K+ K- pi+ ) ]CC')

        c2dtf = DTF_CHI2NDOF(True, strings(['J/psi(1S)', 'D_s+']))
        mfit = DTF_FUN(M, True, strings(['J/psi(1S)', 'D_s+']))
        ctau = DTF_CTAU(0, True, strings(['J/psi(1S)', 'D_s+']))
        ctauC = DTF_CTAU(
            'D_s+' == ABSID, True, strings(['J/psi(1S)', 'D_s+']))
        ctauCs = DTF_CTAUSIGNIFICANCE(
            'D_s+' == ABSID, True, strings(['J/psi(1S)', 'D_s+']))

        miniph = MINTREE(('pi+' == ABSID) | ('K+' == ABSID), BPVIPCHI2())

        doca = DOCA(1, 2, self.distanceCalculator())

        mfit1 = DTF_FUN(M, True, 'J/psi(1S)', self._masses1)
        mfit2 = DTF_FUN(M, True, 'J/psi(1S)', self._masses2)
        mfit3 = DTF_FUN(M, True, 'J/psi(1S)', self._masses3)
        mfit4 = DTF_FUN(M, True, 'J/psi(1S)', self._masses4)

        # NDF : 10
        ## ndf  = DTF_NDOF (     True , strings ( [ 'J/psi(1S)' , 'D_s+' ] ) )

        tup = self.nTuple('BDs')
        for b in bc:
            #

            # print 'NDOF: ', ndf ( b )

            m = M(b) / GeV
            if not m < 7.1:
                continue
            #
            psi = b(1)
            ds = b(2)

            self.treatKine(tup, b, '_bc')
            self.treatKine(tup, psi, '_psi')
            self.treatKine(tup, ds, '_ds')

            tup.column('mass', mfit(b) / GeV)
            tup.column('mass1', mfit1(b) / GeV)
            tup.column('mass2', mfit2(b) / GeV)
            tup.column('mass3', mfit3(b) / GeV)
            tup.column('mass4', mfit4(b) / GeV)
            tup.column('c2dtf', c2dtf(b))
            tup.column('ctau', ctau(b))
            tup.column('ctauC', ctauC(b))
            tup.column('ctauCs', ctauCs(b))
            tup.column('minipH', miniph(b))
            tup.column('mKK', M12(ds) / GeV)

            #
            #
            #
            basic = LHCb.Particle.ConstVector()
            hadrons = LHCb.Particle.ConstVector()

            b.children(LoKi.Child.Selector(
                ISBASIC & HASTRACK), basic)
            b.children(LoKi.Child.Selector(
                ('pi+' == ABSID) | ('K+' == ABSID)), hadrons)

            tup.column('doca', doca.docachi2max(basic))
            tup.column('docaH', doca.docachi2max(hadrons))
            #
            #
            #

            # add the information needed for TisTos
            self.tisTos(psi, tup, 'psi_',
                        self.lines['psi'], self.l0tistos, self.tistos)

            # add the information for Pid efficiency correction
            self.treatPions(tup, b)
            self.treatKaons(tup, b)
            self.treatMuons(tup, b)
            self.treatTracks(tup, b)

            # add some reco-summary information
            self.addRecSummary(tup, rc_summary)
            #
            tup.write()

        return SUCCESS

# =============================================================================
# @class Bc23pi
#  Simple template algorithm to study Bc -> J/psi 3pi
#  @date   2011-05-27
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class Bc23pi (DiCharmAlg):

    """
    Simple algorithm 
    """
    # the only one esential method:

    def analyse(self):
        """
        The only one essential method
        """

        primaries = self.vselect('PVs', ISPRIMARY)

        if primaries.empty():

            self.Warning('No primary vertices are found', SUCCESS)
            mips = -1000 * PONE()

        else:

            mips = MIPCHI2(self.geo(), primaries)

        #
        # rec summary
        rc_summary = self.get('/Event/Rec/Summary').summaryData()
        #

        bc = self.select(
            'bc', '[ B_c+ -> ( J/psi(1S) -> mu+ mu- ) pi+ pi+ pi-]CC')

        c2dtf = DTF_CHI2NDOF(True, 'J/psi(1S)')
        mfit = DTF_FUN(M, True, 'J/psi(1S)')
        ctau = DTF_CTAU(0, True, strings(['J/psi(1S)', 'D_s+']))

        miniph = MINTREE(('pi+' == ABSID) | ('K+' == ABSID), BPVIPCHI2())

        doca = DOCA(1, 2, self.distanceCalculator())

        tup = self.nTuple('B3pi')
        for b in bc:
            #
            m = M(b) / GeV
            if not 5.9 < m < 6.6:
                continue
            #
            psi = b(1)
            pi1 = b(2)
            pi2 = b(3)
            pi3 = b(4)

            self.treatKine(tup, b, '_bc')
            self.treatKine(tup, psi, '_psi')
            self.treatKine(tup, pi1, '_pi1')
            self.treatKine(tup, pi2, '_pi2')
            self.treatKine(tup, pi3, '_pi3')

            tup.column('mass', mfit(b) / GeV)
            tup.column('c2dtf', c2dtf(b))
            tup.column('ctau', ctau(b))
            tup.column('minipH', miniph(b))

            #
            # fictive Ds
            p = LoKi.LorentzVector()
            p += pi1.momentum()
            p += pi2.momentum()
            p += pi3.momentum()

            tup.column('m_ds', p.M() / GeV)
            tup.column('p_ds', p.P() / GeV)
            tup.column('pt_ds', p.Pt() / GeV)

            #
            #
            #
            basic = LHCb.Particle.ConstVector()
            hadrons = LHCb.Particle.ConstVector()

            b.children(LoKi.Child.Selector(
                ISBASIC & HASTRACK), basic)
            b.children(LoKi.Child.Selector(
                ('pi+' == ABSID) | ('K+' == ABSID)), hadrons)

            tup.column('doca', doca.docachi2max(basic))
            tup.column('docaH', doca.docachi2max(hadrons))
            #
            #
            #

            # add the information needed for TisTos
            self.tisTos(psi, tup, 'psi_',
                        self.lines['psi'],
                        self.l0tistos,
                        self.tistos)

            # add the information for Pid efficiency correction
            self.treatPions(tup, b)
            self.treatMuons(tup, b)
            # add the infomation needed for track efficiency correction
            self.treatTracks(tup, b)
            self.treatKaons(tup, b)

            # add some reco-summary information
            self.addRecSummary(tup, rc_summary)
            #
            tup.write()

        return SUCCESS

# =============================================================================
# @class Bc2pi
#  Simple template algorithm to study Bc- >J/psi pi
#  @date   2011-05-27
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class Bc2pi (DiCharmAlg):

    """
    Simple algorithm 
    """
    # the only one esential method:

    def analyse(self):
        """
        The only one essential method
        """

        primaries = self.vselect('PVs', ISPRIMARY)
        if primaries.empty():
            self.Warning('No primary vertices are found', SUCCESS)

        #
        # rec summary
        rc_summary = self.get('/Event/Rec/Summary').summaryData()
        #

        bc = self.select('bc', '[ B_c+ -> ( J/psi(1S) -> mu+ mu- ) pi+ ]CC')

        c2dtf = DTF_CHI2NDOF(True, 'J/psi(1S)')
        mfit = DTF_FUN(M, True, 'J/psi(1S)')
        ctau = DTF_CTAU(0, True, strings(['J/psi(1S)', 'D_s+']))
        miniph = MINTREE(('pi+' == ABSID) | ('K+' == ABSID), BPVIPCHI2())

        doca = DOCA(1, 2, self.distanceCalculator())

        tup = self.nTuple('B1pi')
        for b in bc:
            #
            m = M(b) / GeV
            if not 5.9 < m < 6.6:
                continue
            #
            psi = b(1)
            pi1 = b(2)

            self.treatKine(tup, b, '_bc')
            self.treatKine(tup, psi, '_psi')
            self.treatKine(tup, pi1, '_pi1')

            tup.column('mass', mfit(b) / GeV)
            tup.column('c2dtf', c2dtf(b))
            tup.column('ctau', ctau(b))
            tup.column('minipH', miniph(b))

            #
            # fictive Ds
            #
            tup.column('m_ds', M(pi1) / GeV)
            tup.column('p_ds', P(pi1) / GeV)
            tup.column('pt_ds', PT(pi1) / GeV)

            #
            #
            #
            basic = LHCb.Particle.ConstVector()

            b.children(LoKi.Child.Selector(ISBASIC & HASTRACK), basic)

            tup.column('doca', doca.docachi2max(basic))
            #
            #
            #

            # add the information needed for TisTos
            self.tisTos(psi, tup, 'psi_',
                        self.lines['psi'],
                        self.l0tistos,
                        self.tistos)

            # add the information for Pid efficiency correction
            self.treatPions(tup, b)
            self.treatMuons(tup, b)
            # add the infomation needed for track efficiency correction
            self.treatTracks(tup, b)
            self.treatKaons(tup, b)

            # add some reco-summary information
            self.addRecSummary(tup, rc_summary)
            #
            tup.write()

        return SUCCESS

# =============================================================================
# @class B2D0
#  Simple template algorithm to study
#  @date   2011-05-27
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class B2D0 (DiCharmAlg):

    """
    Simple algorithm 
    """
    # the only one esential method:

    def analyse(self):
        """
        The only one essential method
        """
        primaries = self.vselect('PVs', ISPRIMARY)
        if primaries.empty():
            self.Warning('No primary vertices are found', SUCCESS)

        #
        # rec summary
        #
        rc_summary = self.get('/Event/Rec/Summary').summaryData()
        #

        bc = self.select(
            'bc', '[ Beauty -> ( J/psi(1S) -> mu+ mu- ) ( D0 -> K- pi+ ) ]CC')

        c2dtf = DTF_CHI2NDOF(True, strings(['J/psi(1S)', 'D0']))
        mfit = DTF_FUN(
            M, True, strings(['J/psi(1S)', 'D0']))
        ctau = DTF_CTAU(0, True, strings(['J/psi(1S)', 'D0']))
        ctauC = DTF_CTAU(
            'D0' == ABSID, True, strings(['J/psi(1S)', 'D0']))
        ctauCs = DTF_CTAUSIGNIFICANCE(
            'D0' == ABSID, True, strings(['J/psi(1S)', 'D0']))

        ct1_ = COSPOL('mu+' == ID, 'J/psi(1S)' == ID)
        ct2_ = COSPOL('K+' == ABSID, 'D0' == ABSID)

        ct1 = DTF_FUN(ct1_, True, strings(['J/psi(1S)', 'D0']))
        ct2 = DTF_FUN(ct2_, True, strings(['J/psi(1S)', 'D0']))

        tup = self.nTuple('BD0')
        for b in bc:
            #
            m = M(b) / GeV
            if not m < 7.1:
                continue
            #
            psi = b(1)
            d0 = b(2)

            self.treatKine(tup, b, '_bc')
            self.treatKine(tup, psi, '_psi')
            self.treatKine(tup, d0, '_d0')

            tup.column('mass', mfit(b) / GeV)
            tup.column('c2dtf', c2dtf(b))
            tup.column('ctau', ctau(b))
            tup.column('ctauC', ctauC(b))
            tup.column('ctauCs', ctauCs(b))

            tup.column('ctpsi', ct1(b))
            tup.column('ctd0', ct2(b))

            # add the information needed for TisTos
            self.tisTos(psi, tup, 'psi_',
                        self.lines['psi'], self.l0tistos, self.tistos)

            # add the information for Pid efficiency correction
            self.treatPions(tup, b)
            self.treatKaons(tup, b)
            self.treatMuons(tup, b)
            self.treatTracks(tup, b)

            # add some reco-summary information
            self.addRecSummary(tup, rc_summary)
            #
            tup.write()

        return SUCCESS

# =============================================================================
# @class B2Dp
#  Simple template algorithm to study
#  @date   2011-05-27
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class B2Dp (DiCharmAlg):

    """
    Simple algorithm 
    """
    # the only one esential method:

    def analyse(self):
        """
        The only one essential method
        """
        primaries = self.vselect('PVs', ISPRIMARY)
        if primaries.empty():
            self.Warning('No primary vertices are found', SUCCESS)

        #
        # rec summary
        rc_summary = self.get('/Event/Rec/Summary').summaryData()
        #

        bc = self.select(
            'bc', '[ Beauty -> ( J/psi(1S) -> mu+ mu- ) ( D+ -> K- pi+ pi+ ) ]CC')

        c2dtf = DTF_CHI2NDOF(True, strings(['J/psi(1S)', 'D+']))
        mfit = DTF_FUN(M, True, strings(['J/psi(1S)', 'D+']))
        ctau = DTF_CTAU(0, True, strings(['J/psi(1S)', 'D+']))
        ctauC = DTF_CTAU(
            'D+' == ABSID, True, strings(['J/psi(1S)', 'D+']))
        ctauCs = DTF_CTAUSIGNIFICANCE(
            'D+' == ABSID, True, strings(['J/psi(1S)', 'D+']))

        tup = self.nTuple('BDp')
        for b in bc:
            #
            m = M(b) / GeV
            if not m < 7.1:
                continue
            #
            psi = b(1)
            dp = b(2)

            self.treatKine(tup, b, '_bc')
            self.treatKine(tup, psi, '_psi')
            self.treatKine(tup, dp, '_dp')

            tup.column('mass', mfit(b) / GeV)
            tup.column('c2dtf', c2dtf(b))
            tup.column('ctau', ctau(b))
            tup.column('ctauC', ctauC(b))
            tup.column('ctauCs', ctauCs(b))

            # add the information needed for TisTos
            self.tisTos(psi, tup, 'psi_',
                        self.lines['psi'], self.l0tistos, self.tistos)

            # add the information for Pid efficiency correction
            self.treatPions(tup, b)
            self.treatKaons(tup, b)
            self.treatMuons(tup, b)
            self.treatTracks(tup, b)

            # add some reco-summary information
            self.addRecSummary(tup, rc_summary)
            #
            tup.write()

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

    rootInTES = '/Event/B2PSIC'

    davinci = DaVinci(
        DataType=the_year,
        InputType='MDST',
        RootInTES=rootInTES,
        Simulation=False,
        PrintFreq=1000,
        EvtMax=-1,
        #
        HistogramFile='Bc_Histos.root',
        TupleFile='Bc.root',
        #
        Lumi=True,
        #
    )

    from Configurables import MessageSvc
    msg = MessageSvc()
    msg.setError += ['HcalDet.Quality',
                     'EcalDet.Quality',
                     'MagneticFieldSvc',
                     'PropertyConfigSvc']

    from Configurables import TrackScaleState
    state_scale = TrackScaleState(
        'StateScale',
        ## RootInTES  = rootInTES  ,
    )

    ##from Configurables import TrackSmearState
    # state_smear    = TrackSmearState (
    ##    'StateSmear'            ,
    # RootInTES  = rootInTES  ,
    # )

    davinci.UserAlgorithms = [state_scale,
                              ## state_smear ,
                              'BDs',
                              #'B3pi' ,
                              #'B1pi' ,
                              #'BD0'  ,
                              #'BDp'
                              ]

    from Configurables import CondDB
    #CondDB ( LatestGlobalTagByDataType = the_year )
    ## CondDB ( LatestLocalTagsByDataType = [ the_year ]  )

    from Configurables import CondDBAccessSvc
    CondDB(). addLayer(CondDBAccessSvc(
        "myCond",
        ## ConnectionString = "sqlite_file:SCALE%s.db/LHCBCOND" % the_year,
        ConnectionString="sqlite_file:$HOME/tmp/SCALE.db/LHCBCOND",
        DefaultTAG="HEAD"))

    # ------- decoding set-up start ----------
    #from BenderTools.MicroDST import uDstConf
    #uDstConf ( rootInTES )
    # ------- decoding set-up end  -----------

    #
    # come back to Bender
    #
    setData(datafiles, catalogs, castor)

    #
    # start Gaudi
    #
    gaudi = appMgr()

    #
    alg1 = B2Ds(
        'BDs',  # Algorithm name ,
        RootInTES=rootInTES,
        Inputs=['Phys/PsiDs/Particles']
    )

    alg2 = Bc23pi(
        'B3pi',  # Algorithm name ,
        RootInTES=rootInTES,
        Inputs=['Phys/Psi3pi/Particles']
    )
    #
    alg3 = Bc2pi(
        'B1pi',  # Algorithm name ,
        RootInTES=rootInTES,
        Inputs=['Phys/Psi1Pi/Particles']
    )
    #
    alg4 = B2D0(
        'BD0',  # Algorithm name ,
        RootInTES=rootInTES,
        Inputs=['Phys/PsiD0/Particles']
    )

    alg5 = B2Dp(
        'BDp',  # Algorithm name ,
        RootInTES=rootInTES,
        Inputs=['Phys/PsiDp/Particles']
    )

##     mainSeq = gaudi.algorithm ('GaudiSequencer/DaVinciUserSequence', True )
# mainSeq.Members += [ alg1 . name () ,
##                          alg2 . name () ,
##                          alg3 . name () ,
##                          alg4 . name () ,
# alg5 . name () ]

    return SUCCESS

# =============================================================================
# The actual job steering
if '__main__' == __name__:

    import shelve
    db = shelve.open('/afs/cern.ch/user/i/ibelyaev/public/LFNs/lfns.db')

    bc = []
    # bc += db [ 'J/psiC;MagnetDown;20'  ]
    # bc += db [ 'J/psiC;MagnetUp;20'    ]
    bc += db['J/psiC;MagnetDown;17b']
    bc += db['J/psiC;MagnetUp;17b']

    input = bc

    db.close()

    import DetCond.HistoCond

    configure(input, castor=True)

    run(100)

    gaudi = appMgr()

    det = gaudi.detSvc()
    cond = det['/dd/Conditions/Calibration/LHCb/MomentumScale']
    cpp.ISolid
    DetDesc = cpp.DetDesc

    #h1      = DetDesc.Params.paramAsHisto2D ( cond , 'IdpPlus'  )
    #h2      = DetDesc.Params.paramAsHisto2D ( cond , 'IdpMinus' )
    #offsets = DetDesc.Params.paramAsHisto1D ( cond , 'Offsets'  )

    h1 = cond .paramAsHisto2D('IdpPlus')
    h2 = cond .paramAsHisto2D('IdpMinus')
    offsets = cond .paramAsHisto1D('Offsets')
    delta = cond .paramAsDouble('Delta')

# =============================================================================
# The END
# =============================================================================
