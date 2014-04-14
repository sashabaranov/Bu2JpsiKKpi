#!/usr/bin/env ipython
# ========================================================================
# @file WG1.py
#
#  Access to WG/PsiX production
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
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#
#                    $Revision$
#  Last modification $Date$
#                 by $Author$
# =============================================================================
"""

Simple template algorithm for B -> psi 5h

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
__author__ = 'Vanya BELYAEV  Ivan.Belyaev@itep.ru'
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
# @class  Psi5
#  Simple  template algorithm for B -> psi 5h
#  @date   2011-05-27
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru


class Psi5 (DiCharmAlg):

    """
    Simple  template algorithm for B -> psi 5h 
    """

    def initialize(self):
        """
        initialize it! 
        """
        sc = DiCharmAlg.initialize(self)
        if sc.isFailure():
            return sc

        self.c2dtf = DTF_CHI2NDOF(
            True, strings(['J/psi(1S)', 'psi(2S)']))
        self.mfit = DTF_FUN(
            M, True, strings(['J/psi(1S)', 'psi(2S)']))
        self.ctau = DTF_CTAU(
            0, True, strings(['J/psi(1S)', 'psi(2S)']))

        self.nK = NINTREE('K+' == ABSID)
        self.nPi = NINTREE('pi+' == ABSID)

        self.doca = DOCA(1, 2, self.distanceCalculator())

        return SUCCESS

    def finalize(self):

        self.c2dtf = None
        self.mfit = None
        self.ctau = None
        self.nK = None
        self.nPi = None
        self.doca = None

        return DiCharmAlg.finalize(self)

    # fill infomation about masses
    def _fill_mass_(self, b, tup):
        """
        fill infomation about masses         
        """
        nc = b.nChildren()
        for i in range(1, nc + 1):
            for j in range(i + 1, nc + 1):
                m2 = MASS(i, j)
                tup.column_float('m%s%s' % (i, j), m2(b) / GeV)
                for k in range(j + 1, nc + 1):
                    m3 = MASS(i, j, k)
                    tup.column_float('m%s%s%s' % (i, j, k), m3(b) / GeV)
                    for l in range(k + 1, nc + 1):
                        m4 = MASS(i, j, k, l)
                        tup.column_float(
                            'm%s%s%s%s' % (i, j, k, l), m4(b) / GeV)

        tup.column_int('nK', int(self.nK(b)))
        tup.column_int('nPi', int(self.nPi(b)))

        return SUCCESS

    # fill infomation about docamax and chi2_docamax
    def _fill_doca_(self, b, tup):
        """
        Fill infomation about docamax and chi2_docamax 
        """
        basic = LHCb.Particle.ConstVector()
        hadrons = LHCb.Particle.ConstVector()

        b.children(LoKi.Child.Selector(
            ISBASIC & HASTRACK), basic)
        b.children(LoKi.Child.Selector(
            ('pi+' == ABSID) | ('K+' == ABSID)), hadrons)

        tup.column_float('doca', self.doca.docachi2max(basic))
        tup.column_float('docaH', self.doca.docachi2max(hadrons))

        return SUCCESS

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
        # 1
        B = self.select(
            'b', '[Beauty -> ( Meson -> mu+ mu- ) X+ X+ X+ X- X-]CC')

        tup = self.nTuple('B')
        for b in B:
            #
            m = M(b) / GeV
            if not 5.0 < m < 6.7:
                continue
            #

            psi = b(1)

            mpsi = M(psi)
            if 4.5 * GeV < mpsi:
                continue
            elif 3.3 * GeV < mpsi:
                psi.setParticleID(LHCb.ParticleID(100443))

            p1 = b(2)
            p2 = b(3)
            p3 = b(4)
            n1 = b(5)
            n2 = b(6)

            self.treatKine(tup, b, '_b')
            self.treatKine(tup, psi, '_psi')
            self.treatKine(tup, p1, '_p1')
            self.treatKine(tup, p2, '_p2')
            self.treatKine(tup, p3, '_p3')
            self.treatKine(tup, n1, '_n1')
            self.treatKine(tup, n2, '_n2')

            # add the information for Pid efficiency correction
            self.treatPions(tup, b)
            self.treatKaons(tup, b)
            self.treatMuons(tup, b)
            self.treatTracks(tup, b)

            tup.column('m', M(b) / GeV)
            tup.column('mass', self.mfit(b) / GeV)
            tup.column('c2dtf', self.c2dtf(b))
            tup.column('ctau', self.ctau(b))

            self._fill_doca_(b, tup)
            self._fill_mass_(b, tup)

            # add some reco-summary information
            self.addRecSummary(tup, rc_summary)

            # add the information needed for TisTos
            self.tisTos(psi, tup, 'psi_',
                        self.lines['psi'], self.l0tistos, self.tistos)

            tup.write()

        return SUCCESS

# ========================================================================
# @class  Psi4
#  Simple  template algorithm for B -> psi 4h
#  @date   2011-05-27
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru


class Psi4 (Psi5):

    """
    Simple  template algorithm for B -> psi 4h
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
        # 1
        B = self.select(
            'b', '[Beauty -> ( Meson -> mu+ mu- ) X+ X+ X- X- ]CC')

        tup = self.nTuple('B')
        for b in B:
            #
            m = M(b) / GeV
            if not 5.0 < m < 6.7:
                continue
            #

            psi = b(1)

            mpsi = M(psi)
            if 4.5 * GeV < mpsi:
                continue
            elif 3.3 * GeV < mpsi:
                psi.setParticleID(LHCb.ParticleID(100443))

            p1 = b(2)
            p2 = b(3)
            n1 = b(4)
            n2 = b(5)

            self.treatKine(tup, b, '_b')
            self.treatKine(tup, psi, '_psi')
            self.treatKine(tup, p1, '_p1')
            self.treatKine(tup, p2, '_p2')
            self.treatKine(tup, n1, '_n1')
            self.treatKine(tup, n2, '_n2')

            # add the information for Pid efficiency correction
            self.treatPions(tup, b)
            self.treatKaons(tup, b)
            self.treatMuons(tup, b)
            self.treatTracks(tup, b)

            tup.column('m', M(b) / GeV)
            tup.column('mass', self.mfit(b) / GeV)
            tup.column('c2dtf', self.c2dtf(b))
            tup.column('ctau', self.ctau(b))

            self._fill_doca_(b, tup)
            self._fill_mass_(b, tup)

            # add some reco-summary information
            self.addRecSummary(tup, rc_summary)

            # add the information needed for TisTos
            self.tisTos(psi, tup, 'psi_',
                        self.lines['psi'], self.l0tistos, self.tistos)

            tup.write()

        return SUCCESS

# ========================================================================
# @class  Psi6
#  Simple  template algorithm for B -> psi 6
#  @date   2011-05-27
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru


class Psi6 (Psi5):

    """
    Simple  template algorithm for B -> psi 6
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
        # 1
        B = self.select(
            'b', '[Beauty -> ( Meson -> mu+ mu- ) X+ X+ X+ X- X- X-]CC')

        tup = self.nTuple('B')
        for b in B:
            #
            m = M(b) / GeV
            if not 5.0 < m < 6.7:
                continue
            #

            psi = b(1)

            mpsi = M(psi)
            if 4.5 * GeV < mpsi:
                continue
            elif 3.3 * GeV < mpsi:
                psi.setParticleID(LHCb.ParticleID(100443))

            p1 = b(2)
            p2 = b(3)
            p3 = b(4)
            n1 = b(5)
            n2 = b(6)
            n3 = b(7)

            self.treatKine(tup, b, '_b')
            self.treatKine(tup, psi, '_psi')
            self.treatKine(tup, p1, '_p1')
            self.treatKine(tup, p2, '_p2')
            self.treatKine(tup, p3, '_p3')
            self.treatKine(tup, n1, '_n1')
            self.treatKine(tup, n2, '_n2')
            self.treatKine(tup, n3, '_n3')

            # add the information for Pid efficiency correction
            self.treatPions(tup, b)
            self.treatKaons(tup, b)
            self.treatMuons(tup, b)
            self.treatTracks(tup, b)

            tup.column('m', M(b) / GeV)
            tup.column('mass', self.mfit(b) / GeV)
            tup.column('c2dtf', self.c2dtf(b))
            tup.column('ctau', self.ctau(b))

            self._fill_doca_(b, tup)
            self._fill_mass_(b, tup)

            # add some reco-summary information
            self.addRecSummary(tup, rc_summary)

            # add the information needed for TisTos
            self.tisTos(psi, tup, 'psi_',
                        self.lines['psi'], self.l0tistos, self.tistos)

            tup.write()

        return SUCCESS

# ========================================================================
# @class  Psi3
#  Simple  template algorithm for B -> psi 3h
#  @date   2011-05-27
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru


class Psi3 (Psi5):

    """
    Simple  template algorithm for B -> psi 6
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
        # 1
        B = self.select('b', '[Beauty -> ( Meson -> mu+ mu- ) X+ X+ X-]CC')

        tup = self.nTuple('B')
        for b in B:
            #
            m = M(b) / GeV
            if not 5.0 < m < 6.7:
                continue
            #

            psi = b(1)

            mpsi = M(psi)
            if 4.5 * GeV < mpsi:
                continue
            elif 3.3 * GeV < mpsi:
                psi.setParticleID(LHCb.ParticleID(100443))

            p1 = b(2)
            p2 = b(3)
            n1 = b(4)

            self.treatKine(tup, b, '_b')
            self.treatKine(tup, psi, '_psi')
            self.treatKine(tup, p1, '_p1')
            self.treatKine(tup, p2, '_p2')
            self.treatKine(tup, n1, '_n1')

            # add the information for Pid efficiency correction
            self.treatPions(tup, b)
            self.treatKaons(tup, b)
            self.treatMuons(tup, b)
            self.treatTracks(tup, b)

            tup.column('m', M(b) / GeV)
            tup.column('mass', self.mfit(b) / GeV)
            tup.column('c2dtf', self.c2dtf(b))
            tup.column('ctau', self.ctau(b))

            self._fill_doca_(b, tup)
            self._fill_mass_(b, tup)

            # add some reco-summary information
            self.addRecSummary(tup, rc_summary)

            # add the information needed for TisTos
            self.tisTos(psi, tup, 'psi_',
                        self.lines['psi'], self.l0tistos, self.tistos)

            tup.write()

        return SUCCESS


# ========================================================================
# @class  Psi1
#  Simple  template algorithm for B -> psi K
#  @date   2011-05-27
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
class Psi1 (Psi5):

    """
    Simple  template algorithm for B -> psi K
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
        # 1
        B = self.select('b', '[Beauty -> ( Meson -> mu+ mu- ) X+ ]CC')

        tup = self.nTuple('B')
        for b in B:
            #
            m = M(b) / GeV
            if not 5.0 < m < 6.7:
                continue
            #

            psi = b(1)

            mpsi = M(psi)
            if 4.5 * GeV < mpsi:
                continue
            elif 3.3 * GeV < mpsi:
                psi.setParticleID(LHCb.ParticleID(100443))

            p1 = b(2)

            self.treatKine(tup, b, '_b')
            self.treatKine(tup, psi, '_psi')
            self.treatKine(tup, p1, '_h1')

            # add the information for Pid efficiency correction
            self.treatPions(tup, b)
            self.treatKaons(tup, b)
            self.treatMuons(tup, b)
            self.treatTracks(tup, b)

            tup.column('m', M(b) / GeV)
            tup.column('mass', self.mfit(b) / GeV)
            tup.column('c2dtf', self.c2dtf(b))
            tup.column('ctau', self.ctau(b))

            self._fill_doca_(b, tup)
            self._fill_mass_(b, tup)

            # add some reco-summary information
            self.addRecSummary(tup, rc_summary)

            # add the information needed for TisTos
            self.tisTos(psi, tup, 'psi_',
                        self.lines['psi'], self.l0tistos, self.tistos)

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

    the_year = "2012"

    if params:
        the_year = params['Year']
        logger.info('Year is set from params to be %s ' % the_year)
        logger.info('Year is set from files  to be %s ' % the_year)
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
        elif hasInFile(datafiles, 'Stripping19'):
            the_year = '2012'
        elif hasInFile(datafiles, 'Stripping20'):
            the_year = '2012'
        logger.info('Year is set from files  to be %s ' % the_year)

    #
    # check
    #
    if '2011' == the_year and hasInFile(datafiles, 'Collision12'):
        raise AttributeError, 'Invalid Year %s ' % the_year
    if '2012' == the_year and hasInFile(datafiles, 'Collision11'):
        raise AttributeError, 'Invalid Year %s ' % the_year

    rootInTES = '/Event/PSIX'

    from PhysConf.Filters import LoKi_Filters
    fltr5 = LoKi_Filters(
        VOID_Code="""
        (   0  < CONTAINS ( '/Event/PSIX/Phys/SelPsi5KPiForPsiX/Particles' ) ) 
        | ( 0  < CONTAINS ( '/Event/PSIX/Phys/SelPsi5KForPsiX/Particles'   ) ) 
        | ( 0  < CONTAINS ( '/Event/PSIX/Phys/SelPsi5PiForPsiX/Particles'  ) ) 
        """
    )
    fltr4 = LoKi_Filters(
        VOID_Code="""
        (   0  < CONTAINS ( '/Event/PSIX/Phys/SelPsi4KPiForPsiX/Particles' ) ) 
        | ( 0  < CONTAINS ( '/Event/PSIX/Phys/SelPsi4KForPsiX/Particles'   ) ) 
        | ( 0  < CONTAINS ( '/Event/PSIX/Phys/SelPsi4PiForPsiX/Particles'  ) ) 
        """
    )

    ## filters =fltr5.filters ('Filters')
    # filters.reserse()

    davinci = DaVinci(
        ### EventPreFilters = filters  ,
        DataType=the_year,
        InputType='MDST',
        Simulation=False,
        PrintFreq=1000,
        EvtMax=-1,
        #
        HistogramFile='PsiX_Histos.root',
        TupleFile='PsiX.root',
        #
        Lumi=True,
        #
    )

    #

    from Configurables import CondDB
    CondDB(LatestGlobalTagByDataType=the_year)

    # ------- decoding set-up start ----------
    from BenderTools.MicroDST import uDstConf
    uDstConf(rootInTES)
    # ------- decoding set-up end  -----------

    from Configurables import MessageSvc
    msg = MessageSvc()
    msg.setError += [
        'HcalDet.Quality',
        'EcalDet.Quality',
        'MagneticFieldSvc',
        'PropertyConfigSvc',
        'IntegrateBeamCrossing'
    ]

    # come back to Bender
    #
    setData(datafiles, catalogs, castor, castor)
    #
    # start Gaudi
    #
    gaudi = appMgr()

    alg3 = Psi3(
        'Psi3',  # Algorithm name ,
        RootInTES=rootInTES,
        Inputs=[
            'Phys/SelPsi3KForPsiX/Particles',
            'Phys/SelPsi3KPiForPsiX/Particles',
            'Phys/SelPsi3PiForPsiX/Particles',
        ]
    )

    alg4 = Psi4(
        'Psi4',  # Algorithm name ,
        RootInTES=rootInTES,
        Inputs=[
            'Phys/SelPsi4KForPsiX/Particles',
            'Phys/SelPsi4KPiForPsiX/Particles',
            'Phys/SelPsi4PiForPsiX/Particles',
        ]
    )

    alg5 = Psi5(
        'Psi5',  # Algorithm name ,
        RootInTES=rootInTES,
        Inputs=[
            'Phys/SelPsi5KForPsiX/Particles',
            'Phys/SelPsi5KPiForPsiX/Particles',
            'Phys/SelPsi5PiForPsiX/Particles',
        ]
    )

    alg6 = Psi6(
        'Psi6',  # Algorithm name ,
        RootInTES=rootInTES,
        Inputs=[
            'Phys/SelPsi6KForPsiX/Particles',
            'Phys/SelPsi6KPiForPsiX/Particles',
            'Phys/SelPsi6PiForPsiX/Particles',
        ]
    )

    alg1 = Psi1(
        'Psi1',  # Algorithm name ,
        RootInTES=rootInTES,
        Inputs=[
            'Phys/SelPsiKForPsiX/Particles',
        ]
    )

    mainSeq = gaudi.algorithm('GaudiSequencer/DaVinciUserSequence', True)
    mainSeq.Members += [
        #alg3 . name () ,
        alg4 . name(),
        #alg5 . name () ,
        #alg6 . name () ,
        #alg1 . name ()
    ]

    return SUCCESS

# =============================================================================
# The actual job steering
if '__main__' == __name__:

    input = [
        '/lhcb/LHCb/Collision12/PSIX.MDST/00021093/0000/00021093_000000%02d_1.psix.mdst' % i for i in range(1, 21)
    ]

    configure(input, castor=True)

    run(5000)

    gaudi = appMgr()

# =============================================================================
# The END
# =============================================================================
