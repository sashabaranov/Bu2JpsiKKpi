#!/usr/bin/env python
# ========================================================================
# @file DiCharm/CFake.py
#
#  Helper script to get fakes for popen chamr + W analysis
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

Helper script to get fakes for popen chamr + W analysis a

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
import BenderTools.Fill  # add methods to fill n-tuple info
# =============================================================================

# =============================================================================
# @class CFake
#  Helper script to get fakes for popen chamr + W analysis a
#  @date   2011-05-27
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class CFake (Algo):

    """
    Simple algorithm 
    """

    def initialize(self):

        sc = Algo.       initialize(self)
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
        self.etCone_ = None

        self.  fill_finalize()
        self.tisTos_finalize()

        return Algo.finalize(self)

    # the only one esential method:
    def analyse(self):
        """
        The only one essential method
        """

        #primaries = self.vselect ( 'PVs' , ISPRIMARY )

        # if primaries.empty() :
        #    return self.Warning ( 'No primary vertices are found', SUCCESS )

        #
        # rec summary
        rc_summary = self.get_(
            '/Event/Rec/Summary').summaryData()
        ew_counter = self.get_('/Event/Counters/CharmEW', False)
        #

        d0W = self.select('d0', '[ H_10 -> D0        [mu+]cc ]CC')
        dpW = self.select('dp', '[ H_10 -> D+        [mu+]cc ]CC')
        dsW = self.select('ds', '[ H_10 -> D_s+      [mu+]cc ]CC')
        lcW = self.select('lc', '[ H_10 -> Lambda_c+ [mu+]cc ]CC')

        tupD0 = self.nTuple('D0')
        tupDp = self.nTuple('Dp')
        tupDs = self.nTuple('Ds')
        tupLc = self.nTuple('Lc')

        for p in d0W:

            d0 = p(1)
            mu1 = p(2)

            self.treatKine(tupD0, d0, '_D0')
            self.treatKine(tupD0, mu1, '_mu1')

            tupD0.column_float('mass', self.mfit(d0) / GeV)
            tupD0.column_float('c2dtf', self.c2dtf(p))
            tupD0.column_float('c2dtfC', self.c2dtf(d0))
            tupD0.column_float('c2dtf_mu1', self.c2dtf(mu1))
            tupD0.column_float('ip_mu1', self.ip(mu1))
            tupD0.column_bool('ismuon_mu1', ISMUON(mu1))
            tupD0.column_bool('inmuon_mu1', INMUON(mu1))
            tupD0.column_bool('isloose_mu1', ISMUONLOOSE(mu1))

            pv = self.bestVertex(p)
            if not pv:
                self.Warning('Illegal primary vertex for D0/W!')
                continue

            ip2 = IP(self.geo(), pv)

            tupD0.column_float('ip2_mu1', ip2(mu1))
            tupD0.column_float('perr2_mu1', self.perr2(mu1))
            tupD0.column_float('ptCone_mu1', self.ptCone_(mu1))
            tupD0.column_float('etCone_mu1', self.etCone_(mu1))

            # add the information needed for TisTos
            self.tisTos(mu1, tupD0, 'W1_',
                        self.lines['W+'],
                        self.l0tistos,
                        self.tistos)

            # add the information for Pid efficiency correction
            self.treatMuons(tupD0,  p)
            self.treatKaons(tupD0, d0)
            self.treatPions(tupD0, d0)
            # add the infomation needed for track efficiency correction
            self.treatTracks(tupD0,  p)

            # add some reco-summary information
            self.addRecSummary(tupD0, rc_summary)
            # more counters
            tupD0.column_aux(ew_counter)

            #
            tupD0.write()

        for p in dpW:

            dp = p(1)
            mu1 = p(2)

            self.treatKine(tupDp, dp, '_Dp')
            self.treatKine(tupDp, mu1, '_mu1')

            tupDp.column_float('mass', self.mfit(dp) / GeV)
            tupDp.column_float('c2dtf', self.c2dtf(p))
            tupDp.column_float('c2dtfC', self.c2dtf(dp))
            tupDp.column_float('c2dtf_mu1', self.c2dtf(mu1))
            tupDp.column_float('ip_mu1', self.ip(mu1))
            tupDp.column_bool('ismuon_mu1', ISMUON(mu1))
            tupDp.column_bool('inmuon_mu1', INMUON(mu1))
            tupDp.column_bool('isloose_mu1', ISMUONLOOSE(mu1))

            pv = self.bestVertex(p)
            if not pv:
                self.Warning('Illegal primary vertex for D+/W!')
                continue

            ip2 = IP(self.geo(), pv)

            tupDp.column_float('ip2_mu1', ip2(mu1))
            tupDp.column_float('perr2_mu1', self.perr2(mu1))
            tupDp.column_float('ptCone_mu1', self.ptCone_(mu1))
            tupDp.column_float('etCone_mu1', self.etCone_(mu1))

            # add the information needed for TisTos
            self.tisTos(mu1, tupDp, 'W1_',
                        self.lines['W+'],
                        self.l0tistos,
                        self.tistos)

            # add the information for Pid efficiency correction
            self.treatMuons(tupDp,  p)
            self.treatKaons(tupDp, dp)
            self.treatPions(tupDp, dp)
            # add the infomation needed for track efficiency correction
            self.treatTracks(tupDp,  p)

            # add some reco-summary information
            self.addRecSummary(tupDp, rc_summary)
            # more counters
            tupDp.column_aux(ew_counter)
            #
            tupDp.write()

        for p in dsW:

            ds = p(1)
            mu1 = p(2)

            self.treatKine(tupDs, ds, '_Ds')
            self.treatKine(tupDs, mu1, '_mu1')

            tupDs.column_float('mass', self.mfit(ds) / GeV)
            tupDs.column_float('c2dtf', self.c2dtf(p))
            tupDs.column_float('c2dtfC', self.c2dtf(ds))
            tupDs.column_float('c2dtf_mu1', self.c2dtf(mu1))
            tupDs.column_float('ip_mu1', self.ip(mu1))
            tupDs.column_bool('ismuon_mu1', ISMUON(mu1))
            tupDs.column_bool('inmuon_mu1', INMUON(mu1))
            tupDs.column_bool('isloose_mu1', ISMUONLOOSE(mu1))

            pv = self.bestVertex(p)
            if not pv:
                self.Warning('Illegal primary vertex for D_s+/W!')
                continue

            ip2 = IP(self.geo(), pv)

            tupDs.column_float('ip2_mu1', ip2(mu1))
            tupDs.column_float('perr2_mu1', self.perr2(mu1))
            tupDs.column_float('ptCone_mu1', self.ptCone_(mu1))
            tupDs.column_float('etCone_mu1', self.etCone_(mu1))

            # add the information needed for TisTos
            self.tisTos(mu1, tupDs, 'W1_',
                        self.lines['W+'],
                        self.l0tistos,
                        self.tistos)

            # add the information for Pid efficiency correction
            self.treatMuons(tupDs,  p)
            self.treatKaons(tupDs, ds)
            self.treatPions(tupDs, ds)
            # add the infomation needed for track efficiency correction
            self.treatTracks(tupDs,  p)

            # add some reco-summary information
            self.addRecSummary(tupDs, rc_summary)
            # more counters
            tupDs.column_aux(ew_counter)
            #
            tupDs.write()

        for p in lcW:

            lc = p(1)
            mu1 = p(2)

            self.treatKine(tupLc, lc, '_Lc')
            self.treatKine(tupLc, mu1, '_mu1')

            tupLc.column_float('mass', self.mfit(lc) / GeV)
            tupLc.column_float('c2dtf', self.c2dtf(p))
            tupLc.column_float('c2dtfC', self.c2dtf(lc))
            tupLc.column_float('c2dtf_mu1', self.c2dtf(mu1))
            tupLc.column_float('ip_mu1', self.ip(mu1))
            tupLc.column_bool('ismuon_mu1', ISMUON(mu1))
            tupLc.column_bool('inmuon_mu1', INMUON(mu1))
            tupLc.column_bool('isloose_mu1', ISMUONLOOSE(mu1))

            pv = self.bestVertex(p)
            if not pv:
                self.Warning('Illegal primary vertex for Lambda_c+/W!')
                continue

            ip2 = IP(self.geo(), pv)

            tupLc.column_float('ip2_mu1', ip2(mu1))
            tupLc.column_float('perr2_mu1', self.perr2(mu1))
            tupLc.column_float('ptCone_mu1', self.ptCone_(mu1))
            tupLc.column_float('etCone_mu1', self.etCone_(mu1))

            # add the information needed for TisTos
            self.tisTos(mu1, tupLc, 'W1_',
                        self.lines['W+'],
                        self.l0tistos,
                        self.tistos)

            # add the information for Pid efficiency correction
            self.treatMuons(tupLc,  p)
            self.treatKaons(tupLc, lc)
            self.treatPions(tupLc, lc)
            self.treatProtons(tupLc, lc)
            # add the infomation needed for track efficiency correction
            self.treatTracks(tupLc,  p)

            # add some reco-summary information
            self.addRecSummary(tupLc, rc_summary)
            # more counters
            tupDs.column_aux(ew_counter)
            #
            tupLc.write()

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

    the_year = "2012"

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
        logger.info('Year is set from files  to be %s ' % the_year)

    #
    # check
    #
    if '2011' == the_year and hasInFile(datafiles, 'Collision12'):
        raise AttributeError, 'Invalid Year %s ' % the_year
    if '2012' == the_year and hasInFile(datafiles, 'Collision11'):
        raise AttributeError, 'Invalid Year %s ' % the_year

    logger.info('Use the Year = %s ' % the_year)

    rootInTES = '/Event/CharmFake'

    davinci = DaVinci(
        DataType=the_year,
        InputType='MDST',
        Simulation=False,
        PrintFreq=10000,
        EvtMax=-1,
        #
        HistogramFile='Fake_Histos.root',
        TupleFile='Fake.root',
        #
        Lumi=True,
        #
    )

    from Configurables import CondDB
    CondDB(LatestGlobalTagByDataType=the_year)

    #
    # take care abotu micro-dst
    #

    from BenderTools.MicroDST import uDstConf
    uDstConf(rootInTES)
    #
    # come back to Bender
    #
    setData(datafiles, catalogs, castor)

    #
    # start Gaudi
    #
    gaudi = appMgr()

    algF = CFake(
        'Fake',
        RootInTES=rootInTES,
        Inputs=['/Event/Phys/CharmAndW/Particles']
    )

    mainSeq = gaudi.algorithm('GaudiSequencer/DaVinciUserSequence', True)
    mainSeq.Members += [algF . name()]

    return SUCCESS

# =============================================================================
# The actual job steering
if '__main__' == __name__:

    import shelve
    db = shelve.open('/afs/cern.ch/user/i/ibelyaev/public/LFNs/lfns.db')

    input = db['D+F-mdst:2k+12/Down;May 8']
    input += db['D+F-mdst:2k+12/Up;May 8']

# input = [
# '/lhcb/user/i/ibelyaev/2013_05/48176/48176665/CFake.CharmFake.dst'
# 'CFake.CharmFake.mdst'
# ]

    configure(input, castor=True)

    run(500)

    gaudi = appMgr()

# =============================================================================
# The END
# =============================================================================
