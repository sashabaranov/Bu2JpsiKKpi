#!/usr/bin/env python
# ========================================================================
# @file DiCharm/Z2mumu.py
#
#  Helper script to study Z0 -> mumu kinematic and properties of high-pt muons
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

Helper script to study Z0 -> mumu kinematic and properties of high-pt muons 

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
# @class Z2mumu
#  Helper script to study Z0 -> mumu kinematic and properties of high-pt muons
#  @date   2011-05-27
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class Z2mumu (Algo):

    """
    Simple algorithm 
    """
    # standard constructor

    def __init__(self, name, **args):
        """
        Standard constructor
        """
        Algo.__init__(self, name, **args)
        #

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

        pion_cuts = in_range(300 * MeV, PT, 10 * GeV) & (CLONEDIST > 5000) & (
            TRCHI2DOF < 5) & (TRGHOSTPROB < 0.5) & (PERR2 / P2 < 0.05 ** 2)
        self.ptCone_ = SUMCONE(
            0.25, PT, '/Event/Phys/StdAllLoosePions/Particles')
        self.ptCone_2 = SUMCONE(
            0.25, PT, '/Event/Phys/StdAllLoosePions/Particles', pion_cuts)
        self.etCone_ = SUMCONE(
            0.25, PT, '/Event/Phys/StdLooseAllPhotons/Particles')

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

        primaries = self.vselect('PVs', ISPRIMARY)

        if primaries.empty():
            return self.Warning('No primary vertices are found', SUCCESS)

        #
        # rec summary
        rc_summary = self.get('/Event/Rec/Summary').summaryData()
        #

        Z = self.select('Z', 'Z0 -> mu+ mu-')

        tupZ = self.nTuple('Z')

        for z in Z:
            #
            mZ = M(z)
            if mZ < 40 * GeV:
                continue

            # print z.decay()

            mu1 = z(1)
            mu2 = z(2)

            self.decisions(z, self.triggers['Z0'])
            self.decisions(mu1, self.triggers['W+'])
            self.decisions(mu2, self.triggers['W-'])

            self.treatKine(tupZ, z, '_Z')
            self.treatKine(tupZ, mu1, '_mu1')
            self.treatKine    ( tupZ , mu2   , '_mu2'   )\

            tupZ.column_float('mass', self.mfit(z) / GeV)
            tupZ.column_float('c2dtf', self.c2dtf(z))

            tupZ.column_float('c2dtf_mu1', self.c2dtf(mu1))
            tupZ.column_float('c2dtf_mu2', self.c2dtf(mu2))

            tupZ.column_float('ip_mu1', self.ip(mu1))
            tupZ.column_float('ip_mu2', self.ip(mu2))

            pv = self.bestVertex(z)
            if not pv:
                self.Warning('Illegal primary vertex for Z0!')
                continue

            ip2 = IP(self.geo(), pv)

            tupZ.column_float('ip2_mu1', ip2(mu1))
            tupZ.column_float('ip2_mu2', ip2(mu2))

            tupZ.column_float('perr2_mu1', self.perr2(mu1))
            tupZ.column_float('perr2_mu2', self.perr2(mu2))

            tupZ.column_float('ptCone1_mu1', self.ptCone_(mu1))
            tupZ.column_float('ptCone2_mu1', self.ptCone_2(mu1))
            tupZ.column_float('etCone_mu1', self.etCone_(mu1))
            tupZ.column_float('ptCone1_mu2', self.ptCone_(mu2))
            tupZ.column_float('ptCone2_mu2', self.ptCone_2(mu2))
            tupZ.column_float('etCone_mu2', self.etCone_(mu2))

            # add the information needed for TisTos
            self.tisTos(z, tupZ, 'Z_',
                        self.lines['Z0'],
                        self.l0tistos,
                        self.tistos)

            # add the information needed for TisTos
            self.tisTos(mu1, tupZ, 'W1_',
                        self.lines['W+'],
                        self.l0tistos,
                        self.tistos)

            self.tisTos(mu2, tupZ, 'W2_',
                        self.lines['W-'],
                        self.l0tistos,
                        self.tistos)

            # add the information for Pid efficiency correction
            self.treatMuons(tupZ, z)
            # add the infomation needed for track efficiency correction
            self.treatTracks(tupZ, z)

            # add some reco-summary information
            self.addRecSummary(tupZ, rc_summary)
            #
            tupZ.write()

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

    Z_Location = '/Event/EW/Phys/Z02MuMuLine/Particles'
    from PhysSelPython.Wrappers import AutomaticData
    Z_Data = AutomaticData(Location=Z_Location)

    # Read only fired events to speed up
    from PhysConf.Filters import LoKi_Filters
    fltrs = LoKi_Filters(
        STRIP_Code="HLT_PASS_RE('Stripping.*Z02MuMu*.*')"
    )

    davinci = DaVinci(
        EventPreFilters=fltrs.filters('Filters'),
        DataType=the_year,
        InputType='DST',
        Simulation=False,
        PrintFreq=10000,
        EvtMax=-1,
        #
        HistogramFile='Z0_Histos.root',
        TupleFile='Z0.root',
        #
        Lumi=True,
        #
    )

    from Configurables import CondDB
    CondDB(LatestGlobalTagByDataType=the_year)

    #
    # come back to Bender
    #
    setData(datafiles, catalogs, castor)

    #
    # start Gaudi
    #
    gaudi = appMgr()

    algZ = Z2mumu(
        'Z0',
        Inputs=[Z_Location]
    )

    mainSeq = gaudi.algorithm('GaudiSequencer/DaVinciUserSequence', True)
    mainSeq.Members += [algZ . name()]

    return SUCCESS

# =============================================================================
# The actual job steering
if '__main__' == __name__:

    input = [
        '/lhcb/LHCb/Collision12/EW.DST/00020241/0000/00020241_00000%03d_1.ew.dst' % i for i in range(106, 500)
    ]

    configure(input, castor=True)

    run(50000)

    gaudi = appMgr()

# =============================================================================
# The END
# =============================================================================
