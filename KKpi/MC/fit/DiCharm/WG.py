#!/usr/bin/env ipython
# ========================================================================
# @file WG.py
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
# @class  PsiX
#  Simple  template algorithm for B -> psi X
#  @date   2011-05-27
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class PsiX (Algo):

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

        B = self.select('b', 'Beauty -> ( Meson -> mu+ mu- ) ...')

        c2dtf = DTF_CHI2NDOF(True, strings(['J/psi(1S)', 'psi(2S)']))
        mfit = DTF_FUN(
            M, True, strings(['J/psi(1S)', 'psi(2S)']))
        ctau = DTF_CTAU(
            0, True, strings(['J/psi(1S)', 'psi(2S)']))

        maxTrChi2 = MAXTREE(ISBASIC & HASTRACK, TRCHI2DOF)
        maxTrGhost = MAXTREE(ISBASIC & HASTRACK, TRGHOSTPROB)

        minDllK = MINTREE('K+' == ABSID, PIDK - PIDpi)
        minDllPi = MINTREE('pi+' == ABSID, PIDpi - PIDK)

        nK = NINTREE('K+' == ABSID)
        nPi = NINTREE('pi+' == ABSID)

        tup = self.nTuple('B')
        for b in B:
            #
            m = M(b) / GeV
            if not 5.0 < m < 7.1:
                continue
            #

            psi = b(1)

            mpsi = M(psi)
            if 4.5 * GeV < mpsi:
                continue
            elif 3.3 * GeV < mpsi:
                psi.setParticleID(LHCb.ParticleID(100443))

            tup.column('m', M(b) / GeV)
            tup.column('mass', mfit(b) / GeV)
            tup.column('c2dtf', c2dtf(b))
            tup.column('ctau', ctau(b))

            tup.column('trChi2', maxTrChi2(b))
            tup.column('trGhPr', maxTrGhost(b))

            tup.column('mdimu',       M(psi) / GeV)
            tup.column('id_psi', int(ID(psi)))
            tup.column('id_b', int(ID(b)))

            tup.column('pid_K', minDllK(b))
            tup.column('pid_pi', minDllPi(b))

            tup.column('n_K', int(nK(b)))
            tup.column('n_Pi', int(nPi(b)))

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

    rootInTES = '/Event/PSIX'

    davinci = DaVinci(
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

    from Configurables import CondDB
    CondDB(LatestGlobalTagByDataType=the_year)

    # ------- decoding set-up start ----------
    from BenderTools.MicroDST import uDstConf
    uDstConf(rootInTES)
    # ------- decoding set-up end  -----------

    #
    # come back to Bender
    #
    setData(datafiles, catalogs, castor, castor)

    #
    # start Gaudi
    #
    gaudi = appMgr()

    alg = PsiX(
        'PsiX',  # Algorithm name ,
        RootInTES=rootInTES,
        Inputs=[
            'Phys/SelPsi%sForPsiX/Particles' % i for i in ['K',
                                                           'Pi',
                                                           #
                                                           '2K',
                                                           'KPi',
                                                           '2Pi',
                                                           #
                                                           '3K',
                                                           '3KPi',
                                                           '3Pi',
                                                           #
                                                           '4K',
                                                           '4KPi',
                                                           '4Pi',
                                                           #
                                                           '5K',
                                                           '5KPi',
                                                           '5Pi']
        ]
    )

    mainSeq = gaudi.algorithm('GaudiSequencer/DaVinciUserSequence', True)
    mainSeq.Members += [alg . name()]

    return SUCCESS

# =============================================================================
# The actual job steering
if '__main__' == __name__:

    input = [
        '/lhcb/LHCb/Collision12/PSIX.MDST/00020584/0000/00020584_00000003_1.psix.mdst'
    ]

    configure(input, castor=True)

    run(5000)

    gaudi = appMgr()

# =============================================================================
# The END
# =============================================================================
