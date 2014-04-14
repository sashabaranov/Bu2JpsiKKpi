#!/usr/bin/env ipython
# ========================================================================
# @file DiCharm/Bc2.py
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
# @class Bc2Bp
#  Simple template algorithm to study Charm & Y
#  @date   2011-05-27
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
class Bc2Bp (DiCharmAlg):

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

        bc = self.select('bc', '[ Beauty -> Beauty K- pi+]CC')

        c2dtf = DTF_CHI2NDOF(True, strings(['J/psi(1S)', 'B+']))
        mfit = DTF_FUN(M, True, strings(['J/psi(1S)', 'B+']))
        m1fit = DTF_FUN(M1, True, strings(['J/psi(1S)']))
        ctau = DTF_CTAU(0, True, strings(['J/psi(1S)', 'B+']))

        ctauC = DTF_CTAU(
            'B+' == ABSID, True, strings(['J/psi(1S)', 'B+']))
        ctauCs = DTF_CTAUSIGNIFICANCE(
            'B+' == ABSID, True, strings(['J/psi(1S)', 'B+']))

        miniph = MINTREE(('pi+' == ABSID) | ('K+' == ABSID), BPVIPCHI2())

        doca = DOCA(1, 2, self.distanceCalculator())

        tup = self.nTuple('B2B')
        for b in bc:
            #
            m = M(b) / GeV
            if not m < 7.1:
                continue
            #
            bp = b(1)
            kaon = b(2)
            pion = b(3)

            psi = bp.child(1)

            self.treatKine(tup, b, '_bc')
            self.treatKine(tup, bp, '_bp')
            self.treatKine(tup, kaon, '_k')
            self.treatKine(tup, pion, '_pi')

            psi = bp.child(1)
            bk = bp.child(2)

            self.treatKine(tup, psi, '_psi')
            self.treatKine(tup, bk, '_bk')

            tup.column('mass', mfit(b) / GeV)
            tup.column('m1', m1fit(b) / GeV)
            tup.column('c2dtf', c2dtf(b))
            tup.column('ctau', ctau(b))
            tup.column('ctauC', ctauC(b))
            tup.column('ctauCs', ctauCs(b))
            tup.column('minipH', miniph(b))
            tup.column('mKpi', M23(b) / GeV)

            # fictive K*
            p = LoKi.LorentzVector()
            p += kaon.momentum()
            p += pion.momentum()

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
                        self.lines['psi'], self.l0tistos, self.tistos)

            #
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

    logger.info('Use the Year = %s ' % the_year)

    rootInTES = '/Event/BC'

    davinci = DaVinci(
        DataType=the_year,
        InputType='MDST',
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

    from Configurables import CondDB
    CondDB(LatestGlobalTagByDataType=the_year)

    # ------- decoding set-up start ----------
    from BenderTools.MicroDST import uDstConf
    uDstConf(rootInTES)
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
    alg = Bc2Bp(
        'Bc2Bp',  # Algorithm name ,
        RootInTES=rootInTES,
        Inputs=['Phys/Bc/Particles']
    )

    mainSeq = gaudi.algorithm('GaudiSequencer/DaVinciUserSequence', True)
    mainSeq.Members += [alg . name()]

    return SUCCESS

# =============================================================================
# The actual job steering
if '__main__' == __name__:

    import shelve
    db = shelve.open('/afs/cern.ch/user/i/ibelyaev/public/LFNs/lfns.db')

    bc = db['Bc2DpK;Stripping20,MagDown']
    bc += db['Bc2DpK;Stripping20,MagUp']

    input = bc

    db.close()

    configure(input, castor=True)

    run(5000)

    gaudi = appMgr()

# =============================================================================
# The END
# =============================================================================
