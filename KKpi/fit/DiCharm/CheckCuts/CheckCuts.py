#!/usr/bin/env python
# =============================================================================
# @file
#
# Simple module to check/validate the stirpping cuts for charmed particles
#
# @author Vanya  BELYAEV Ivan.Belyaev@cern.ch
# @date   2011-06-19
#
#                   $Revision$
# Last modification $Date$
#                by $Author$
# =============================================================================
"""

Simple module to check/validate the stirpping cuts for charmed particles 


"""
# =============================================================================
__author__ = " Vanya BELYAEV Ivan.Belyaev@cern.ch "
__date__ = " 2011-06-19 "
__version__ = " $Revision: 1.2 $ "
# =============================================================================
# import everything from bender
from Bender.Main import *
from GaudiKernel.SystemOfUnits import GeV
from GaudiKernel.PhysicalConstants import c_light
# =============================================================================

# ========================================================================
# @class Dplus
#  Simple class to check striping cuts
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date 2011-06-17


class Dplus(Algo):

    """
    Simple class to check striping cuts
    """
    # standard method for analyses

    def initialize(self):
        """
        Initialization of algorithm
        """

        sc = Algo.initialize(self)
        if sc.isFailure():
            return sc

        ## self._decay = '[X -> ^( D+ --> K- pi+ pi+ )  Charm ]CC'
        self._decay = '[D+ --> K- pi+ pi+]CC'

        self._tupvars = {
            'mass':                       M / GeV,
            'pt':                       PT / GeV,
            'trChi2': MAXTREE(HASTRACK,  TRCHI2DOF),
            'minPT': MINTREE(HASTRACK,  PT) / GeV,
            #
            'DLLK': MINTREE('K+' == ABSID, PIDK - PIDpi),
            'DLLpi': MINTREE('pi+' == ABSID, PIDpi - PIDK),
            #
            'chi2vx': VFASPF(VCHI2),
            'ctau': BPVLTIME() * c_light,
            'ipchi2': BPVLTFITCHI2()
        }

        return SUCCESS

    # the standard method for analyse
    def analyse(self):
        """
        Standard method for analyses
        """
        # get particles
        parts = self.select('parts', self._decay)

        tup = self.nTuple('Cuts')

        for p in parts:

            for k, v in self._tupvars.iteritems():
                tup.column(k, v(p))

            tup.write()

        return SUCCESS                                           # RETURN

    # standard method for finalization
    def finalize(self):
        """
        print cuts
        """

        print 100 * '*'
        print ' DECAY  : ', self._decay
        print ' VARS   : ', self._tupvars.keys()
        print 100 * '*'

        self._tupvars = None

        return Algo.finalize(self)

# ========================================================================
# @class D0
#  Simple class to check stripping cuts for D0
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date 2011-06-17


class D0(Dplus):

    """
    Simple class to check stripping cuts for D0
    """
    # standard method for analyses

    def initialize(self):
        """
        Initialization of algorithm
        """

        sc = Dplus.initialize(self)
        if sc.isFailure():
            return sc

        self._decay = '[ D0 -> K- pi+ ]CC'
        ## self._decay = '[X -> ^( D0 -> K- pi+)  Charm]CC'

        return SUCCESS


# ========================================================================
# @class Dimu
#  Simple class to check stripping cuts for dimuons
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date 2011-06-17
class Dimu(Dplus):

    """
    Simple class to check stripping cuts for dimuons
    """
    # standard method for analyses

    def initialize(self):
        """
        Initialization of algorithm
        """
        sc = Dplus.initialize(self)
        if sc.isFailure():
            return sc

        self._decay = '[X -> ^( J/psi(1S) -> mu+ mu-)  Charm]CC'

        self._tupvars = {
            #
            'mass':            M / GeV,
            'pt':            PT / GeV,
            #
            'trChi2': MAXTREE(HASTRACK,  TRCHI2DOF),
            'minPT': MINTREE(HASTRACK,  PT) / GeV,
            #
            'DLLmu': MINTREE('mu+' == ABSID, PIDmu - PIDpi),
            #
            'chi2vx':             VFASPF(VCHI2),
        }

        return SUCCESS

# ========================================================================
# @class Dstar
#  Simple class to check stripping cuts for D*+
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date 2011-06-17


class Dstar(D0):

    """
    Simple class to check stripping cuts for D*+
    """
    # standard method for analyses

    def initialize(self):
        """
        Initialization of algorithm
        """

        sc = D0.initialize(self)
        if sc.isFailure():
            return sc

        ## self._decay = '[X -> ^( D*(2010)+ -> ( D0 -> K- pi+) pi+ ) Charm]CC'
        self._decay = '[D*(2010)+ -> ( D0 -> K- pi+) pi+]CC'

        self._tupvars = {
            #
            'mass':        M / GeV,
            'deltaM':      (M - M1) / GeV,
            'pt':            PT / GeV,
            'pt_D0': CHILD(1, PT) / GeV,
            #
            'trChi2_D0': CHILD(1, MAXTREE(HASTRACK,  TRCHI2DOF)),
            'trChi2_pi': CHILD(2,                       TRCHI2DOF),
            'minPT_D0': CHILD(1, MINTREE(HASTRACK,  PT)) / GeV,
            'minPT_pi': CHILD(2,                       PT) / GeV,
            #
            'DLLK': CHILD(1, MINTREE('K+' == ABSID, PIDK - PIDpi)),
            'DLLpi': CHILD(1, MINTREE('pi+' == ABSID, PIDpi - PIDK)),
            #
            'DLLpi_pi': CHILD(2, PIDpi - PIDK),
            #
            'chi2vx_D0': CHILD(1, VFASPF(VCHI2)),
            'chi2vx':             VFASPF(VCHI2),
            'ctau': CHILD(1, BPVLTIME()) * c_light,
            'ipchi2': CHILD(1, BPVLTFITCHI2())
        }

        return SUCCESS

# ========================================================================
# commmon configuration


def configure(datafiles, catalogs=[]):
    """
    Configure the job: common part
    """

    from PhysConf.Filters import LoKi_Filters
    fltrs = LoKi_Filters(
        VOID_Code="""
        ( 0 < CONTAINS ( '/Event/Charm/Phys/CharmAndDiMuonForPromptCharm/Particles') )
        |
        ( 0 < CONTAINS ( '/Event/Charm/Phys/DiCharmForPromptCharm/Particles'       ) ) 
        | 
        ( 0 < CONTAINS ( '/Event/Charm/Phys/D02HHForPromptCharm/Particles' ) )
        |
        ( 0 < CONTAINS ( '/Event/Charm/Phys/DForPromptCharm/Particles'     ) )
        |
        ( 0 < CONTAINS ( '/Event/Charm/Phys/DsForPromptCharm/Particles'    ) )
        """
    )
    filters = fltrs.filters('Filters')
    filters.reverse()

    from Configurables import DaVinci
    daVinci = DaVinci(
        DataType='2011',
        InputType='MDST',
        Simulation=False,
        PrintFreq=1000,
        #
        HistogramFile='CheckCuts_Histos.root',
        TupleFile='CheckCuts.root',
        #
        EventPreFilters=filters,
        Lumi=True
    )

    from Configurables import CondDB
    CondDB(UseLatestTags=["2011"],
           ##UseOracle     = True
           )

    # define/set the input data
    setData(datafiles, catalogs)

    # get the actual application manager (create if needed)
    gaudi = appMgr()

    # create local algorithm:
    dplus = Dplus(
        'Dplus',
        PropertiesPrint=True,
        HistoPrint=True,
        #
        RootInTES='/Event/Charm',
        Inputs=['Phys/DForPromptCharm/Particles']
        #
    )
    # create local algorithm:
    d0 = D0(
        'D0',
        PropertiesPrint=True,
        HistoPrint=True,
        #
        RootInTES='/Event/Charm',
        Inputs=['Phys/D02HHForPromptCharm/Particles']
        #
    )
    # create local algorithm:
    dstar = Dstar(
        'Dstar',
        PropertiesPrint=True,
        HistoPrint=True,
        #
        RootInTES='/Event/Charm',
        Inputs=['Phys/DstarForPromptCharm/Particles']
        #
    )
    # create local algorithm:
    dimu = Dimu(
        'Dimu',
        PropertiesPrint=True,
        HistoPrint=True,
        #
        RootInTES='/Event/Charm',
        Inputs=['Phys/CharmAndDiMuonForPromptCharm/Particles']
    )

    # finally inform Application Manager about our algorithm
    mainSeq = gaudi.algorithm('GaudiSequencer/DaVinciUserSequence', True)
    mainSeq.Members += [dplus  . name(),
                        d0     . name(),
                        dstar  . name(),
                        dimu   . name()]

    return SUCCESS

# =============================================================================
# job steering
if __name__ == '__main__':

    # make printout of the own documentations
    print '*' * 120
    print __doc__
    print ' Author  : %s ' % __author__
    print ' Version : %s ' % __version__
    print ' Date    : %s ' % __date__
    print '*' * 120

    # configure the job:
    # Stripping 13 Charm Micro-DST
    inputdata = [
        "/castor/cern.ch/grid/lhcb/LHCb/Collision11/CHARM.MDST/00010646/0000/00010646_000000%02d_1.charm.mdst" % i for i in range(1, 20)
    ]

    configure(inputdata)

    # run the job
    run(1000)

# =============================================================================
# The END
# =============================================================================
