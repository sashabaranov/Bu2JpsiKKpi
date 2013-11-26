#!/usr/bin/env python
# =============================================================================
# @file
#
# Simple script to get the generator cuts acceptance for prompt charm particles
#
# The actual code has been kindly provided by Greig Cowan
#
# @author Vanya  BELYAEV Ivan.Belyaev@cern.ch
# @date   2011-06-18
#
#                   $Revision$
# Last modification $Date$
#                by $Author$
# =============================================================================
"""

Simple script to get the generator cuts acceptance for D+

The actual code has been kindly provided by Greig Cowan

"""
# =============================================================================
__author__ = " Vanya BELYAEV Ivan.Belyaev@cern.ch "
__date__ = " 2011-06-18 "
__version__ = " $Revision: 1.2 $ "
# =============================================================================
# import everything from bender
from Bender.MainMC import *
from GaudiKernel.SystemOfUnits import GeV
from GaudiKernel.PhysicalConstants import c_light
# =============================================================================
HepMC.GenParticle.Vector = std.vector('HepMC::GenParticle*')

# ========================================================================
# @class Dplus
#  Simple class to check generator cuts efficiency
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date 2011-06-17


class Dplus(AlgoMC):

    """
    Simple class to check Generator Level Cuts for D+ -> K- pi+ pi+ decay
    """
    # standard method for analyses

    def initialize(self):
        """
        Initialization of algorithm
        """

        sc = AlgoMC.initialize(self)
        if sc.isFailure():
            return sc

        self.inLHCb = self.tool(cpp.IGenCutTool, 'DaughtersInLHCb')
        self._part = ('D+' == GABSID) & in_range(
            2.0, GY, 4.5)
        self._part = self._part & in_range(
            0 * GeV, GPT, 12 * GeV)

        self._decay = "[ D+ ==> K- pi+ pi+ ]CC"

        self._cnt = self.counter("multiplicity")

        self._gv = HepMC.GenEvent()
        self._gc = LHCb.GenCollision()

        return SUCCESS

    # the standard method for analyse
    def analyse(self):
        """
        Standard method for analyses
        """
        # skip B-decays
        gb = self.gselect('beauty', GBEAUTY)

        gen0 = self.gselect('gen0',   self._decay)
        if gen0.empty():
            return self.Warning('No proper MC-decays are found', SUCCESS)

        gen1 = self.gselect('gen1', gen0, self._part)
        if gen1.empty():
            return SUCCESS

        self._cnt += len(gen1)

        self.setFilterPassed(not gen1.empty())
        tup = self.nTuple('Gen')

        n = 0
        for p in gen1:
            v = HepMC.GenParticle.Vector()
            v.push_back(p)
            n += 1

            tup.column('y',  GY(p))
            tup.column('pt', GPT(p) / GeV)
            tup.column('B', not gb.empty(), 0, 1)
            tup.column('n', n, 1, 7)
            acc = self.inLHCb.applyCut(v, self._gv, self._gc)
            tup.column('acc', acc, 0, 1)

            tup.write()

        return SUCCESS                                           # RETURN

    # standard method for finalization
    def finalize(self):
        """
        print cuts
        """

        print 100 * '*'

        print 'DECAY  : ', self._decay
        print 'CUTS   : ', self._part

        print 100 * '*'

        self.inLHCb = None

        return AlgoMC.finalize(self)

# ========================================================================
# @class Dstar
#  Simple class to check generator cuts efficiency for D*+
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date 2011-06-17


class Dstar(Dplus):

    """
    Simple class to check Generator Level Cuts for D*+ -> ( D0 -> K pi ) pi+ decays
    """
    # standard method for analyses

    def initialize(self):
        """
        Initialization of algorithm
        """

        sc = Dplus.initialize(self)
        if sc.isFailure():
            return sc

        self._part = ('D*(2010)+' == GABSID) & in_range(
            2.0, GY, 4.5)
        self._part = self._part & in_range(
            0 * GeV, GPT, 12 * GeV)

        self._decay = "[ D*(2010)+ -> ( D0 => K- pi+) pi+ ]CC"

        return SUCCESS

# ========================================================================
# @class Ds
#  Simple class to check generator cuts efficiency for Ds
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date 2011-06-17


class Ds(Dplus):

    """
    Simple class to check Generator Level Cuts for Ds+ -> K+ K- pi+ decays
    """
    # standard method for analyses

    def initialize(self):
        """
        Initialization of algorithm
        """

        sc = Dplus.initialize(self)
        if sc.isFailure():
            return sc

        self._part = ('D_s+' == GABSID) & in_range(
            2.0, GY, 4.5)
        self._part = self._part & in_range(
            0 * GeV, GPT, 12 * GeV)

        self._decay = "[ D_s+ ==> K+ K- pi+ ]CC"

        return SUCCESS

# ========================================================================
# @class Lc
#  Simple class to check generator cuts efficiency for Lc
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date 2011-06-17


class Lc(Dplus):

    """
    Simple class to check Generator Level Cuts for Lc+ -> p+ K- pi+ decays
    """
    # standard method for analyses

    def initialize(self):
        """
        Initialization of algorithm
        """

        sc = Dplus.initialize(self)
        if sc.isFailure():
            return sc

        self._part = ('Lambda_c+' == GABSID) & in_range(
            2.0, GY, 4.5)
        self._part = self._part & in_range(
            0 * GeV, GPT, 12 * GeV)

        self._decay = "[ Lambda_c+ ==> p+ K- pi+ ]CC"

        return SUCCESS

# ========================================================================
# @class D0
#  Simple class to check generator cuts efficiency for D0 form D*+ decays
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date 2011-06-17


class D0(Dstar):

    """
    Simple class to check Generator Level Cuts for D*+ -> ( D0 -> K pi ) pi+ decays
    """
    # standard method for analyses

    def initialize(self):
        """
        Initialization of algorithm
        """

        sc = Dstar.initialize(self)
        if sc.isFailure():
            return sc

        self._part = ('D0' == GABSID) & in_range(
            2.0, GY, 4.5)
        self._part = self._part & in_range(
            0 * GeV, GPT, 12 * GeV)

        self._decay = "[ D*(2010)+ -> ^( D0 => K- pi+) pi+ ]CC"

        return SUCCESS

    # the standard method for analyse
    def analyse(self):
        """
        Standard method for analyses
        """
        # skip B-decays
        gb = self.gselect('beauty', GBEAUTY)

        gen0 = self.gselect('gen0',   self._decay)
        if gen0.empty():
            return self.Warning('No proper MC-decays are found', SUCCESS)

        gen1 = self.gselect('gen1', gen0, self._part)
        if gen1.empty():
            return SUCCESS

        self._cnt += len(gen1)

        self.setFilterPassed(not gen1.empty())
        tup = self.nTuple('Gen')

        n = 0
        dstar = 'D*(2010)+' == GABSID
        for p in gen1:

            v = HepMC.GenParticle.Vector()
            m = None
            for a in p.ancestors():
                if dstar(a):
                    m = a
                    break
            if not m:
                continue

            v.push_back(m)
            n += 1

            tup.column('y',  GY(p))
            tup.column('pt', GPT(p) / GeV)
            tup.column('B', not gb.empty(), 0, 1)
            tup.column('n', n, 1, 7)
            acc = self.inLHCb.applyCut(v, self._gv, self._gc)
            tup.column('acc', acc, 0, 1)

            tup.write()

        return SUCCESS                                           # RETURN


# ========================================================================
# @class Jpsi
#  Simple class to check generator cuts efficiency for Jpsi
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date 2011-06-17
class Jpsi(Dplus):

    """
    Simple class to check Generator Level Cuts for Jpsi
    """
    # standard method for analyses

    def initialize(self):
        """
        Initialization of algorithm
        """

        sc = Dplus.initialize(self)
        if sc.isFailure():
            return sc

        self._part = ('J/psi(1S)' == GABSID) & in_range(
            2.0, GY, 4.5)
        self._part = self._part & in_range(
            0 * GeV, GPT, 12 * GeV)

        self._decay = "J/psi(1S) => mu+ mu-"

        return SUCCESS


# ========================================================================
# commmon configuration
def configure_alg(alg_type, datafiles, catalogs=[]):
    """
    Configure the job: common part
    """
    from Configurables import DaVinci
    daVinci = DaVinci(
        DataType='2010',
        Simulation=True
    )

    from Configurables import EventClockSvc
    EventClockSvc(
        EventTimeDecoder="FakeEventTime"
    )

    from Configurables import HistogramPersistencySvc
    HistogramPersistencySvc(OutputFile='generator_histos.root')

    from Configurables import NTupleSvc
    NTupleSvc(Output=[
        "CHECK DATAFILE='generator_tuple.root' TYPE='ROOT' OPT='NEW'"]
    )

    # define the proper postconfig action
    def postconfig():

        from Configurables import EventClockSvc, CondDB
        EventClockSvc(
            EventTimeDecoder="FakeEventTime"
        )
        CondDB(IgnoreHeartBeat=True)

    # Important: use Post Config action!
    from Gaudi.Configuration import appendPostConfigAction
    appendPostConfigAction(postconfig)

    # define/set the input data
    setData(datafiles, catalogs)

    # get the actual application manager (create if needed)
    gaudi = appMgr()

    # create local algorithm:
    alg = alg_type(
        'Acceptance',
        PropertiesPrint=True,
        HistoPrint=True,
        NTupleLUN="CHECK",
        PP2MCs=[],
        InputPrimaryVertices="None"
    )

    # finally inform Application Manager about our algorithm
    gaudi.setAlgorithms([
        alg
    ])

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
    inputdata = [
        # '/castor/cern.ch/grid/lhcb/user/i/ibelyaev/2011_06/21695/21695982/Gauss-21263010-1000ev-20110617.sim'
        # '/castor/cern.ch/grid/lhcb/user/i/ibelyaev/2011_06/21695/21695982/Gauss-21263010-1000ev-20110617.sim'
        # D*?
        '/castor/cern.ch/grid/lhcb/user/i/ibelyaev/2011_06/21691/21691479/Gauss-27163003-1000ev-20110617.sim'
    ]

    configure_alg(Dplus, inputdata)

    # run the job
    run(100)

# =============================================================================
# The END
# =============================================================================
