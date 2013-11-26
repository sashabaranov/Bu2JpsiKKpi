#!/usr/bin/env ipython
# =============================================================================
# @file DiCharm/Alg.py
#
#  Get from MC the background from B-decays
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

Get from MC the background from B-decays 

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
#
from GaudiKernel.SystemOfUnits import GeV, MeV, mm
from GaudiKernel.PhysicalConstants import c_light
import math
# =============================================================================


def dphi(c1, c2):

    d = float(MCPHI(c1) - MCPHI(c2))
    d /= math.pi

    while d > 1.:
        d -= 2.
    while d < -1.:
        d += 2.

    return d


def bkey(part):

    if BEAUTY(part):
        return part.key()
    mother = part.mother()
    while mother:
        if BEAUTY(mother):
            return mother.key()
        mother = mother.mother()
    return -1

import DiCharm.Efficiency as Eff
# =============================================================================
# @class MCB
#  Simple template algorithm to get background from B-decays
#  @date   2011-05-27
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch


class MCB2(AlgoMC):

    """
    Simple algorithm to study feeddown from B-decays 
    """
    # the only one esential method:

    def analyse(self):
        """
        The only one essential method
        """

        from_B = MCINANCESTORS(BEAUTY)
        pt_d_cut = in_range(2 * GeV, MCPT, 12 * GeV)
        y_d_cut = in_range(2, MCY, 4.0)
        pt_j_cut = MCPT < 12 * GeV
        y_j_cut = in_range(2, MCY, 4.5)
        # MCDECTREE ( '[ D0         -> K- pi+     ]CC') " ,
        d0 = 'D0' == MCABSID
        # MCDECTREE ( '[ D+        --> K- pi+ pi+ ]CC') " ,
        dp = 'D+' == MCABSID
        # MCDECTREE ( '[ D_s+      --> K- K+  pi+ ]CC') " ,
        ds = 'D_s+' == MCABSID
        # MCDECTREE ( '[ Lambda_c+ --> p+ K-  pi+ ]CC') " ,
        lc = 'Lambda_c+' == MCABSID
        # MCDECTREE ( '  J/psi(1S)  => mu+ mu-'       ) " ,
        psi = 'J/psi(1S)' == MCABSID
        psi_B = psi & from_B
        charm = (d0 | dp | ds | lc)
        charm_B = charm & from_B

        pv_key = MCVPXFUN(MCVKEY)

        charm = self. mcselect('charm', charm_B | psi_B)

        tup1 = self.nTuple('C1')
        tup = self.nTuple('C2')

        len_charm = len(charm)

        for i in range(0, len_charm):

            # loop over the second charm
            for j in range(i + 1, len_charm):

                c1 = charm[i]
                c2 = charm[j]

                if psi(c2):
                    c1 = charm[j]
                    c2 = charm[i]

                id1 = int(MCID(c1))
                pt1 = MCPT(c1) / GeV
                y1 = MCY(c1)

                id2 = int(MCID(c2))
                pt2 = MCPT(c2) / GeV
                y2 = MCY(c2)

                key1 = int(pv_key(c1))
                key2 = int(pv_key(c2))

                tup.column('pid1',  id1)
                tup.column('pid2',  id2)

                tup.column('pv1', key1)
                tup.column('pv2', key2)

                tup.column('bk1', bkey(c1))
                tup.column('bk2', bkey(c2))

                tup.column('pt1', pt1)
                tup.column('pt2', pt2)
                tup.column('y1', y1)
                tup.column('y2', y2)

                trg1 = Eff.VE(0.0, 0.0)
                rec1 = Eff.VE(0.0, 0.0)
                ok1 = False
                if psi(c1):
                    trg1 = Eff.tosEff_Jpsi(pt1, y1, 'L0xL1xL2')
                    rec1 = Eff.selEff_Jpsi_2D(pt1, y1)
                    rec1 *= Eff.br_Jpsi
                    ok1 = 0 <= pt1 <= 12 and 2 <= y1 < 4.5
                elif d0(c1):
                    trg1 = Eff.tosEff_D0(pt1, y1, 'L0xL1xL2')
                    rec1 = Eff.selEff_D0(pt1, y1)
                    rec1 *= Eff.br_D0
                    ok1 = 2 <= pt1 <= 12 and 2 <= y1 < 4.0
                elif dp(c1):
                    trg1 = Eff.tosEff_Dp(pt1, y1, 'L0xL1xL2')
                    rec1 = Eff.selEff_Dp(pt1, y1)
                    rec1 *= Eff.br_Dp
                    ok1 = 2 <= pt1 <= 12 and 2 <= y1 < 4.0
                elif ds(c1):
                    trg1 = Eff.tosEff_Ds(pt1, y1, 'L0xL1xL2')
                    rec1 = Eff.selEff_Ds(pt1, y1)
                    rec1 *= Eff.br_Ds
                    ok1 = 2 <= pt1 <= 12 and 2 <= y1 < 4.0
                elif lc(c1):
                    trg1 = Eff.tosEff_Lc(pt1, y1, 'L0xL1xL2')
                    rec1 = Eff.selEff_Lc(pt1, y1)
                    rec1 *= Eff.br_Lc
                    ok1 = 2 <= pt1 <= 12 and 2 <= y1 < 4.0

                trg2 = Eff.VE(0.0, 0.0)
                rec2 = Eff.VE(0.0, 0.0)
                ok2 = False
                if psi(c2):
                    trg2 = Eff.tosEff_Jpsi(pt2, y2, 'L0xL1xL2')
                    rec2 = Eff.selEff_Jpsi_2D(pt2, y2)
                    rec2 *= Eff.br_Jpsi
                    ok2 = 0 <= pt2 <= 12 and 2 <= y2 < 4.5
                elif d0(c2):
                    trg2 = Eff.tosEff_D0(pt2, y2, 'L0xL1xL2')
                    rec2 = Eff.selEff_D0(pt2, y2)
                    rec2 *= Eff.br_D0
                    ok2 = 2 <= pt2 <= 12 and 2 <= y2 < 4.0
                elif dp(c2):
                    trg2 = Eff.tosEff_Dp(pt2, y2, 'L0xL1xL2')
                    rec2 = Eff.selEff_Dp(pt2, y2)
                    rec2 *= Eff.br_Dp
                    ok2 = 2 <= pt2 <= 12 and 2 <= y2 < 4.0
                elif ds(c2):
                    trg2 = Eff.tosEff_Ds(pt2, y2, 'L0xL1xL2')
                    rec2 = Eff.selEff_Ds(pt2, y2)
                    rec2 *= Eff.br_Ds
                    ok2 = 2 <= pt2 <= 12 and 2 <= y2 < 4.0
                elif lc(c2):
                    trg2 = Eff.tosEff_Lc(pt2, y2, 'L0xL1xL2')
                    rec2 = Eff.selEff_Lc(pt2, y2)
                    rec2 *= Eff.br_Lc
                    ok2 = 2 <= pt2 <= 12 and 2 <= y2 < 4.0

                tup.column('trg1', trg1 . value())
                tup.column('trg2', trg2 . value())
                tup.column('rec1', rec1 . value())
                tup.column('rec2', rec2 . value())

                tup.column('ok1', ok1, 0, 1)
                tup.column('ok2', ok2, 0, 1)

                v2 = c1.momentum() + c2.momentum()

                tup.column('m2c', v2.M() / GeV)

                tup.column('cc', id1 * id2 > 0, 0, 1)
                tup.column('psi', psi(c1) or psi(c2), 0, 1)

                tup.column('dphi', dphi(c1, c2))

                tup.write()

        self.setFilterPassed(False)

        return SUCCESS

# =============================================================================
# configure the job


def configure(datafiles, catalogs=[]):
    """
    Job configuration 
    """

    ## needed for job configuration
    from Configurables import DaVinci
    ## needed for job configuration
    from Configurables import EventSelector
    ## needed for job configuration
    from GaudiConf.Configuration import FileCatalog
    ## needed for job configuration
    from GaudiConf.Configuration import NTupleSvc

    from PhysConf.Filters import LoKi_Filters

    fltrs = LoKi_Filters(
        #
        MC_Preambulo=[
            "from LoKiCore.functions import in_range, count, has ",
            "from GaudiKernel.SystemOfUnits import GeV",
            "from_B   = MCINANCESTORS ( BEAUTY )",
            "pt_d_cut = in_range ( 2 * GeV , MCPT , 12 * GeV )  ",
            "y_d_cut  = in_range ( 2       , MCY  , 4.0      )  ",
            "pt_j_cut =                      MCPT < 12 * GeV    ",
            "y_j_cut  = in_range ( 2       , MCY  , 4.5      )  ",
            "d0       = 'D0'        == MCABSID ",
            # MCDECTREE ( '[ D0         -> K- pi+     ]CC') " ,
            "dp       = 'D+'        == MCABSID ",
            # MCDECTREE ( '[ D+        --> K- pi+ pi+ ]CC') " ,
            "ds       = 'D_s+'      == MCABSID ",
            # MCDECTREE ( '[ D_s+      --> K- K+  pi+ ]CC') " ,
            "lc       = 'Lambda_c+' == MCABSID ",
            # MCDECTREE ( '[ Lambda_c+ --> p+ K-  pi+ ]CC') " ,
            # MCDECTREE ( '  J/psi(1S)  => mu+ mu-'       ) " ,
            "psi      = 'J/psi(1S)' == MCABSID ",
            "psi_B    = psi   & from_B ",
            "charm    = ( d0 | dp | ds | lc ) ",
            "charm_B  = charm & from_B ",
        ],
        #
        # di-charm mode
        #
# MC_Code    = """
        ##         has ( BEAUTY )
        # &
        ##         ( 1.5 < count ( y_d_cut & pt_d_cut & charm_B ) )
# """
        #
        # dimuon & charm
        #
# MC_Code    = """
        ##         has ( BEAUTY )
        # &
        ##         ( 0.5 < count ( y_d_cut & pt_d_cut & charm_B ) )
        # &
        ##         ( 0.5 < count ( y_j_cut & pt_j_cut & psi_B   ) )
# """
        #
        # alltogather
        #
        MC_Code="""
        has ( BEAUTY )
        &
        ( 1.5 < count (
        ( y_d_cut & pt_d_cut & charm_B ) |
        ( y_j_cut & pt_j_cut &   psi_B ) ) ) 
        """
    )

    filters = fltrs.filters('Filters')
    filters.reverse()

    davinci = DaVinci(
        DataType='2011',
        InputType='DST',
        Simulation=True,
        PrintFreq=1000,
        #EventPreFilters = filters ,
        EvtMax=-1,
        #
        HistogramFile='2xCharm_fromB_Histos.root',
        TupleFile='2xCharm_fromB.root',
        #
        Lumi=False
        #
    )

    # come back to Bender
    setData(datafiles, catalogs)

    gaudi = appMgr()

    alg = MCB2(
        #
        'MCB2',  # Algorithm name
        #
        Inputs=[],
        #
        PP2MCs=[]
    )
    #
    # mainSeq = gaudi.algorithm ('GaudiSequencer/DaVinciUserSequence', True )
    # mainSeq.Members += [ alg.name() ]
    gaudi.setAlgorithms([alg])

    return SUCCESS
# =============================================================================
# The actual job steering
if '__main__' == __name__:

    files = [
        # J/psi->mu+mu-
        "/castor/cern.ch/grid/lhcb/MC/MC10/ALLSTREAMS.DST/00009117/0000/00009117_00000%03d_1.allstreams.dst" % i for i in range(1, 720)
        # [ D*+ -> ( D0 => K- pi+ ) pi+ ]CC
        #"/castor/cern.ch/grid/lhcb/MC/MC10/ALLSTREAMS.DST/00008554/0000/00008554_00000%03d_1.allstreams.dst" % i for i in range(1,300)
        # D+ -> K- pi+ pi+  21263010
        #"/castor/cern.ch/grid/lhcb/MC/MC10/ALLSTREAMS.DST/00009404/0000/00009404_00000%03d_1.allstreams.dst" % i for i in range(1,400)
        # D_s+ -> K+ K- pi+   23263001
        #"/castor/cern.ch/grid/lhcb/MC/MC10/ALLSTREAMS.DST/00011124/0000/00011124_00000%03d_1.allstreams.dst" % i for i in range(1,600)
        # Lambda_c+ -> p K- pi+ 25103000
        #"/castor/cern.ch/grid/lhcb/MC/MC10/ALLSTREAMS.DST/00011118/0000/00011118_00000%03d_1.allstreams.dst" % i for i in range(1,255)
    ]

    configure(files)

    run(1000)

    gaudi = appMgr()

# =============================================================================
# The END
# =============================================================================
