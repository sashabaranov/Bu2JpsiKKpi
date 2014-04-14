#!/usr/bin/env python

__author__ = "Alexander Baranov a.baranov@cern.ch"
__date__ = "Summer'13"
__version__ = " "

from Bender.All import *
from Gaudi.Configuration import *
from GaudiKernel.SystemOfUnits import GeV, MeV, mm, micrometer
from GaudiKernel.PhysicalConstants import c_light
from math import *
from LoKiTracks.decorators import *  # needed for TrKEY work

import BenderTools.Fill
import BenderTools.TisTos


class B2JpsiKKpi(Algo):

    def initialize(self):

        triggers = {}
        triggers ['psi'] = {}

        lines = {}
        lines [ "psi"      ] = {}
        lines [ "psi1"     ] = {} ## clean-selection: dimuon
        lines [ "psi2"     ] = {} ## clean-selection: detached
        lines [ "psi3"     ] = {} ## clean-selection: unbiased

        #
        ## J/psi
        #
        lines [ "psi" ][   'L0TOS' ] = 'L0(DiMuon|Muon)Decision'
        lines [ "psi" ][   'L0TIS' ] = 'L0(Hadron|DiMuon|Muon|Electron|Photon)Decision'
        lines [ "psi" ][ 'Hlt1TOS' ] = 'Hlt1(DiMuon|SingleMuon|TrackMuon).*Decision'
        lines [ "psi" ][ 'Hlt1TIS' ] = 'Hlt1(DiMuon|SingleMuon|Track).*Decision'
        lines [ "psi" ][ 'Hlt2TOS' ] = 'Hlt2(DiMuon|ExpressJPsi|SingleMuon).*Decision'
        lines [ "psi" ][ 'Hlt2TIS' ] = 'Hlt2(Charm|Topo|DiMuon|Single).*Decision'
        #
        ## J/psi ("clean-1")
        #
        lines [ "psi1" ][   'L0TOS' ] = 'L0(DiMuon|Muon)Decision'
        lines [ "psi1" ][   'L0TIS' ] = 'L0(Hadron|DiMuon|Muon|Electron|Photon)Decision'
        lines [ "psi1" ][ 'Hlt1TOS' ] = 'Hlt1(DiMuon|TrackMuon).*Decision'
        lines [ "psi1" ][ 'Hlt1TIS' ] = 'Hlt1(DiMuon|SingleMuon|Track).*Decision'
        lines [ "psi1" ][ 'Hlt2TOS' ] = 'Hlt2DiMuon.*Decision'
        lines [ "psi1" ][ 'Hlt2TIS' ] = 'Hlt2(Charm|Topo|DiMuon|Single).*Decision'
        #
        ## J/psi ("clean-2")
        #
        lines [ "psi2" ][   'L0TOS' ] = 'L0(DiMuon|Muon)Decision'
        lines [ "psi2" ][   'L0TIS' ] = 'L0(Hadron|DiMuon|Muon|Electron|Photon)Decision'
        lines [ "psi2" ][ 'Hlt1TOS' ] = 'Hlt1DiMuonHighMass.*Decision'
        lines [ "psi2" ][ 'Hlt1TIS' ] = 'Hlt1(DiMuon|SingleMuon|Track).*Decision'
        lines [ "psi2" ][ 'Hlt2TOS' ] = 'Hlt2DiMuonDetached(Heavy|JPsi)Decision'
        lines [ "psi2" ][ 'Hlt2TIS' ] = 'Hlt2(Charm|Topo|DiMuon|Single).*Decision'
        #
        ## J/psi ("unbiased")
        #
        lines [ "psi3" ][   'L0TOS' ] = 'L0(DiMuon|Muon)Decision'
        lines [ "psi3" ][   'L0TIS' ] = 'L0(Hadron|DiMuon|Muon|Electron|Photon)Decision'
        lines [ "psi3" ][ 'Hlt1TOS' ] = 'Hlt1DiMuonHighMass.*Decision'
        lines [ "psi3" ][ 'Hlt1TIS' ] = 'Hlt1(DiMuon|SingleMuon|Track).*Decision'
        lines [ "psi3" ][ 'Hlt2TOS' ] = 'Hlt2DiMuonJPsiHighPT.*Decision'
        lines [ "psi3" ][ 'Hlt2TIS' ] = 'Hlt2(Charm|Topo|DiMuon|Single).*Decision'


        sc = Algo.initialize(self)
        if sc.isFailure():
            return sc

        sc = self.fill_initialize()
        if sc.isFailure():
            return sc

        sc = self.tisTos_initialize ( triggers , lines )
        if sc.isFailure () : return sc

        return SUCCESS

    # finalize & print histos
    def finalize(self):
        self.fill_finalize()
        self.tisTos_finalize ()
        return Algo.finalize(self)

    def analyse(self):
        MyB = self.select(
            'psi', '[( B+ ->  ( J/psi(1S) ->  mu+  mu-  )  K+  K-  pi+ )]CC')
        if MyB.empty():
            return self.Warning("No B+ are found!", SUCCESS)

        # Run and Event Numbers definition using ODIN
        # odin = self.get('DAQ/ODIN')
        # runNum = odin.runNumber()
        # evNum = odin.eventNumber()

        # PV
        prims = self.vselect('prims', ISPRIMARY)
        if prims.empty():
            return self.Warning("No primary vertices", SUCCESS)

        # Minimal impact parameter chi2
        mipFun = MIP(prims, self.geo())
        CHI2mipFun = MIPCHI2(prims, self.geo())
        MIPCHI2DVfun= MIPCHI2DV()

        DLLKpi = PIDK - PIDpi
        DLLKp = PIDK - PIDp

        # Constrains
        dtffun_ctau = DTF_CTAU(0, True)
        dtffun_chi2 = DTF_CHI2NDOF(True, "J/psi(1S)")
        dtffun_m = DTF_FUN(M, True, "J/psi(1S)")

        dtffun_m123 = DTF_FUN(MASS(1, 2, 3), True, "J/psi(1S)")
        dtffun_m23 = DTF_FUN(M23, True, "J/psi(1S)")
        dtffun_m34 = DTF_FUN(M34, True, "J/psi(1S)")
        dtffun_m234 = DTF_FUN(MASS(2, 3, 4), True, "J/psi(1S)")

        nt = self.nTuple("t")

        for myb in MyB:

            b, jpsi, k1, k2, pi = tuple(myb(i) for i in xrange(5))
            #b    = myb(0)
            #jpsi = myb(1)
            #k1   = myb(2)
            #k2   = myb(3)
            #pi   = myb(4)

            self.treatKine(nt, b, '_b')
            self.treatKine(nt, jpsi, '_jpsi')

            # add the information for Pid efficiency correction
            self.treatPions(nt, b)
            self.treatKaons(nt, b)
            self.treatMuons(nt, b)
            self.treatTracks(nt, b)

            # add DTF-applied information
            nt.column('DTFctau', dtffun_ctau(myb))
            nt.column('DTFchi2ndof', dtffun_chi2(myb))
            nt.column('DTFm_b', dtffun_m(myb) / GeV)

            nt.column('DTFm_jpsikk', dtffun_m123(myb) / GeV)
            nt.column('DTFm_kk', dtffun_m23(myb) / GeV)
            nt.column('DTFm_kpi', dtffun_m34(myb) / GeV)
            nt.column('DTFm_kkpi', dtffun_m234(myb) / GeV)

            nt.column('MIPCHI2DV_k1', MIPCHI2DVfun(k1))
            nt.column('MIPCHI2DV_k2', MIPCHI2DVfun(k2))
            nt.column('MIPCHI2DV_pi', MIPCHI2DVfun(pi))

            # add the information needed for TisTos
            self.tisTos ( jpsi  , nt  , 'psi_' ,
                          self.lines [ 'psi' ] , self.l0tistos , self.l1tistos , self.l2tistos )

            self.tisTos ( jpsi  , nt  , 'psi1_' ,
                          self.lines [ 'psi1' ] , self.l0tistos , self.l1tistos , self.l2tistos )

            self.tisTos ( jpsi  , nt  , 'psi2_' ,
                          self.lines [ 'psi2' ] , self.l0tistos , self.l1tistos , self.l2tistos )

            self.tisTos ( jpsi  , nt  , 'psi3_' ,
                          self.lines [ 'psi3' ] , self.l0tistos , self.l1tistos , self.l2tistos )


            nt.write()

        self.setFilterPassed(True)

        return SUCCESS

# =============================================================================
# configure the job


def configure(datafiles, catalogs=[], params={}, castor=False):
    from Configurables import DaVinci
    from Configurables import EventSelector

    from PhysConf.Filters import LoKi_Filters

    fltrs = LoKi_Filters(
        VOID_Code="""
        0 < CONTAINS ('/Event/PSIX/Phys/SelPsi3KPiForPsiX/Particles') 
        """
    )
    filters = fltrs.filters('Filters')

    daVinci = DaVinci(
        DataType=params['year'],
        PrintFreq=1000,
        EvtMax=-1,
        #
        TupleFile='B2JpsiKKpi.root',
        EventPreFilters=filters,
        #
        Lumi=True,
        #RootInTES = 'PSIX'
        #
    )

    from Configurables import Gaudi__IODataManager as IODataManager
    IODataManager().AgeLimit = 2

    # define input data
    setData(datafiles, catalogs, castor)

    #
    # Dynamic Configuration: Jump into the wonderful world of GaudiPython
    #

    # get the actual application manager (create if needed)
    gaudi = appMgr()

    # create local algorithm:
    alg = B2JpsiKKpi(
        'JpsiKKpi',
        # input particles :
        RootInTES='/Event/PSIX',
        Inputs=['Phys/SelPsi3KPiForPsiX/Particles']
    )
    mainSeq = gaudi.algorithm('GaudiSequencer/DaVinciUserSequence', True)
    mainSeq.Members += [alg.name()]

    # gaudi.addAlgorithm ( alg )

    return SUCCESS


if __name__ == '__main__':

    # make printout of the own documentations
    print '*' * 120
    print __doc__
    print ' Author  : %s ' % __author__
    print ' Version : %s ' % __version__
    print ' Date    : %s ' % __date__
    print '*' * 120

    inputdata = ['/lhcb/LHCb/Collision12/PSIX.MDST/00033675/0000/00033675_00000001_1.psix.mdst']
    params = {'year' : '2012'}

    configure(inputdata, params=params, castor=True)
    run(200)
