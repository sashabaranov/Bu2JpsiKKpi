import ROOT
from AnalysisPython.PyRoUts import *
from AnalysisPython.Utils import timing
from AnalysisPython.Utils import rooSilent
from AnalysisPython.Logger import getLogger
logger = getLogger(__name__)


from cuts import m_Bu, nbin_Bu
from cuts import m_Kstar, nbin_Kstar
from cuts import m_Phi, nbin_Phi

from cuts import cuts_Bu, prntCuts
from cuts import h1
logger.info('Import models/PDFs')
from model import model_Bu


logger.info('Import data (for the first time it could take some time)')
from data_Bc import tBu2 as tBu


logger.info('DATA chain name is %s ' % (tBu.GetName()))

data = []

for cutval in xrange(1,9):
    data_cuts_Bu = cuts_Bu + "&& c2ip_b    < %d " % cutval

    # logger.info('Fill control B+  histogram (takes some time)')
    with timing():
        tBu.Project(h1.GetName(), 'DTFm_b', data_cuts_Bu)

    with rooSilent():
        # logger.info('Fit Bc+ & B+ histogram (check the model)')
        r, f = model_Bu.fitHisto(h1)
        model_Bu.signal.mean .release()
        model_Bu.signal.sigma.release()
        r, f = model_Bu.fitHisto(h1)
        r, f = model_Bu.fitHisto(h1, draw=True)

    from Selectors import SBT

    sel_Bu = SBT(m_Bu, cuts_Bu)


    # logger.info(
    #     'Build RooFit dataset for B+ , it could take as long as 3-5 minutes')
    with timing():
        ## logger.warning('almost skip B+')
        ## with rooSilent () : tBu.process ( sel_Bu , 100000 )
        with rooSilent():
            tBu.process(sel_Bu)

    ds_Bu = sel_Bu.dataset()
    ds_Bu.Print('v')

    with rooSilent():
        logger.info('Make unbinned fit for B+')
        #
        model_Bu.s.setMax(1.2 * len(ds_Bu))
        ru, fu = model_Bu.fitTo(ds_Bu, draw=True, nbins=nbin_Bu)
        model_Bu.signal.mean .fix(ru('mean_Bu1')[0])
        model_Bu.signal.sigma.fix(ru('sigma_Bu1')[0])
        
        ru, fu = model_Bu.fitTo(ds_Bu, draw=True, nbins=nbin_Bu)

    print "c2ip < %d  |  %4.4e +/- %4.2e  |  %4.5f" % (cutval , ru("SBu")[0].value(), ru("SBu")[0].error(), ru("SBu")[0].prec().value())
    data.append((cutval, ru("SBu")[0].prec().value()))

    # logger.info('running sPlot')
    # model_Bu.sPlot(ds_Bu)

      
    logger.info('end of the module')
