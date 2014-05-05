import ROOT
from AnalysisPython.PyRoot import *
from AnalysisPython.PyRoUts import *
from AnalysisPython.Utils import timing
from AnalysisPython.Utils import rooSilent
from AnalysisPython.Logger import getLogger
logger = getLogger(__name__)


from Variables import m_Bu, nbin_Bu
from cuts import cuts_Bu, prntCuts
from cuts import h1


from model import model_Bu


from data import selection7
tBu = selection7.data


logger.info('DATA chain name is %s ' % (tBu.GetName()))


for i in prntCuts(cuts_Bu, "  CUTS B+  "):
    logger.info(i)


# logger.info('Fill control B+ histogram (takes some time)')
# with timing():
#     tBu.Project(h1.GetName(), 'DTFm_b', cuts_Bu)


from PyPAW.Selectors import SelectorWithVars
from Variables import selector_variables

sel_Bu = SelectorWithVars(
    variables=selector_variables,
    selection=cuts_Bu
)

logger.info('Build RooFit dataset for B+ , it could take as long as 3-5 minutes')

tBu.process(sel_Bu)


ds_Bu = sel_Bu.dataset()
ds_Bu.Print('v')


logger.info('Make unbinned fit for B+')
#
model_Bu.s.setMax(1.2 * len(ds_Bu))
ru, fu = model_Bu.fitTo(ds_Bu, draw=True, nbins=nbin_Bu)

model_Bu.signal.sigma.release()
model_Bu.signal.mean.release()

ru, fu = model_Bu.fitTo(ds_Bu, draw=True, nbins=nbin_Bu)



print ru
print 'FIT results for B+  ', ru(model_Bu.s_name)[0]
print 'FIT precision:', ru("SBu")[0].prec()

# print 'FIT#2 results for B+  ', ru(model_Bu.s_name)[0]
# print 'FIT#2 precision:', ru("SBu")[0].prec()


logger.info('running sPlot')
model_Bu.sPlot(ds_Bu)


def make_hist(name, variable, cuts):
    h = ROOT.TH1F(name, '', 30, 5.16, 5.45)
    h.Sumw2()

    ds_Bu.project(h, variable, cuts)

    for j in xrange(0, h.GetNbinsX()):
        if h.GetBinContent(j) < 0:
            h.SetBinContent(j, 0)

    return h



hists = [
    ("pi1", "mass_pi1ask", "SBu_sw"),
    ("pi1_cuts", "mass_pi1ask", "SBu_sw && ann_pion_K > 0.1"),
]

logger.info('Writing histos')
db = shelve.open('$KKpidir/fit/histos.shelve')

db['Kpipi'] = {
    'RD': {param[0]: make_hist(*param) for param in hists}
}

db.close()

