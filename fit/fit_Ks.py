import ROOT
from AnalysisPython.PyRoUts import *
from AnalysisPython.Utils import timing
from AnalysisPython.Utils import rooSilent
from AnalysisPython.Logger import getLogger
logger = getLogger(__name__)

from cuts import h2
from fits import ds_Bu
from model import model_Kstar

# ds_Bu.draw("DTFm_kpi >> %s" % h2.GetName() , "SBu_sw")
ds_Bu.project(h2, "DTFm_kpi", "SBu_sw")
h2.Draw()


with rooSilent():
    logger.info('Fit Kstar histogram (check the model)')
    with timing():
        r, f = model_Kstar.fitHisto(h2)
    r, f = model_Kstar.fitHisto(h2, draw=True)
