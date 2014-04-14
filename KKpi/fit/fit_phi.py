import ROOT
from AnalysisPython.PyRoUts import *
from AnalysisPython.Utils import timing
from AnalysisPython.Utils import rooSilent
from AnalysisPython.Logger import getLogger
logger = getLogger(__name__)

from cuts import h3
from fits import ds_Bu
from model import model_Phi

# ds_Bu.draw("DTFm_kpi >> %s" % h2.GetName() , "SBu_sw")
ds_Bu.project(h3, "DTFm_kk", "SBu_sw")
h3.Draw()


with rooSilent():
    logger.info('Fit Phi histogram (check the model)')
    r, f = model_Phi.fitHisto(h3)
    model_Phi.signal.mean.release()
    r, f = model_Phi.fitHisto(h3, draw=True)