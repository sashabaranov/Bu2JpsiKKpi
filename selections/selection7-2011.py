#-- GAUDI jobOptions generated on Fri Apr 18 12:40:49 2014
#-- Contains event types : 
#--   90000000 - 58 files - 7885234 events - 256.27 GBytes


#--  Extra information about the data processing phases:


#--  Processing Pass Step-126370 

#--  StepId : 126370 
#--  StepName : Merging for WGBandQSelection 
#--  ApplicationName : DaVinci 
#--  ApplicationVersion : v34r0 
#--  OptionFiles : $APPCONFIGOPTS/Merging/DV-Stripping-Merging.py 
#--  DDDB : None 
#--  CONDDB : None 
#--  ExtraPackages : AppConfig.v3r193 
#--  Visible : N 

from Gaudi.Configuration import * 
from GaudiConf import IOHelper
IOHelper('ROOT').inputFiles(['LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035294/0000/00035294_00000064_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035294/0000/00035294_00000065_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035294/0000/00035294_00000066_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035294/0000/00035294_00000068_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035294/0000/00035294_00000071_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035294/0000/00035294_00000075_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035294/0000/00035294_00000077_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035294/0000/00035294_00000078_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035294/0000/00035294_00000079_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035294/0000/00035294_00000080_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035294/0000/00035294_00000081_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035294/0000/00035294_00000082_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035294/0000/00035294_00000083_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035294/0000/00035294_00000084_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035294/0000/00035294_00000085_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035294/0000/00035294_00000086_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035294/0000/00035294_00000087_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035294/0000/00035294_00000088_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035294/0000/00035294_00000090_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035294/0000/00035294_00000091_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035294/0000/00035294_00000092_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035294/0000/00035294_00000093_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035294/0000/00035294_00000094_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035294/0000/00035294_00000098_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035298/0000/00035298_00000087_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035298/0000/00035298_00000088_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035298/0000/00035298_00000089_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035298/0000/00035298_00000091_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035298/0000/00035298_00000093_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035298/0000/00035298_00000094_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035298/0000/00035298_00000096_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035298/0000/00035298_00000097_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035298/0000/00035298_00000100_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035298/0000/00035298_00000101_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035298/0000/00035298_00000102_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035298/0000/00035298_00000110_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035298/0000/00035298_00000111_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035298/0000/00035298_00000112_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035298/0000/00035298_00000113_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035298/0000/00035298_00000115_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035298/0000/00035298_00000116_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035298/0000/00035298_00000117_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035298/0000/00035298_00000118_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035298/0000/00035298_00000119_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035298/0000/00035298_00000120_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035298/0000/00035298_00000121_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035298/0000/00035298_00000122_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035298/0000/00035298_00000123_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035298/0000/00035298_00000128_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035298/0000/00035298_00000130_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035298/0000/00035298_00000131_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035298/0000/00035298_00000132_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035298/0000/00035298_00000133_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035298/0000/00035298_00000134_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035298/0000/00035298_00000135_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035298/0000/00035298_00000136_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035298/0000/00035298_00000137_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision11/PSIX.MDST/00035298/0000/00035298_00000140_1.psix.mdst'
], clear=True)
