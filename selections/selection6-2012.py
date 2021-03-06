#-- GAUDI jobOptions generated on Mon Jan 20 16:05:52 2014
#-- Contains event types : 
#--   90000000 - 77 files - 11883916 events - 343.12 GBytes


#--  Extra information about the data processing phases:


#--  Processing Pass Step-125603 

#--  StepId : 125603 
#--  StepName : Merging for WGBandQSelection 
#--  ApplicationName : DaVinci 
#--  ApplicationVersion : v33r6p1 
#--  OptionFiles : $APPCONFIGOPTS/Merging/DV-Stripping-Merging.py 
#--  DDDB : None 
#--  CONDDB : None 
#--  ExtraPackages : AppConfig.v3r177 
#--  Visible : N 

from Gaudi.Configuration import * 
from GaudiConf import IOHelper
IOHelper('ROOT').inputFiles(['LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033675/0000/00033675_00000001_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033675/0000/00033675_00000002_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033675/0000/00033675_00000003_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033675/0000/00033675_00000004_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033675/0000/00033675_00000005_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033675/0000/00033675_00000006_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033675/0000/00033675_00000007_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033675/0000/00033675_00000008_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033675/0000/00033675_00000009_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033675/0000/00033675_00000010_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033675/0000/00033675_00000011_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033675/0000/00033675_00000012_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033675/0000/00033675_00000013_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033675/0000/00033675_00000014_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033675/0000/00033675_00000015_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033675/0000/00033675_00000016_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033675/0000/00033675_00000017_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033675/0000/00033675_00000018_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033675/0000/00033675_00000019_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033675/0000/00033675_00000020_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033675/0000/00033675_00000021_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033675/0000/00033675_00000022_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033675/0000/00033675_00000023_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033675/0000/00033675_00000024_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033675/0000/00033675_00000025_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033675/0000/00033675_00000026_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033675/0000/00033675_00000027_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033675/0000/00033675_00000028_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033675/0000/00033675_00000029_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033675/0000/00033675_00000030_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033675/0000/00033675_00000031_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033675/0000/00033675_00000032_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033675/0000/00033675_00000033_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033675/0000/00033675_00000034_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033675/0000/00033675_00000035_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033675/0000/00033675_00000036_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033675/0000/00033675_00000037_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033675/0000/00033675_00000038_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033675/0000/00033675_00000039_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033679/0000/00033679_00000001_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033679/0000/00033679_00000002_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033679/0000/00033679_00000003_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033679/0000/00033679_00000004_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033679/0000/00033679_00000005_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033679/0000/00033679_00000006_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033679/0000/00033679_00000007_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033679/0000/00033679_00000008_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033679/0000/00033679_00000009_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033679/0000/00033679_00000010_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033679/0000/00033679_00000011_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033679/0000/00033679_00000012_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033679/0000/00033679_00000013_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033679/0000/00033679_00000014_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033679/0000/00033679_00000015_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033679/0000/00033679_00000016_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033679/0000/00033679_00000017_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033679/0000/00033679_00000018_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033679/0000/00033679_00000019_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033679/0000/00033679_00000020_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033679/0000/00033679_00000021_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033679/0000/00033679_00000022_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033679/0000/00033679_00000023_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033679/0000/00033679_00000024_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033679/0000/00033679_00000025_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033679/0000/00033679_00000026_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033679/0000/00033679_00000027_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033679/0000/00033679_00000028_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033679/0000/00033679_00000029_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033679/0000/00033679_00000030_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033679/0000/00033679_00000031_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033679/0000/00033679_00000032_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033679/0000/00033679_00000033_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033679/0000/00033679_00000034_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033679/0000/00033679_00000035_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033679/0000/00033679_00000036_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033679/0000/00033679_00000037_1.psix.mdst',
'LFN:/lhcb/LHCb/Collision12/PSIX.MDST/00033679/0000/00033679_00000038_1.psix.mdst'
], clear=True)
