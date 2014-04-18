#!/usr/bin/env ganga
# Bu -> J/psi K K K: 12245001

my_area = '/afs/cern.ch/user/a/albarano/cmtuser'

mode = 'B+ -> J/psi K+ K+ K-'

config = {
    "KKK-MC11a" : [
        ('/MC/MC11a/Beam3500GeV-2011-MagUp-Nu2-EmNoCuts/Sim05d/Trig0x40760037Flagged/Reco12a/Stripping17NoPrescalingFlagged/12245001/ALLSTREAMS.DST', 'MC2011-20120727', 'MC2011-20120727-vc-mu100'),
        ('/MC/MC11a/Beam3500GeV-2011-MagDown-Nu2-EmNoCuts/Sim05d/Trig0x40760037Flagged/Reco12a/Stripping17NoPrescalingFlagged/12245001/ALLSTREAMS.DST', 'MC2011-20120727', 'MC2011-20120727-vc-md100'),
    ],

    "KKK-2012-Pythia6": [
        ('/MC/2012/Beam4000GeV-2012-MagDown-Nu2.5-Pythia6/Sim08b/Digi13/Trig0x409f0045/Reco14a/Stripping20NoPrescalingFlagged/12245001/ALLSTREAMS.DST', 'Sim08-20130503-1', 'Sim08-20130503-1-vc-md100'),
        ('/MC/2012/Beam4000GeV-2012-MagUp-Nu2.5-Pythia6/Sim08b/Digi13/Trig0x409f0045/Reco14a/Stripping20NoPrescalingFlagged/12245001/ALLSTREAMS.DST', 'Sim08-20130503-1', 'Sim08-20130503-1-vc-mu100'),
    ],

    "KKK-2012-Pythia8": [
        ('/MC/2012/Beam4000GeV-2012-MagDown-Nu2.5-Pythia8/Sim08b/Digi13/Trig0x409f0045/Reco14a/Stripping20NoPrescalingFlagged/12245001/ALLSTREAMS.DST', 'Sim08-20130503-1', 'Sim08-20130503-1-vc-md100'),
        ('/MC/2012/Beam4000GeV-2012-MagUp-Nu2.5-Pythia8/Sim08b/Digi13/Trig0x409f0045/Reco14a/Stripping20NoPrescalingFlagged/12245001/ALLSTREAMS.DST', 'Sim08-20130503-1', 'Sim08-20130503-1-vc-mu100'),
    ],

    "KKK-2011-Pythia6": [
        ('/MC/2011/Beam3500GeV-2011-MagDown-Nu2-Pythia6/Sim08b/Digi13/Trig0x40760037/Reco14a/Stripping20r1NoPrescalingFlagged/12245001/ALLSTREAMS.DST', 'Sim08-20130503', 'Sim08-20130503-vc-md100'),
        ('/MC/2011/Beam3500GeV-2011-MagUp-Nu2-Pythia6/Sim08b/Digi13/Trig0x40760037/Reco14a/Stripping20r1NoPrescalingFlagged/12245001/ALLSTREAMS.DST', 'Sim08-20130503', 'Sim08-20130503-vc-mu100'),
    ],

    "KKK-2011-Pythia8": [
        ('/MC/2011/Beam3500GeV-2011-MagDown-Nu2-Pythia8/Sim08b/Digi13/Trig0x40760037/Reco14a/Stripping20r1NoPrescalingFlagged/12245001/ALLSTREAMS.DST', 'Sim08-20130503', 'Sim08-20130503-vc-md100'),
        ('/MC/2011/Beam3500GeV-2011-MagUp-Nu2-Pythia8/Sim08b/Digi13/Trig0x40760037/Reco14a/Stripping20r1NoPrescalingFlagged/12245001/ALLSTREAMS.DST', 'Sim08-20130503', 'Sim08-20130503-vc-mu100'),
    ],
}

output =  [ DiracFile("*.mdst") ,
            SandboxFile("*.xml") ,
            SandboxFile("*db") ,
            SandboxFile("*.txt") ,
            SandboxFile("*.root")
]

j = Job(
    application=Bender(
        events=-1,
        version='v23r3',
        user_release_area=my_area,
        module='$KKKdir/MC/Bu2JpsiKKK.py'
    ),
    outputfiles=output,
    backend=Dirac(),
    splitter=SplitByFiles(filesPerJob=3, maxFiles=-1, bulksubmit = True)
)


for job_prefix, job_configs in config.items():
    year = "2011"
    if "2012" in job_prefix:
        year = "2012"

    for path, dddb, cond in job_configs:
        mag = "-U"
        if "Down" in path:
            mag = "-D"

        if j.status != 'new':
            j = j.copy()

        j.name = job_prefix + mag,
        j.splitter.ignoremissing = False

        j.application.params = {
            'Mode': mode,
            'Year': year,
            'DDDB': dddb,
            'SIMCOND': cond,
        }

        j.inputdata = BKQuery(path).getDataset()

        j.comment = str((path, dddb, cond))
        j.submit()
