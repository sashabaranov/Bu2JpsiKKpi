#!/usr/bin/env ganga

my_area = '/afs/cern.ch/user/a/albarano/cmtuser'

mode = 'B+ -> J/psi K+ K+ K-'



config = {
    "Kpipi-MC-2012-Pythia6": [
        ('/MC/2012/Beam4000GeV-2012-MagDown-Nu2.5-Pythia6/Sim08b/Digi13/Trig0x409f0045/Reco14a/Stripping20NoPrescalingFlagged/12245021/ALLSTREAMS.DST', 'Sim08-20130503-1', 'Sim08-20130503-1-vc-md100'),
        ('/MC/2012/Beam4000GeV-2012-MagUp-Nu2.5-Pythia6/Sim08b/Digi13/Trig0x409f0045/Reco14a/Stripping20NoPrescalingFlagged/12245021/ALLSTREAMS.DST', 'Sim08-20130503-1', 'Sim08-20130503-1-vc-mu100'),

    ],

    "Kpipi-MC-2012-Pythia8": [
        ('/MC/2012/Beam4000GeV-2012-MagDown-Nu2.5-Pythia8/Sim08b/Digi13/Trig0x409f0045/Reco14a/Stripping20NoPrescalingFlagged/12245021/ALLSTREAMS.DST', 'Sim08-20130503-1', 'Sim08-20130503-1-vc-md100'),
        ('/MC/2012/Beam4000GeV-2012-MagUp-Nu2.5-Pythia8/Sim08b/Digi13/Trig0x409f0045/Reco14a/Stripping20NoPrescalingFlagged/12245021/ALLSTREAMS.DST', 'Sim08-20130503-1', 'Sim08-20130503-1-vc-mu100'),
    ],

    "Kpipi-MC-2011-Pythia6": [
        ('/MC/2011/Beam3500GeV-2011-MagDown-Nu2-Pythia6/Sim08b/Digi13/Trig0x40760037/Reco14a/Stripping20r1NoPrescalingFlagged/12245021/ALLSTREAMS.DST', 'Sim08-20130503', 'Sim08-20130503-vc-md100'),
        ('/MC/2011/Beam3500GeV-2011-MagUp-Nu2-Pythia6/Sim08b/Digi13/Trig0x40760037/Reco14a/Stripping20r1NoPrescalingFlagged/12245021/ALLSTREAMS.DST', 'Sim08-20130503', 'Sim08-20130503-vc-mu100'),

    ],

    "Kpipi-MC-2011-Pythia8": [
        ('/MC/2011/Beam3500GeV-2011-MagDown-Nu2-Pythia8/Sim08b/Digi13/Trig0x40760037/Reco14a/Stripping20r1NoPrescalingFlagged/12245021/ALLSTREAMS.DST', 'Sim08-20130503', 'Sim08-20130503-vc-md100'),
        ('/MC/2011/Beam3500GeV-2011-MagUp-Nu2-Pythia8/Sim08b/Digi13/Trig0x40760037/Reco14a/Stripping20r1NoPrescalingFlagged/12245021/ALLSTREAMS.DST', 'Sim08-20130503', 'Sim08-20130503-vc-mu100'),
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
        version='v24r2',
        user_release_area=my_area,
        module='$Kpipidir/MC/Bu2JpsiKpipi.py'
    ),
    outputfiles=output,
    backend=Dirac(),
    splitter=SplitByFiles(filesPerJob=1, maxFiles=-1, bulksubmit = True)
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

        j.name = job_prefix + mag
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
