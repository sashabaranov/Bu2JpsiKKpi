my_area = '/afs/cern.ch/user/a/albarano/cmtuser'


output =  [ DiracFile("*.mdst") , 
            SandboxFile("*.xml") , 
            SandboxFile("*db") , 
            SandboxFile("*.txt") , 
            SandboxFile("*.root") 
]

j2011 = Job(
    name='B2JpsiKKpi',
    application=Bender(
        events=-1,
        version='v23r3',
        user_release_area=my_area,
        module='$KKpidir/B2JpsiKKpi.py'
    ),
    outputfiles=output,
    backend=Dirac(),
    splitter=SplitByFiles(filesPerJob=1, maxFiles=-1, bulksubmit = True),
)

j2011.splitter.ignoremissing = False

j2011.name = 'RD-KKpi-11'
j2011.inputdata = DaVinci(version='v33r9').readInputData('$ANALYSIS/selections/selection7-2011.py')
j2011.application.params = {'year': '2011'}

j2011.submit()

j2012 = Job(
    name='B2JpsiKKpi',
    application=Bender(
        events=-1,
        version='v23r3',
        user_release_area=my_area,
        module='$KKpidir/B2JpsiKKpi.py'
    ),
    outputfiles=output,
    backend=Dirac(),
    splitter=SplitByFiles(filesPerJob=1, maxFiles=-1, bulksubmit = True)
)

j2012.splitter.ignoremissing = False

j2012.name = 'RD-KKpi-12'
j2012.inputdata = DaVinci(version='v33r9').readInputData('$ANALYSIS/selections/selection7-2012.py')
j2012.application.params = {'year': '2012'}


j2012.submit()
