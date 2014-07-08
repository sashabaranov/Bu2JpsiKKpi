my_area = '/afs/cern.ch/user/a/albarano/cmtuser'


output =  [ DiracFile("*.mdst") ,
            SandboxFile("*.xml") ,
            SandboxFile("*db") ,
            SandboxFile("*.txt") ,
            SandboxFile("*.root")
]

j2011 = Job(
    name='RD-KKK-11',
    application=Bender(
        events=-1,
        version='v24r2',
        user_release_area=my_area,
        module='$KKKdir/Bu2JpsiKKK.py'
    ),
    outputfiles=output,
    backend=Dirac(),
    splitter=SplitByFiles(filesPerJob=5, maxFiles=-1, bulksubmit = True),
)

j2011.splitter.ignoremissing = False

j2011.inputdata = DaVinci(version='v33r9').readInputData('$ANALYSIS/selections/selection7-2011.py')
j2011.application.params = {'year': '2011'}

j2011.submit()

j2012 = Job(
    name='RD-KKK-12',
    application=Bender(
        events=-1,
        version='v24r2',
        user_release_area=my_area,
        module='$KKKdir/Bu2JpsiKKK.py'
    ),
    outputfiles=output,
    backend=Dirac(),
    splitter=SplitByFiles(filesPerJob=5, maxFiles=-1, bulksubmit = True)
)

j2012.splitter.ignoremissing = False

j2012.inputdata = DaVinci(version='v33r9').readInputData('$ANALYSIS/selections/selection7-2012.py')
j2012.application.params = {'year': '2012'}


j2012.submit()
