my_area = '/afs/cern.ch/user/a/albarano/cmtuser'

j = Job(
    name='B2JpsiKKpi',
    application=Bender(
        events=-1,
        version='v22r8',
        user_release_area=my_area,
        module=my_area +
        '/Bender_v22r8/Scripts/testing/B2JpsiKKpi.py'
    ),
    outputdata=['B2JpsiKKpi.root'],
    backend=Dirac(),
    splitter=SplitByFiles(filesPerJob=50, maxFiles=-1),
    inputdata=DaVinci(version='v32r2p1').readInputData('selection5-2011.py'),
)
j.splitter.ignoremissing = True
j.submit()
