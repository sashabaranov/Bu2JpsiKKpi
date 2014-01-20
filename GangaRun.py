my_area = '/afs/cern.ch/user/a/albarano/cmtuser'

j2011 = Job(
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
    splitter=SplitByFiles(filesPerJob=10, maxFiles=-1)
)

j2011.splitter.ignoremissing = True

j2011.name = 'RD-KKpi-2011'
j2011.inputdata = DaVinci(version='v32r2p1').readInputData('selection6-2011.py')
j2011.application.params = {'year': '2011'}

j2011.submit()

j2012 = Job(
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
    splitter=SplitByFiles(filesPerJob=10, maxFiles=-1)
)

j2012.name = 'RD-KKpi-2012'
j2012.inputdata = DaVinci(version='v32r2p1').readInputData('selection6-2012.py')
j2012.application.params = {'year': '2012'}


j2012.submit()
