my_area = '/afs/cern.ch/user/a/albarano/cmtuser'

j2011 = Job(
    name='RD-Kpipi-11',
    application=Bender(
        events=-1,
        version='v22r8',
        user_release_area=my_area,
        module=my_area +
        '/Bender_v22r8/Scripts/Bu2JpsiKpipi/Bu2JpsiKpipi.py'
    ),
    outputdata=['Bu2JpsiKpipi.root'],
    backend=Dirac(),
    splitter=SplitByFiles(filesPerJob=5, maxFiles=-1)
)

j2011.splitter.ignoremissing = True
j2011.inputdata = DaVinci(version='v32r2p1').readInputData('selection6-2011.py')
j2011.application.params = {'year': '2011'}

#j2011.submit()



j2012 = Job(
    name='RD-Kpipi-12',
    application=Bender(
        events=-1,
        version='v22r8',
        user_release_area=my_area,
        module=my_area +
        '/Bender_v22r8/Scripts/Bu2JpsiKpipi/Bu2JpsiKpipi.py'
    ),
    outputdata=['Bu2JpsiKpipi.root'],
    backend=Dirac(),
    splitter=SplitByFiles(filesPerJob=1, maxFiles=-1)
)

j2012.splitter.ignoremissing = True
j2012.inputdata = DaVinci(version='v32r2p1').readInputData('selection6-2012.py')
j2012.application.params = {'year': '2012'}

j2012.submit()
