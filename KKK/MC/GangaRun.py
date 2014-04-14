#!/usr/bin/env ganga

my_area = '/afs/cern.ch/user/a/albarano/cmtuser'


j = Job(
    application=Bender(
        events=-1,
        version='v22r8',
        user_release_area=my_area,
        module=my_area +
        '/Bender_v22r8/Scripts/Bu2JpsiKKK/MC/Bu2JpsiKKK.py'
    ),
    outputdata=['DVNtuples.root'],
    backend=Dirac(),
    splitter=SplitByFiles(filesPerJob=3, maxFiles=-1),
)


# Bu -> J/psi K K K
bu_entries = getBKInfo(12245001)

tasks = bu_entries

keys = tasks.keys()
keys.sort()

for key in keys:
    if 'new' != j.status:
        j = j.copy()

    j.name = 'MC-KKK'

    task = tasks[key]

    num_files = task[2]
    num_events = task[3]

    if num_files < 1:
        continue
    if num_events < 1:
        continue

    path = key

    
    mode = 'B+ -> J/psi K+ K+ K-'

    year = '2011'

    if 0 <= path.find('/MC/2012'):
        year = '2012'

    j.name = 'KKK/MC/' + year

    dddb = task[0]
    cond = task[1]

    if not dddb:
        print 'INVALID DDDB-tag', task
        continue

    if not cond:
        print 'INVALID CONDDB-tag', task
        continue

    params = {
        'Mode': mode,
        'Year': year,
        'DDDB': dddb,
        'SIMCOND': cond,
    }

    j.application.params = params

    query = BKQuery(path)
    j.inputdata = query.getDataset()

    print  key, task, len(j.inputdata)

    j.comment = str([path] + list(task))
    print 'COMMENT:', j.comment
    j.submit()

jobs

# =============================================================================
# The END
# =============================================================================
