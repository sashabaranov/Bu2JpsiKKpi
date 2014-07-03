#!/bin/bash
PORT=8001

ROOT_VERSION=67
check() {
    which runpy > /dev/null
    if [ $? -eq 1 ] ; then
        echo "Did you run `source /afs/cern.ch/lhcb/software/releases/LBSCRIPTS/LBSCRIPTS_v7r10p1/InstallArea/scripts/LbLogin.sh`?"
        exit
    fi
}

check
hostname -f
uptime
free -m

export ENV4MU=/afs/cern.ch/work/a/albarano/public/4mu_env

source SetupProject.sh Bender v24r2
source $ENV4MU/bin/activate

# Dirty hack to make things work
export PYTHONPATH=$ENV4MU/lib64/python2.7/site-packages/:$PYTHONPATH
export PYTHONPATH=$ANALYSIS/Ostap:$PYTHONPATH

ipython notebook --ip='*' --pylab inline --no-browser --port=$PORT
