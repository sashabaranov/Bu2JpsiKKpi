Bu → J/ψ K K π
===========

LHCb physics analysis

To setup use next commands:

```
git clone git@github.com:scr4t/Bu2JpsiKKpi.git bu_jpsikkpi_analysis
cd bu_jpsikkpi_analysis
git submodule init
git submodule update
export ANALYSIS=$(pwd)
export KKKdir=$ANALYSIS/KKK
export KKpidir=$ANALYSIS/KKpi
export Kpipidir=$ANALYSIS/Kpipi
export PYTHONPATH=$PYTHONPATH:$ANALYSIS/Ostap # Needed for some scripts
```

In order to use data nTuples, you should set `$WORKDIR` variable, for afs it could be:

```
export WORKDIR=/afs/cern.ch/work/a/albarano/public
```


(For everyday use better place `export`'s to `.bash_profile`)


Also, I use use library named [Ostap](https://github.com/scr4t/Ostap) for common task, and Bender as analysis framework.


