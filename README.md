Bu → J/ψ K K π
===========

LHCb physics analysis

To setup use next commands:

```
git clone git@github.com:scr4t/Bu2JpsiKKpi.git bu_jpsikkpi_analysis
cd bu_jpsikkpi_analysis
export ANALYSIS=$(pwd)
export KKKdir=$ANALYSIS/KKK
export KKpidir=$ANALYSIS/KKpi
export Kpipidir=$ANALYSIS/Kpipi
```


(For everyday use better place `export`'s to `.bash_profile`)

To use data(on lxplus):

```
export DATADIR=/afs/cern.ch/work/a/albarano/data
ln -s $DATADIR/KKK $KKKdir/output 
ln -s $DATADIR/KKpi $KKpidir/output
ln -s $DATADIR/Kpipi $Kpipidir/output
```
