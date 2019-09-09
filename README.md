Update prediction for altitude at perigee
==========================================


Default:
--------
```
cd src/
./runall
```

(Uses orbut_ap.dat that on disk, but this file is not used in the modelling or prediction)


Step by Step:
-------------

1. (Optional step) Download orbut_ap.txt (contains the orbital parameters of Integral)

  * The data comes from the file called orbut_ap.txt sent around by MOC with the bubble plots.
  * Download by hand from the most recent email with subject "New INTEGRAL Mission Planning Constraint Plots: ..."

2. Execute getdatafiles : downloads and copies other required files

  * rad (wget directly from ISDC) : Contains the perigeen entry and exit altitudes
  * revno : Contains the start and end time of each revolution; copied from ops system.
  * Note: this script will also modify orbut_ap.txt by running orbut2orbit (in $IDX_DIR defined in getdatafiles)
  	  to get rid of the headers and black lines in the file.
```
  cd src/
  ./getdatafiles
```

3. Compile the java classes (just to be safe)

```
javac *.java
```

4. Run

```
java GetPerigeeHeights
```

5. Make the plots

```
./makeplots
```

Warning:
-------
For list of required packages and classes see [DEPENDENCIES](https://github.com/gbelanger/gb-esac-ops-perigee/blob/master/DEPENDENCIES.md).
