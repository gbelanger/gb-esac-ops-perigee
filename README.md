How to get an updated prediction for the evolution of the perigee height at entry and exit.


DEFAULT:

cd src/
./runall

(Uses orbut_ap.dat that on disk, but this file is not used in the modelling or prediction):


STEP by STEP:

OPTIONAL 1) Download orbut_ap.txt (contains the orbital parameters of Integral)

  The data comes from the file called orbut_ap.txt sent around by MOC with the bubble plots.
  - Download by hand from the most recent email with subject "New INTEGRAL Mission Planning Constraint Plots: ..."

2) Execute getdatafiles : downloads and copies other required files

  - rad (wget directly from ISDC) : Contains the perigeen entry and exit altitudes
  - revno : Contains the start and end time of each revolution; copied from ops system.

  - Note: this script will also modify orbut_ap.txt by running orbut2orbit (in $IDX_DIR defined in getdatafiles)
  	  to get rid of the headers and black lines in the file.

  cd src/
  ./getdatafiles

3) Compile the java classes (just to be safe)

javac *.java

4) Run

java GetPerigeeHeights

  Optional arguments: Can run with 1 (startRev) or 2 (startRev yes|no) or 3 (startRev yes|no model)
  1) startRev : default is (lastRev-274), which is 2 years
  2) yes or no : default is 'no' - specify that the orbit.dat file should be read and plots made from it
  3) QDP model file (default is "modelUpTo1145-minusFirstFourYears.mod")

5) Combined the output files

./makeplots
