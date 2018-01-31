
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.Date;

import cern.colt.list.DoubleArrayList;
import hep.aida.IFunction;
import org.apache.log4j.Logger;

import gb.esac.aida.functions.SineFunction;
import gb.esac.io.AsciiDataFileFormatException;
import gb.esac.io.AsciiDataFileReader;
import gb.esac.io.AsciiDataFileWriter;
import gb.esac.timeseries.TimeSeries;
import gb.esac.timeseries.TimeSeriesMaker;
import gb.esac.timeseries.TimeSeriesOperations;
import gb.esac.tools.BasicStats;

/**

This programme calculates predictions for the evolution of the safe perigee entry and exit altitudes
using a simple sinusoidal model based on the satellite's past altitude measurements (2 years).

The altitudes are those measured at perigee entry and exit where the IREM TC3 (soft electrons) rate reads 600 counts.

@author G. Belanger, ESA, ESAC
@version 2017 July

July 3:
- Changed phase offset to year/4 in method getInitParValues, 
  because the fit on the exit data was failing using year/2.

Feb 2:
- Added safetyFactorInSigmas as an argument to the method calculateModelPredictions
- Moved prediction calculation to class ModelPredictionCalculator

**/

public class GetPerigeeHeights {

    private static Logger logger  = Logger.getLogger(GetPerigeeHeights.class);
    private static double year = 31556908.8;     
    private static DecimalFormat threeDigits = new DecimalFormat("#.000");
    private static DecimalFormat sci = new DecimalFormat("0.00E00");
    private static int nRevs;
    static int nRevsPerYear = 137;
    private static int startRev;
    private static int lastRev;
    
    private static File modelFile;
    private static double[] revNum;
    private static double[] exitAltitude;
    private static double[] entryAltitude;
    private static long[] timeInSecAtStartOfRev;
    private static long[] timeInSecAtEndOfRev;
    private static IFunction fittedFunctionExit;
    private static IFunction fittedFunctionEntry;

    // Directories for reading the input data, and writing the results
    static String data = "../data/";
    static String results = "../results/current/";

    public static void main(String[] args) throws Exception {
	
	// This method also defines lastRev and startRev=(lastRev - nRevs) with nRevs=2*nRevsPerYear;
	logger.info("Reading rad and revno files");
	readRadAndRevnoFiles();
	
	//  Handle arguments (if any are supplied)
	if ( args.length > 0 ) {
	    handleArguments(args);
	}

	//  Make complete time series of entry and exit altitudes
	TimeSeries[] tsExitAndEntry = makeTimeSeries();
	TimeSeries tsExit = tsExitAndEntry[0];
	TimeSeries tsEntry = tsExitAndEntry[1];

	// Get the segments from startRev to lastRev
	logger.info("Extracting the segments from stratRev to lastRev");
	int startIndex = startRev-1;
	int lastIndex = nRevs-2;
	TimeSeries segmentExit = TimeSeriesOperations.getSegment(tsExit, startIndex, lastIndex);
	TimeSeries segmentEntry = TimeSeriesOperations.getSegment(tsEntry, startIndex, lastIndex);
	segmentExit.writeCountsAsQDP(results+"ts_exit_"+startRev+"-"+nRevs+".qdp");
	segmentEntry.writeCountsAsQDP(results+"ts_entry_"+startRev+"-"+nRevs+".qdp");

	// Get model
	logger.info("Contructing sinusoidal models");
	double[] binCentres = segmentExit.getBinCentres();
	double[] modelExit = new double[binCentres.length];
	double[] modelEntry = new double[binCentres.length];
	boolean thisWorks = false;
	if ( thisWorks ) {
	    if ( args.length == 3 ) {
		//  Read model from QDP fit result written to .mod file
		modelExit = getModelFromFile(binCentres, modelFile);
		modelEntry = getModelFromFile(binCentres, modelFile);
		//  Perform the fit to get:
		//fittedFunctionExit =
		//fittedFunctinEntry =
	    }
	}
	else {
	    //  Get model from fitting the data
	    String label = "EXIT";	    
	    double[] initParValues = getInitParValues(segmentExit);
	    double[][] parBounds = getParBounds(segmentExit);
	    fittedFunctionExit = SinusoidModelProvider.getModel(binCentres, segmentExit.getBinHeights(), initParValues, parBounds, label);
	    // ENTRY
	    label = "ENTRY";	    
	    initParValues = getInitParValues(segmentEntry);
	    parBounds = getParBounds(segmentEntry);
	    fittedFunctionEntry = SinusoidModelProvider.getModel(binCentres, segmentEntry.getBinHeights(), initParValues, parBounds, label);
	    // Get models from fitted functions
	    for ( int i=0; i < binCentres.length; i++ ) {
		modelExit[i] = fittedFunctionExit.value(new double[] {binCentres[i]});
		modelEntry[i] = fittedFunctionEntry.value(new double[] {binCentres[i]});
	    }
	}

	//  Calculate the model prediction
	double safetyFactorInSigmas = 2;
	ModelPredictionCalculator.calculateModelPredictions(startRev, nRevs, segmentExit, segmentEntry, fittedFunctionExit, fittedFunctionEntry, modelExit, modelEntry, safetyFactorInSigmas);
    }
    // END main

    
    private static double[] getInitParValues(TimeSeries ts) {
	// parameters: period, phase, amplitude, yOffset
	double[] binHeights = ts.getBinHeights();
	double meanHeight = ts.meanBinHeight();
	double meanDev = Math.sqrt(ts.varianceInBinHeights());
	double amplitude = meanDev; 
	return new double[] {year, year/4., amplitude, meanHeight}; 
    }
    
    private static double[][] getParBounds(TimeSeries ts) {
	double quarter = year/4;
	double third = year/3;
	double[] periodBounds = new double[] {year-third, year+third};
	double[] phaseBounds = new double[] {0, year};
	double[] amplBounds = new double[] {3, 10*Math.sqrt(ts.varianceInBinHeights())};
	double[] yOffsetBounds = new double[] {ts.minBinHeight(), ts.maxBinHeight()};
	return new double[][]{periodBounds, phaseBounds, amplBounds, yOffsetBounds};
    }
    
    
    private static TimeSeries[] makeTimeSeries() throws Exception {
	logger.info("Making time series of exit and entry heights");
	long tZero = timeInSecAtStartOfRev[0];
	double[] durations = new double[nRevs];
	double[] binEdges = new double[2*nRevs];
	DoubleArrayList binEdgesList_minBins = new DoubleArrayList();
	DoubleArrayList binEdgesList_allBins = new DoubleArrayList();
	double[] binHeightsExit = new double[nRevs];
	DoubleArrayList binHeightsExitList_minBins = new DoubleArrayList();
	DoubleArrayList binHeightsExitList_allBins = new DoubleArrayList();
	double[] binHeightsEntry = new double[nRevs];
	DoubleArrayList binHeightsEntryList_minBins = new DoubleArrayList();
	DoubleArrayList binHeightsEntryList_allBins = new DoubleArrayList();
	for ( int i=0; i < nRevs; i++ ) {
	    timeInSecAtStartOfRev[i] -= tZero;
	    timeInSecAtEndOfRev[i] -= tZero;
	    durations[i] = timeInSecAtEndOfRev[i] - timeInSecAtStartOfRev[i];
	    binEdges[2*i] = timeInSecAtStartOfRev[i];
	    binEdges[2*i+1] = timeInSecAtEndOfRev[i];
	    binEdgesList_allBins.add(binEdges[2*i]);
	    binEdgesList_allBins.add(binEdges[2*i+1]);
	    binHeightsExit[i] = exitAltitude[i];
	    binHeightsEntry[i] = entryAltitude[i];
	    binHeightsExitList_allBins.add(binHeightsExit[i]);
	    binHeightsEntryList_allBins.add(binHeightsEntry[i]);
	    // Exclude bins for which the intensity is NaN
	    // This introduces gaps in sampling but eliminates NaNs from the rates
	    if ( !Double.isNaN(binHeightsExit[i]) && !Double.isNaN(binHeightsEntry[i]) ) {
		binHeightsExitList_minBins.add(binHeightsExit[i]);
		binHeightsEntryList_minBins.add(binHeightsEntry[i]);
		binEdgesList_minBins.add(binEdges[2*i]);
		binEdgesList_minBins.add(binEdges[2*i+1]);
	    }
	}
	binHeightsExitList_allBins.trimToSize();
	binHeightsEntryList_allBins.trimToSize();
	binEdgesList_allBins.trimToSize();
	binHeightsExitList_minBins.trimToSize();
	binHeightsEntryList_minBins.trimToSize();
	binEdgesList_minBins.trimToSize();

	////  Write Exit TS
	TimeSeries tsExit_minBins = TimeSeriesMaker.makeTimeSeries(binEdgesList_minBins.elements(), binHeightsExitList_minBins.elements());
	tsExit_minBins.writeCountsAsQDP(results+"ts_exit_noNaNs_gaps_all.qdp");
	TimeSeries tsExit_allBins = TimeSeriesMaker.makeTimeSeries(binEdgesList_allBins.elements(), binHeightsExitList_allBins.elements());
	tsExit_allBins.writeCountsAsQDP(results+"ts_exit_withNaNs_noGaps_all.qdp");

	////  Write Entry TS
	
	TimeSeries tsEntry_minBins = TimeSeriesMaker.makeTimeSeries(binEdgesList_minBins.elements(), binHeightsEntryList_minBins.elements());
	tsEntry_minBins.writeCountsAsQDP(results+"ts_entry_noNaNs_gaps_all.qdp");
	TimeSeries tsEntry_allBins = TimeSeriesMaker.makeTimeSeries(binEdgesList_allBins.elements(), binHeightsEntryList_allBins.elements());
	tsEntry_allBins.writeCountsAsQDP(results+"ts_entry_withNaNs_noGaps_all.qdp");

	return new TimeSeries[]{tsExit_allBins, tsEntry_allBins};
    }
    

    private static void readOrbitFileAndMakePlots() throws Exception {
 	AsciiDataFileReader orbut = new AsciiDataFileReader(data+"orbit.dat");
	String[] perDates = orbut.getStrCol(1);
	double[] perigee = orbut.getDblCol(2);
	String[] apoDates = orbut.getStrCol(3);
	double[] apogee = orbut.getDblCol(4);
	double[] eccen = orbut.getDblCol(5);
	double[] inclAngle = orbut.getDblCol(6);
	double[] argPer = orbut.getDblCol(7);
	SimpleDateFormat perFormat = new SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ss");
	double[] timeInSecAtPerigee = new double[nRevs];
	double[] timeInSecAtApogee = new double[nRevs];
	double[] passageTimes = new double[nRevs];
	for ( int i=0; i < nRevs; i++ ) {
	    Date per = perFormat.parse(perDates[i]);
	    Calendar cal = Calendar.getInstance();
	    cal.setTime(per);
	    String year = String.valueOf(cal.get(Calendar.YEAR));
	    Date apo = perFormat.parse(year+"-"+apoDates[i]);
	    long diff = apo.getTime() - per.getTime();
	    if ( diff < 0 ) {
		year = String.valueOf(cal.get(Calendar.YEAR)+1);
		apo = perFormat.parse(year+"-"+apoDates[i]);
	    }
	    timeInSecAtPerigee[i] = per.getTime()/1000;
	    timeInSecAtApogee[i] = apo.getTime()/1000;
	    passageTimes[i] = timeInSecAtApogee[i] - timeInSecAtPerigee[i];
	}
	double  zero = timeInSecAtPerigee[0];
	for ( int i=0; i < nRevs; i++ ) {
	    timeInSecAtPerigee[i] -= zero;
	    timeInSecAtApogee[i] -= zero;
	}
	AsciiDataFileWriter out;
	String[] header = new String[] {
	    "DEV /XS",
	    "LAB F",
	    "TIME OFF",
	    "VIEW 0.2 0.1 0.8 0.9",
	    "CS 1.5",
	    "LW 4",
	};
	out = new AsciiDataFileWriter(results+"perigee.qdp");
	out.writeData(header, timeInSecAtPerigee, perigee);
	out = new AsciiDataFileWriter(results+"apogee.qdp");
	out.writeData(header, timeInSecAtApogee, apogee);
	out = new AsciiDataFileWriter(results+"inclAngle.qdp");
	out.writeData(header, timeInSecAtPerigee, inclAngle);
	out = new AsciiDataFileWriter(results+"eccen.qdp");
 	out.writeData(header, timeInSecAtPerigee, eccen);
 	out = new AsciiDataFileWriter(results+"argPerigee.qdp");
 	out.writeData(header, timeInSecAtPerigee, argPer);
	//  Correlations
	double[] correlCoef = BasicStats.getCorrelationCoefficient(exitAltitude, inclAngle);
	logger.info("Correlation coefficients between:");
	logger.info("   inclination angle and exit altitude = "+threeDigits.format(correlCoef[0])+ "\t"+ threeDigits.format(correlCoef[1]));
	out = new AsciiDataFileWriter(results+"exitAltitudeVsInclAngle.qdp");
	out.writeData(header, exitAltitude, inclAngle);
	correlCoef = BasicStats.getCorrelationCoefficient(argPer, exitAltitude);
	logger.info("   argument of perigee and exit altitude = "+threeDigits.format(correlCoef[0])+ "\t"+ threeDigits.format(correlCoef[1]));
	out = new AsciiDataFileWriter(results+"exitAltitudeVsArgPerigree.qdp");
	out.writeData(header, argPer, exitAltitude);
    }

    
    //  This hasn't been used in a long time and will probably not work
    private static double[] getModelFromFile(double[] binCentres, File modelFile) throws AsciiDataFileFormatException, IOException  {
	//  Combined model cons sin sin sin
	AsciiDataFileReader mod = new AsciiDataFileReader(modelFile);
	double[] parValues = mod.getDblCol(0);
	String[] parNames = mod.getStrCol(mod.getNDataCols()-1);
	//  Parameters for each sine are: 
	//  CO PE PH SN
	//  constant, period, phase, norm/amplitude

	//  Define the 4 par values for each of the 3 sine components
	int nSines = 3;
	IFunction[] sines = new IFunction[nSines];
	for ( int i=0; i < nSines; i++ ) {
	    sines[i] = new SineFunction("sine");
	    sines[i].setParameter("yOffset", parValues[4*i]);
	    sines[i].setParameter("period", parValues[4*i+1]);
	    sines[i].setParameter("xOffset", parValues[4*i+2]);
	    sines[i].setParameter("amplitude", parValues[4*i+3]);
	}
	//  Calculate combined model for each rev
	int n = binCentres.length;
	double[] model = new double[n];
	for ( int i=0; i < n; i++ ) {
	    for ( int j=0; j < nSines; j++ ) {
		model[i] += sines[j].value(new double[] {binCentres[i]});
	    }	    
	}
	return model;
    }


    private static void readRadAndRevnoFiles() throws Exception {
	// read radiation data file downloaded from IDSC webpage
	AsciiDataFileReader in = new AsciiDataFileReader(data+"rad.dat");
	revNum = in.getDblCol(0);
	nRevs = revNum.length-1;
	lastRev = (int)revNum[nRevs-1];
	int nRevsToModel = 2*nRevsPerYear;
	startRev  = lastRev - nRevsToModel;
	exitAltitude = in.getDblCol(6);
	entryAltitude = in.getDblCol(7);
	// read revno file from MOC to get correct timeline
 	AsciiDataFileReader revno = new AsciiDataFileReader(data+"revno");
 	String[] dates = revno.getStrCol(3);
	SimpleDateFormat revnoFormat = new SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ss'Z'");
	timeInSecAtStartOfRev = new long[nRevs];
	timeInSecAtEndOfRev = new long[nRevs];
	for ( int i=0; i < nRevs; i++ ) {
	    Date start = revnoFormat.parse(dates[i]);
	    Date end = revnoFormat.parse(dates[i+1]);
	    timeInSecAtStartOfRev[i] = start.getTime()/1000;
	    timeInSecAtEndOfRev[i] = end.getTime()/1000;
	}
    }

    private static void handleArguments(String[] args) throws Exception {
	// If a third argument is supplied it is interpreted as a model file
	String modelFilename = "modelUpTo1053/modelUpTo1053.mod";
	modelFilename = "modelUpTo1095/modelUpTo1095-minusFirstFourYears.mod";
	modelFilename = "modelUpTo1145/modelUpTo1145-minusFirstFourYears.mod"; // last is default model
	if ( args.length == 3 ) {
	    startRev = (Integer.valueOf(args[0])).intValue();
	    if ( args[1].equals("yes") ) {
		readOrbitFileAndMakePlots();
	    }
	    else {
		logger.info("Orbit file will not be read");
	    }
	    modelFile = new File(args[2]);
	    if ( !modelFile.exists() ) {
		logger.error("Usage: java GetPerigeeHeights (startRev) (readOrbitFile (yes|no)) (modelFilename.mod)");
		logger.error("  File not found: "+modelFile.getName());
		System.exit(-1);
	    }
	}
	else if ( args.length == 2 ) {
	    startRev = (Integer.valueOf(args[0])).intValue();
	    if ( args[1].equals("yes") ) {
		readOrbitFileAndMakePlots();
	    }
	    else {
		logger.info("Orbit file will not be read");
	    }
	}
	else if ( args.length == 1 ) {
	    startRev = (Integer.valueOf(args[0])).intValue();
	    if ( startRev >= lastRev ) {
		logger.error("Usage: java GetPerigeeHeights (startRev) (readOrbitFile (yes|no)) (modelFilename.mod)");
		logger.error("  startRev ("+startRev+") must be smaller than lastRev ("+lastRev+")");
		System.exit(-1);
	    }
	}
	else {
	    logger.error("Usage: java GetPerigeeHeights (startRev) (readOrbitFile (yes|no)) (modelFilename.mod)");	    
	}
    }

}
