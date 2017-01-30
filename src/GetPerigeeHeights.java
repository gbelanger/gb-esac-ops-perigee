
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.Date;

import cern.colt.list.DoubleArrayList;
import hep.aida.IAnalysisFactory;
import hep.aida.IFitFactory;
import hep.aida.IFitResult;
import hep.aida.IFitter;
import hep.aida.IFunction;
import hep.aida.IHistogram1D;
import org.apache.log4j.Logger;

import gb.esac.aida.functions.SineFunction;
import gb.esac.binner.Binner;
import gb.esac.io.AsciiDataFileFormatException;
import gb.esac.io.AsciiDataFileReader;
import gb.esac.io.AsciiDataFileWriter;
import gb.esac.periodogram.FFTPeriodogram;
import gb.esac.periodogram.ModifiedRayleighPeriodogram;
import gb.esac.periodogram.PeriodogramMaker;
import gb.esac.timeseries.TimeSeries;
import gb.esac.timeseries.TimeSeriesMaker;
import gb.esac.timeseries.TimeSeriesOperations;
import gb.esac.tools.BasicStats;

/**

This programme calculates predictions for the evolution of the safe perigee entry and exit altitudes
using a simple sinusoidal model based on the satellite's past altitude measurements (2 years).

The altitudes are those measured at perigee entry and exit where the radiation reaches a certain
level (of ...) in the radiation monitor.

@author G. Belanger, ESA, ESAC
@version 2017 January

**/

public class GetPerigeeHeights {

    private static Logger logger  = Logger.getLogger(GetPerigeeHeights.class);
    private static double year = 31556908.8;     
    private static DecimalFormat threeDigits = new DecimalFormat("#.000");
    private static DecimalFormat sci = new DecimalFormat("0.00E00");
    private static int nRevs;
    private static int nRevsPerYear = 137;
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
    private static String data = "../data/";
    private static String results = "../results/current/";

    public static void main(String[] args) throws Exception {
	
	// This method also defines lastRev and startRev=(lastRev - nRevs) with nRevs=2*nRevsPerYear;
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
	int startIndex = startRev-1;
	int lastIndex = nRevs-2;
	TimeSeries segmentExit = TimeSeriesOperations.getSegment(tsExit, startIndex, lastIndex);
	TimeSeries segmentEntry = TimeSeriesOperations.getSegment(tsEntry, startIndex, lastIndex);
	segmentExit.writeCountsAsQDP(results+"ts_exit_"+startRev+"-"+nRevs+".qdp");
	segmentEntry.writeCountsAsQDP(results+"ts_entry_"+startRev+"-"+nRevs+".qdp");

	// Get model
	double[] modelExit;
	double[] modelEntry;
	double[] binCentres = segmentExit.getBinCentres();
	if ( args.length == 3 ) {
	    //  Read model from QDP fit result written to .mod file
	    modelExit = getModelFromFile(binCentres, modelFile);
	    modelEntry = getModelFromFile(binCentres, modelFile);
	}
	else {
	    //  Get model from fitting the data
	    String label = "EXIT";	    
	    double[] initParValues = getInitParValues(segmentExit);
	    double[][] parBounds = getParBounds(segmentExit);
	    fittedFunctionExit = SinusoidModelProvider.getModel(binCentres, segmentExit.getBinHeights(), initParValues, parBounds, label);
	    modelExit = new double[binCentres.length];
	    for ( int i=0; i < binCentres.length; i++ ) {
		modelExit[i] = fittedFunctionExit.value(new double[] {binCentres[i]});
	    }
	    // ENTRY
	    label = "ENTRY";	    
	    initParValues = getInitParValues(segmentEntry);
	    parBounds = getParBounds(segmentEntry);
	    fittedFunctionEntry = SinusoidModelProvider.getModel(binCentres, segmentEntry.getBinHeights(), initParValues, parBounds, label);
	    modelEntry = new double[binCentres.length];
	    for ( int i=0; i < binCentres.length; i++ ) {
		modelEntry[i] = fittedFunctionEntry.value(new double[] {binCentres[i]});
	    }
	}
	
	//  Calculate the model prediction
	calculateModelPredictions(segmentExit, modelExit, segmentEntry, modelEntry);
    }
    // END main

    private static double[] getInitParValues(TimeSeries ts) {
	// parameters: period, phase, amplitude, yOffset
	double[] binHeights = ts.getBinHeights();
	double meanHeight = ts.meanBinHeight();
	double meanDev = Math.sqrt(ts.varianceInBinHeights());
	double amplitude = meanDev; 
	return new double[] {year, year/2, amplitude, meanHeight}; 
    }
    
    private static double[][] getParBounds(TimeSeries ts) {
	double quarter = year/4;
	double[] periodBounds = new double[] {year-quarter, year+quarter};
	double[] phaseBounds = new double[] {0, year};
	double[] amplBounds = new double[] {0, 5*Math.sqrt(ts.varianceInBinHeights())};
	double[] yOffsetBounds = new double[] {ts.minBinHeight(), ts.maxBinHeight()};
	return new double[][]{periodBounds, phaseBounds, amplBounds, yOffsetBounds};
    }
    
    private static void calculateModelPredictions(TimeSeries segmentExit, double[] modelExit, TimeSeries segmentEntry, double[] modelEntry) throws Exception {
	//  Compute envelopes
	double sigmaExit = getSigmaFromResiduals(segmentExit.getBinHeights(), modelExit);
	double sigmaEntry = getSigmaFromResiduals(segmentEntry.getBinHeights(), modelEntry);
	int nFittedRevs = modelExit.length;
	double[] revNum = new double[nFittedRevs];
	double[] safeExit = new double[nFittedRevs];
	double avgExit = BasicStats.getMean(modelExit);
	double sigmaExitFraction = sigmaExit/avgExit;
	double[] safeEntry = new double[nFittedRevs];
	double avgEntry = BasicStats.getMean(modelEntry);
	double sigmaEntryFraction = sigmaEntry/avgEntry;
	for ( int i=0; i < nFittedRevs; i++ ) {
	    revNum[i] = (double)startRev + i;
	    safeExit[i] = modelExit[i] + modelExit[i]*1.5*sigmaExitFraction;
	    safeEntry[i] = modelEntry[i] + modelEntry[i]*1.5*sigmaEntryFraction;
	}
	double lastRev = revNum[nFittedRevs-1];
	String[] header = new String[] {
	    "DEV /XS",
	    "LAB F",
	    "TIME OFF",
	    "LINE STEP",
	    "VIEW 0.1 0.2 0.9 0.8",
	    "CS 1.3",
	    "LW 4",
	    "LAB X Revolution Number",
	    "LAB Y Altitude (km)",
	    "R Y 2.5e4 8.2e4"
	};
	AsciiDataFileWriter out = new AsciiDataFileWriter(results+"envelopeExit_"+startRev+"-"+nRevs+".qdp");
	out.writeData(header, revNum, segmentExit.getBinHeights(), modelExit, safeExit);
	out = new AsciiDataFileWriter(results+"envelopeEntry_"+startRev+"-"+nRevs+".qdp");
	out.writeData(header, revNum, segmentEntry.getBinHeights(), modelEntry, safeEntry);

	//  Compute predictions
	int nChangesPerYear = 137;
	int timeInRevsBetweenChanges = nRevsPerYear/nChangesPerYear;
	int nRevsToPredict = (int)1.5*nRevsPerYear;
	int nPredictionPoints = nRevsToPredict/timeInRevsBetweenChanges;
	double[] predictionExit = new double[nPredictionPoints];
	double[] predictionEntry = new double[nPredictionPoints];
	double[] safePredictionExit = new double[nPredictionPoints];
	double[] safePredictionEntry = new double[nPredictionPoints];
	double[] binCentres = segmentExit.getBinCentres();
	double lastRevTime = binCentres[binCentres.length-1];
	double dt = BasicStats.getMean(segmentExit.getBinWidths());
	double timeBetweenChanges = dt*timeInRevsBetweenChanges;
	double[] halfWidths = new double[nPredictionPoints];
	double[] times = new double[nPredictionPoints];
	double[] revTimes = new double[nPredictionPoints];
	double[] nan = new double[nPredictionPoints];
	int i=0;
	
	//  Print out model/safe predictions to file
	int bufferSize = 256000;
	String filename = results+"predicted_heights.txt";
  	PrintWriter printWriter = new PrintWriter(new BufferedWriter(new FileWriter(filename), bufferSize));
	printWriter.println("REV_RANGE ENTRY/SAFE EXIT/SAFE");
	while ( i < predictionExit.length ) {
	    times[i] = lastRevTime + i*timeBetweenChanges;
	    //halfWidths[i] = dt/2;
	    revTimes[i] = lastRev + i*timeInRevsBetweenChanges;
	    nan[i] = Double.NaN;
	    predictionExit[i] = fittedFunctionExit.value(new double[] {times[i]});
	    safePredictionExit[i] = predictionExit[i] + predictionExit[i]*1.5*sigmaExitFraction;
	    predictionEntry[i] = fittedFunctionEntry.value(new double[] {times[i]});
	    safePredictionEntry[i] = predictionEntry[i] + predictionEntry[i]*1.5*sigmaEntryFraction;
	    double start = revTimes[i] - timeInRevsBetweenChanges/2;
	    double end = revTimes[i] + timeInRevsBetweenChanges/2;
	    printWriter.println((int)start+"-"+(int)end+" "+(int)predictionEntry[i]+"/"+(int)safePredictionEntry[i]
				+" "+(int)predictionExit[i]+"/"+(int)safePredictionExit[i]);
	    i++;
	}
	printWriter.close();

	// EXIT MODEL
	int minY = (int)Math.floor(BasicStats.getMedian(modelExit) - 4*sigmaExit);
	int maxY = (int)Math.ceil(BasicStats.getMedian(modelExit) + 4*sigmaExit);
	int minX = startRev - 10;
	int maxX = (int)lastRev + nRevsToPredict;
	header = new String[] {
	    "DEV /XS",
	    "LAB F",
	    "TIME OFF",
	    "LINE STEP",
	    "VIEW 0.1 0.2 0.9 0.8",
	    "CS 1.3",
	    "LW 4",
	    "LAB X Revolution Number",
	    "LAB Y Altitude (km)",
	    "R Y "+minY+" "+maxY,
	    "R X "+minX+" "+maxX
	};
	out = new AsciiDataFileWriter(results+"predictionExit_"+minX+"-"+maxX+".qdp");
	out.writeData(header, revTimes, nan, predictionExit, safePredictionExit);

	// ENTRY MODEL
	minY = (int)Math.floor(BasicStats.getMedian(modelEntry) - 4*sigmaEntry);
	maxY = (int)Math.ceil(BasicStats.getMedian(modelEntry) + 4*sigmaEntry);
	minX = startRev - 10;
	maxX = (int)lastRev + nRevsToPredict;
	header = new String[] {
	    "DEV /XS",
	    "LAB F",
	    "TIME OFF",
	    "LINE STEP",
	    "VIEW 0.1 0.2 0.9 0.8",
	    "CS 1.3",
	    "LW 4",
	    "LAB X Revolution Number",
	    "LAB Y Altitude (km)",
	    "R Y "+minY+" "+maxY,
	    "R X "+minX+" "+maxX
	};
	out = new AsciiDataFileWriter(results+"predictionEntry_"+minX+"-"+maxX+".qdp");
	out.writeData(header, revTimes, nan, predictionEntry, safePredictionEntry);
    }
    
    private static TimeSeries[] makeTimeSeries() throws Exception {
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
	    if ( !Double.isNaN(binHeightsExit[i]) ) {
		binHeightsExitList_minBins.add(binHeightsExit[i]);
		binEdgesList_minBins.add(binEdges[2*i]);
		binEdgesList_minBins.add(binEdges[2*i+1]);
	    }
	    if ( !Double.isNaN(binHeightsEntry[i]) ) {
		binHeightsEntryList_minBins.add(binHeightsEntry[i]);
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

    private static double getSigmaFromResiduals(double[] data, double[] model) throws Exception {
	int nFittedRevs = model.length;
	double[] residuals = new double[nFittedRevs];
	double[] fractionalResiduals = new double[nFittedRevs];
	for ( int i=0;  i < nFittedRevs; i++ ) {
	    residuals[i] = data[i] - model[i];
	    fractionalResiduals[i] = 100.*(data[i] - model[i])/data[i];
	}
	int nBins = residuals.length/5;
	IHistogram1D histoOfRes = Binner.makeHisto(residuals, nBins);
	IHistogram1D histoOfFracRes = Binner.makeHisto(fractionalResiduals, nBins);
	AsciiDataFileWriter out = new AsciiDataFileWriter(results+"residuals_"+startRev+"-"+nRevs+".qdp");
	out.writeHisto(histoOfFracRes, "Fractional Residuals (%)", true);
	//  Fit and get sigma
 	IAnalysisFactory af = IAnalysisFactory.create();
  	IFitFactory fitF   = af.createFitFactory();
 	IFitter fitter = fitF.createFitter("Chi2", "jminuit");
	IFitResult gaussFitResult = fitter.fit(histoOfRes, "g");
	return gaussFitResult.fittedParameter("sigma");
    }

    private static void readRadAndRevnoFiles() throws Exception {
	// read radiation data file downloaded from IDSC webpage
	AsciiDataFileReader in = new AsciiDataFileReader(data+"rad.dat");
	revNum = in.getDblCol(0);
	nRevs = revNum.length;
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
