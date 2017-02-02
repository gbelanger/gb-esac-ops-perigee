
import java.io.FileWriter;
import java.io.BufferedWriter;
import java.io.PrintWriter;
import gb.esac.io.AsciiDataFileWriter;
import gb.esac.tools.MinMax;
import gb.esac.tools.BasicStats;
import gb.esac.timeseries.TimeSeries;
import hep.aida.IFitResult;
import hep.aida.IFitter;
import hep.aida.IFitFactory;
import hep.aida.IAnalysisFactory;
import hep.aida.IFunction;
import hep.aida.IHistogram1D;
import gb.esac.binner.Binner;
import hep.aida.IAxis;

/**

This class calculates the model predictions for the exit and entry altitudes.
It also produces the plots of the safe envelopes and prediction, the
histogram of residuals, and the predictions in a text file.

@author G. Belanger, ESA, ESAC
@version 2017 Feb

 **/

public class ModelPredictionCalculator {
    
    private static String results = "../results/current/";    

    static void calculateModelPredictions(int startRev, int nRevs, TimeSeries segmentExit,  TimeSeries segmentEntry, IFunction fittedFunctionExit, IFunction fittedFunctionEntry, double[] modelExit, double[] modelEntry, double safetyFactorInSigmas) throws Exception {
	//  Compute envelopes
	double sigmaExit = getSigmaFromResiduals(startRev, nRevs, segmentExit.getBinHeights(), modelExit);
	double sigmaEntry = getSigmaFromResiduals(startRev, nRevs, segmentEntry.getBinHeights(), modelEntry);
	int nFittedRevs = modelExit.length;
	double[] revNum = new double[nFittedRevs];
	double[] safeExit = new double[nFittedRevs];
	double[] safeEntry = new double[nFittedRevs];
	double normFactorExit = BasicStats.getMedian(modelExit);
	double normFactorEntry = BasicStats.getMedian(modelEntry);
	for ( int i=0; i < nFittedRevs; i++ ) {
	    revNum[i] = (double)startRev + i;
	    safeExit[i] = modelExit[i] + safetyFactorInSigmas * sigmaExit * (modelExit[i]/normFactorExit); 
	    safeEntry[i] = modelEntry[i] + safetyFactorInSigmas * sigmaEntry * (modelEntry[i]/normFactorEntry);
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
	int nRevsPerYear = GetPerigeeHeights.nRevsPerYear;
	int timeInRevsBetweenChanges = nRevsPerYear/nChangesPerYear;
	int nRevsToPredict = (int)1.5*nRevsPerYear;
	int nPredictionPoints = nRevsToPredict/timeInRevsBetweenChanges;
	double[] predictedExit = new double[nPredictionPoints];
	double[] predictedEntry = new double[nPredictionPoints];
	double[] predictedSafeExit = new double[nPredictionPoints];
	double[] predictedSafeEntry = new double[nPredictionPoints];
	double[] binCentres = segmentExit.getBinCentres();
	double lastRevTime = binCentres[binCentres.length-1];
	double dt = BasicStats.getMean(segmentExit.getBinWidths());
	double timeBetweenChanges = dt*timeInRevsBetweenChanges;
	double[] times = new double[nPredictionPoints];
	double[] revTimes = new double[nPredictionPoints];
	double[] nan = new double[nPredictionPoints];
	int i=0;
	
	//  Print out model/safe predictions to file
	int bufferSize = 256000;
	String filename = results+"predicted_heights.txt";
  	PrintWriter printWriter = new PrintWriter(new BufferedWriter(new FileWriter(filename), bufferSize));
	printWriter.println("REV_RANGE ENTRY/SAFE EXIT/SAFE");
	while ( i < predictedExit.length ) {
	    times[i] = lastRevTime + i*timeBetweenChanges;
	    revTimes[i] = lastRev + i*timeInRevsBetweenChanges;
	    nan[i] = Double.NaN;
	    //
	    predictedExit[i] = fittedFunctionExit.value(new double[] {times[i]});
	    predictedEntry[i] = fittedFunctionEntry.value(new double[] {times[i]});
	    //
	    predictedSafeExit[i] = predictedExit[i] + safetyFactorInSigmas * sigmaExit * (predictedExit[i]/normFactorExit);
	    predictedSafeEntry[i] = predictedEntry[i] + safetyFactorInSigmas * sigmaEntry * (predictedEntry[i]/normFactorEntry);
	    //
	    double start = revTimes[i] - timeInRevsBetweenChanges/2;
	    double end = revTimes[i] + timeInRevsBetweenChanges/2;
	    printWriter.println((int)start+"-"+(int)end+" "+(int)predictedEntry[i]+"/"+(int)predictedSafeEntry[i]+" "+(int)predictedExit[i]+"/"+(int)predictedSafeExit[i]);
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
	out = new AsciiDataFileWriter(results+"predictedExit_"+minX+"-"+maxX+".qdp");
	out.writeData(header, revTimes, nan, predictedExit, predictedSafeExit);

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
	out = new AsciiDataFileWriter(results+"predictedEntry_"+minX+"-"+maxX+".qdp");
	out.writeData(header, revTimes, nan, predictedEntry, predictedSafeEntry);
    }

    private static double getSigmaFromResiduals(int startRev, int nRevs, double[] data, double[] model) throws Exception {
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
	//  Fit and get sigma
 	IAnalysisFactory af = IAnalysisFactory.create();
  	IFitFactory fitF   = af.createFitFactory();
 	IFitter fitter = fitF.createFitter("Chi2", "jminuit");
	IFitResult fitResult_fracRes = fitter.fit(histoOfFracRes, "g");
	IFitResult fitResult_res = fitter.fit(histoOfRes, "g");
	//  Plot histogram and fitted function
	IFunction fittedFunction = fitResult_fracRes.fittedFunction();
	IAxis axis = histoOfFracRes.axis();
	nBins = axis.bins();
	double[] fittedCurve = new double[nBins];
	for ( int i=0; i < nBins; i++ ) {
	    double binCentre = axis.binCenter(i);
	    fittedCurve[i] = fittedFunction.value(new double[] {binCentre});
	}
	AsciiDataFileWriter out = new AsciiDataFileWriter(results+"residuals_"+startRev+"-"+nRevs+".qdp");
	out.writeHisto(histoOfFracRes, fittedCurve, "Fractional Residuals (%)", true);

	return fitResult_res.fittedParameter("sigma");
    }

}
