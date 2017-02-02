
import hep.aida.IFitResult;
import hep.aida.IFitData;
import hep.aida.IFitter;
import hep.aida.IDataPointSet;
import hep.aida.IFitFactory;
import hep.aida.IDataPointSetFactory;
import hep.aida.ITree;
import hep.aida.IAnalysisFactory;
import gb.esac.aida.functions.SineFunction;
import hep.aida.IFunction;

import org.apache.log4j.Logger;

/**

This class is used to fit sinusoidally shapped functions to data.

@author G. Belanger, ESA, ESAC
@version 2017 Feb

 **/

public class SinusoidModelProvider {

    private static Logger logger  = Logger.getLogger(SinusoidModelProvider.class);

    public static IFunction getModel(double[] x, double[] y, double[] initialParValues, double[][] bounds, String label) {
	
	logger.info("- "+label+" -");

	// Define sine function with initial parameter values
	String[] parNames = new String[] {"period", "phase", "amplitude", "yOffset"};
	IFunction function = new SineFunction("sine");

	//  Print out the parameter names and values
	logger.info("Initial:");
	for ( int i=0; i < parNames.length; i++ ) {
	    logger.info("  "+parNames[i]+" = "+initialParValues[i]);
	    function.setParameter(parNames[i], initialParValues[i]);
	}	    

	//  Set up factories for fitting the data
	IAnalysisFactory af = IAnalysisFactory.create();
	ITree tree = af.createTreeFactory().create();
	IDataPointSetFactory dpsf = af.createDataPointSetFactory(tree);
  	IFitFactory fitF   = af.createFitFactory();

	//  Print available fit methods and engines
	boolean showFitMethods = false;
	if ( showFitMethods ) {
	    String[] fitMethods = fitF.availableFitMethods();
	    logger.info("Fit Methods");
	    for ( int i=0; i < fitMethods.length; i++ ) {
		logger.info(fitMethods[i]);
	    }
	    String[] fitEngines = fitF.availableFitEngines();
	    logger.info("Fit Engines");
	    for ( int i=0; i < fitEngines.length; i++ ) {
		logger.info(fitEngines[i]);
	    }
	}
	
	//  Create 2D IDataPointSet
	int n = x.length;
 	IDataPointSet dps = dpsf.create("dps", "Data", 2);
	int k = 0;
	for ( int i = 0; i < n; i++ ) {
	    if ( !Double.isNaN(x[i]) && !Double.isNaN(y[i]) ) {
		dps.addPoint();
		dps.point(k).coordinate(0).setValue(x[i]);
		dps.point(k).coordinate(1).setValue(y[i]);
		k++;
	    }
	}

	//  Create fitter
	IFitter fitter = fitF.createFitter("leastsquares", "jminuit");  //  Better
 	IFitData data = fitF.createFitData();
	int xCoordIndex = 0;
	int yCoordIndex = 1;
 	data.create1DConnection(dps, xCoordIndex, yCoordIndex);

	//  Define the bounds to constrain the parameter space
	fitter.fitParameterSettings("period").setBounds(bounds[0][0], bounds[0][1]);
	fitter.fitParameterSettings("phase").setBounds(bounds[1][0], bounds[1][1]);
	fitter.fitParameterSettings("amplitude").setBounds(bounds[2][0], bounds[2][1]);
	fitter.fitParameterSettings("yOffset").setBounds(bounds[3][0], bounds[3][1]);
	
	//  Do the fit
 	IFitResult fitResult = fitter.fit(data, function);
 	IFunction fittedFunction = fitResult.fittedFunction();
 	double[] fittedParValues = fitResult.fittedParameters();
	if ( fitResult.isValid() ) {
	    logger.info("Fitted:");
	    for ( int i=0; i < fittedParValues.length; i++ ) {
		logger.info("  "+parNames[i]+" = "+fittedParValues[i]);
	    }
	}
	else {
	    logger.warn("Fit result is NOT valid");
	}
	fitter.resetParameterSettings();

	return fittedFunction;
    }


}
