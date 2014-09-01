package AnalysisAndVideoBackend;
/**
 *         LinearRegression
 */
import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;
/** 
 *  Takes a set of data values and computes the
 *  best fit (least squares) line y  = mx + b through the set of values.
 *  Also computes the correlation coefficient and the standard deviation
 *  of the regression coefficients. This class is used in a trigger hence
 *  a trigger time function.
 * 
 * @author      James Silva <jbsilva @ bu.edu>                 
 * @since       2011-12
 */

public class LinearRegression { 

    private double[] dataX;
    private double[] dataY;
    private double m; 
    private double b;
    private double rsq;
    private double merr;
    private double berr;
    private boolean dataXset=false;
    private boolean dataYset=false;
    
    
    /**
    *       getParameter1 returns the m (slope) of linear fit
    * 
    * @return slope of linear fit
    */
    public double getParameter1(){return m;}
    
    /**
    *       getParameter2 returns the b (intercept) of linear fit
    * 
    * @return intercept of linear fit
    */
    public double getParameter2(){return b;}
    
    /**
    *       getRsquared returns the r squared value of fit
    * 
    * @return r squared of fit
    */
    public double getRsquared(){return rsq;}
    
    /**
    *       getParameter1Error returns the m (slope) error of linear fit
    * 
    * @return slope error of linear fit
    */
    public double getParameter1Error(){return merr;}
    
    /**
    *       getParameter2 returns the b (intercept) error of linear fit
    * 
    * @return error on intercept of linear fit
    */
    public double getParameter2Error(){return berr;}
    
    
    /**
    *       setDataX sets a list of doubles as the x coordinate data.
    * 
    * @param que - x coordinate data as list of doubles
    */
    public void setDataX( List<Double>  que){dataX=processList(que);dataXset=true;}
    /**
    *       setDataY sets a list of doubles as the y coordinate data.
    * 
    * @param que - y coordinate data as list of doubles
    */
    public void setDataY( List<Double>  que){dataY=processList(que);dataYset=true;}
    
    
    public void readInData(int size){
        
        // Your standard, unsynchronized list
        List<Double> xtemp = new ArrayList<Double>();
        List<Double> ytemp = new ArrayList<Double>();
        
        Scanner scanner;
        int i=0;
	try {
            scanner = new Scanner(new File("x.txt"));

	    while(scanner.hasNextDouble()) {
	    double x = scanner.nextDouble();
            xtemp.add(x);
           
            }
      
	} catch (FileNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
	}
         
        i=0;
	try {
            scanner = new Scanner(new File("y.txt"));

	    while(scanner.hasNextDouble()) {
	    double y = scanner.nextDouble();
            ytemp.add(y);
            }
            
            
            setDataX(xtemp);
            setDataY(ytemp);      
	} catch (FileNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
	}
        
        
        
    
    }
    
    /**
    *       calculateTriggerTime calculates the time where the fluctuations 
    *   where stable as x for y=0.5 in the fit
    * 
    *   @return time resulting from fit for y=0.5
    */
    public double calculateTriggerTime(){
    return ((0.5-b)/m);
    }
    
    /**
    *      processList takes a list of doubles and makes it into a double array
    * 
    * @param que - list of doubles to process
    * @return 
    */
    private double[] processList(List<Double>  que){
        double[] arr = new double[que.size()];
        int arrSize = que.size();
        for(int i=0;i<arrSize;i++){
        arr[i] = que.remove(0);
        }
        return arr;
    }
    
    
    /**
    *       performRegression performs linear regression on the array of data 
    */
    public void performRegression() { 

        if(!dataXset){
            System.out.println("SET YOUR X COORDINATE DATA USING SETDATAX");
            return;}
        
        if(!dataYset){
            System.out.println("SET YOUR Y COORDINATE DATA USING SETDATAY");
            return;}
        
        int sizeProper = 0;
        //  compute xMean and yMean
        double sumX = 0.0, sumY = 0.0, sumX2 = 0.0,sumY2=0.0;
        for(int i=0;i<dataX.length;i++) {
            sumX  += dataX[i];
            sumX2 += dataX[i] * dataX[i];
            sumY  += dataY[i];
            sumY2 += dataY[i] * dataY[i];
        }
        
        double xMean = sumX / dataX.length;
        double yMean = sumY / dataX.length;

        
        // second pass for ijMeans
        double xxMean = 0.0, yyMean = 0.0, xyMean = 0.0;
        for (int i = 0; i < dataX.length; i++) {
            xxMean += (dataX[i] - xMean) * (dataX[i] - xMean);
            yyMean += (dataY[i] - yMean) * (dataY[i] - yMean);
            xyMean += (dataX[i] - xMean) * (dataY[i] - yMean);
        }
        m = xyMean / xxMean;
        b = yMean - m * xMean;
        
        // print results
        //System.out.println("y   = " + m+ " * x + " + b);

        // analyze results
        int df = dataX.length - 2;
        double resSumSq = 0.0;      // residual sum of squares
        double regSumSq = 0.0;      // regression sum of squares
        
        for (int i = 0; i < dataX.length; i++) {
            double fit =m*dataX[i] + b;
            resSumSq += (fit - dataY[i]) * (fit - dataY[i]);
            regSumSq += (fit - yMean) * (fit - yMean);
        }
        
        
        
        double rSquared    = regSumSq / yyMean;
        double std  = resSumSq / df;
        double stdErrM = std / xxMean;
        double stdErrB = std/dataX.length + xMean*xMean*stdErrM;
        
        rsq = rSquared;
        merr = stdErrM;
        berr = stdErrB;
        
    }
 
    // test the class
    public static void main(String[] args) {
        LinearRegression reg = new LinearRegression();
        reg.readInData(20);
        reg.performRegression();
        System.out.println("y = 1 at x = "+reg.calculateTriggerTime());
    }
}
 
