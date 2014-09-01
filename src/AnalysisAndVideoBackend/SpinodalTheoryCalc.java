package AnalysisAndVideoBackend;

/**
*      @(#)     SpinodalTheoryCalc 
*/


/**
*       SpinodalTheoryCalc calculates some basic values from mean field spinodal theory.
* 
* <br>
* 
* @author James B. Silva <jbsilva @ bu.edu>
* @since 2013-11
*/
public class SpinodalTheoryCalc {
    public static double calcSpinodalField(double temperature, int range ){
        double Q = Math.pow(2*((double)range)+1,2)-1.0;
        double beta = 1.0/temperature; double J = 4.0/Q;
        return (temperature)*acosh(Math.sqrt(beta*J*Q))-Q*J*Math.sqrt((beta*Q*J-1.0)/(beta*Q*J));
    }

    private static double acosh(double x){ 
        return Math.log(x + Math.sqrt(x*x - 1.0)); 
    } 
    
    public static double calcDelHsScaled(double field,double  temperature,int range, double xper, boolean stable){
        return Math.abs((Math.abs(field) + ((stable) ? 4.0*xper: -4.0*xper) -
                (  Math.abs(calcSpinodalField(temperature,range))) )/Math.abs(calcSpinodalField(temperature,range)));
    }
    
    public static double calcDelHsRaw(double field,double  temperature,int range, double xper, boolean stable){
        return Math.abs(Math.abs(field) + ((stable) ? 4.0*xper: -4.0*xper) - (  Math.abs(calcSpinodalField(temperature,range))));
    }
    
    public static double calcGinzbergParam(double field,double  temperature,int range,double xper, boolean stable){
        return Math.pow(range, 2.0)*Math.pow(calcDelHsRaw(field,temperature,range,xper,stable), 2.0/2.0);
    }
    
    public static double calcCorrelationLength(double field,double  temperature,int range,double xper, boolean stable){
        return range*Math.pow(calcDelHsRaw(field,temperature,range,xper,stable), -1.0/4.0);
    }
    
    public static void main(String[] args){
        System.out.println("dh : "+SpinodalTheoryCalc.calcDelHsScaled(0.775, 1.8, 8, 0.05, true));
        System.out.println("Hs : "+SpinodalTheoryCalc.calcSpinodalField(1.8, 8));
        System.out.println("Correlation Length : "+SpinodalTheoryCalc.calcCorrelationLength(0.775, 1.8, 8, 0.05, true));
        //System.out.println("dh : "+SpinodalTheoryCalc.calcDelHsScaled(1.14, 1.8, 8, 0.0, true));
    }
    
}