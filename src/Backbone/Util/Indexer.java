package Backbone.Util;

import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;

/**
*    @(#)   Indexer 
*/  

/**
*       Indexer holds functions reused for calculating indexes.
* 
* <br>
* 
* @author      James B. Silva <jbsilva @ bu.edu>                 
* @since       2013-10
*/
public class Indexer {
    public static int getX(int loc, int L){
        return (loc%L);
    }
    public static int getY(int loc, int L){
        return ((int)((double)loc/((double)L))%L);
    }
    public static int getZ(int loc, int L){
        return ((int)((double)loc/((double)(L*L)))%L);
    }
    public static int getIndex(int x, int y, int L){
        return (x+y*L);
    }
    public static int getIndex(int x, int y,int z, int L){
        return (x+y*L+z*L*L);
    }
    
    public static String getSignificant(double value, int sigFigs) {
        MathContext mc = new MathContext(sigFigs, RoundingMode.DOWN);
        BigDecimal bigDecimal = new BigDecimal(value, mc);
        return bigDecimal.toPlainString();
    }
}