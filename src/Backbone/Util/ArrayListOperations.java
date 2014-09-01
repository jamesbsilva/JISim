package Backbone.Util;

/**
*      @(#)     ArrayListOperations
*/

import java.util.ArrayList;

/**
*       ArrayListOperations calculates some basic values for arraylist.
* 
* <br>
* 
* @author James B. Silva <jbsilva @ bu.edu>
* @since 2013-11
*/
public class ArrayListOperations {
    public static double calculateMeanIntAL(ArrayList<Integer> arr){
        double sum = 0.0;
        for(int u = 0; u < arr.size();u++){
            sum += arr.get(u);
        }
        return (sum/((double)arr.size()));
    }
    public static double calculateStdMeanIntAL(ArrayList<Integer> arr){
        double sum = 0.0; double mean  = calculateMeanIntAL(arr);
        for(int u = 0; u < arr.size();u++){
            sum += Math.pow((arr.get(u) - mean),2.0);
        }
        sum = (sum)/((double)arr.size());
        sum = Math.sqrt(sum);
        sum = sum/(Math.sqrt(arr.size()));
        return sum;
    }
    public static double calculateMeanDblAL(ArrayList<Double> arr){
        double sum = 0.0;
        for(int u = 0; u < arr.size();u++){
            sum += arr.get(u);
        }
        return (sum/((double)arr.size()));
    }
    public static double calculateStdMeanDblAL(ArrayList<Double> arr){
        double sum = 0.0; double mean  = calculateMeanDblAL(arr);
        for(int u = 0; u < arr.size();u++){
            sum += Math.pow((arr.get(u) - mean),2.0);
        }
        sum = (sum)/((double)arr.size());
        sum = Math.sqrt(sum);
        sum = sum/(Math.sqrt(arr.size()));
        return sum;
    }
}