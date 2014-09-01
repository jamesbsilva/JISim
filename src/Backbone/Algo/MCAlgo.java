/**
 * 
 *   @(#) MCAlgo
 */  

package Backbone.Algo;

import Backbone.System.SimSystem;

 /**  
 *      Random MC algo general version.
 *  <br>
 * 
 * @author      James Silva <jbsilva @ bu.edu>                 
 * @since       2013-04    
 */


public interface MCAlgo {

     /**
    *        doOneStep should do a single monte carlo step.
    */ 
    public void doOneStep();
    /**
    *         setSeed should set the random number seed to the value given
    * 
    *  @param seed1 - random number seed
    */ 
    public void setSeed(int seed1);
    /**
    *         setRun should set the number of the current run to the value given.
    *
    *  @param cr - current run value 
    */
    public void setRun(int cr);
    /**
    *         resetSimulation should reset all simulation parameters.
    */
    public void resetSimulation();
    /**
    *         resetSimulation should reset all simulation parameters but 
    *  set seed to the given value.
    * 
    *  @param seed1 - random number seed
    */ 
    public void resetSimulation(int seed1);
    /**
    *         resetSimulation should reset all simulation parameters, set
    *   all lattice values to lattice at the time given, and set the 
    *   seed to the given value
    * 
    *  @param time - time to set lattice to
    *  @param seed1 - random number seed
    */ 
    public void resetSimulation(int time , int seed1);
     /**
    *   getSeed should return the current random number seed
    * 
    */
    public int getSeed();
    /**
    *   getRun should return the current run number
    * 
    */
    public int getRun();
    /**
    *   getSimSystem should return the system object
    * 
    */
    public SimSystem getSimSystem();
    
    /**
    *         getConfigRange should save the configurations ins the given range
    * 
    *  @param tInitial - initial time in range
    *  @param tFinal - final time in range
    */ 
    public void getConfigRange(int tInitial, int tFinal);

    /**
    *         setMeasuringStartTime sets the time to begin all measurement of all variables
    *   by a system measurer.
    * 
    *  @param tin - initial measurement time of all data
    */ 
    public void setMeasuringStartTime(int tin);
    
}
