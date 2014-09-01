/**
 * 
 *   @(#) IsingMC
 */  

package Backbone.Algo;

import Backbone.System.LatticeMagInt;
import Backbone.System.SimSystem;
import Backbone.Util.MeasureIsingSystem;
import java.util.ArrayList;

 /**  
 *      Ising monte carlo simulation interface to be implemented 
 *   depending on which algorithm one decides to use.
 *  <br>
 * 
 * @author      James Silva <jbsilva @ bu.edu>                 
 * @since       2011-09    
 */


public interface IsingMC extends MCAlgo{

     /**
    *         setTrigger should initialize the trigger to be used in the simulation
    */ 
    public void setTrigger();
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
    *   getM should return the current magnetization.
    * 
    */
    public int getM();
    /**
    *   getMStag should return the current staggered magnetization.
    * 
    */
    public int getMStag();
    /**
    *   getEnergy should return the current energy.
    * 
    */
    public double getEnergy();
    /**
    *   getSpecificHeat should return the current specific heat.
    * 
    */
    public double getSpecificHeat();
    /**
    *   getSusceptibility should return the current susceptibility
    * 
    */
    public double getSusceptibility(); 
    /**
    *   getSusceptibility should return the current susceptibility
    * 
    */
    public ArrayList<Double> getSusceptibilityMult(); 
    /**
    *   getTriggerTime should return the time the trigger went off.
    * 
    */
    public int getTriggerTime();
    /**
    *         setTriggerOnOff allows turning triggers off
    * 
    * @param tr - on is true
    */ 
    public void setTriggerOnOff(boolean tr);
    /**
    *   flipField should flip the magnetic field
    * 
    */
    public void flipField();
    /**
    *   alignField should align the magnetic field with current state
    * 
    */
    public void alignField();
    /**
    *         changetT should set the temperature to the given value.
    * 
    *  @param temp - temperature to set simulation to
    */ 
    public void changeT(double temp);
    /**
    *         changetH should set the magnetic field to the given value.
    * 
    *  @param hnow - new magnetic field
    */ 
    public void changeH(double hnow);
    /**
    *         changetTandH should set the temperature to the given value and
    *   the magnetic field.
    * 
    * 
    *  @param temp - new temperature
    *  @param hnow - new field
    */ 
    public void changeTandH(double temp,double hnow);
    /**
    *         changetTFlipField should flip the field and set the temperature to
    *   the given value
    * 
    *  @param temp - new temperature
    */ 
    public void changeTFlipField(double temp);
    
    /**
    *   nucleated should return true if nucleation has occurred.
    * 
    */
    public boolean nucleated();
    /**
    *   getHfield should return the magnetic field value.
    * 
    */
    public double getHfield();
    /**
    *   getTemp should return the temperature value.
    * 
    */
    public double getTemperature();
    /**
    *       getSystemMeasurer should return the system measurer which contains system data
    *   and accumulators.
    * 
    */
    public MeasureIsingSystem getSystemMeasurer();
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
    /**
    *   getJinteraction should return the current strength of interaction.
    * 
    */
    public double getJinteraction();
        
}
