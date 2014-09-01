/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package Backbone.System;

/**
 * 
 *   @(#) Network
 */  
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.atomic.AtomicIntegerArray;

/**
 *      Network interface to be implemented with all compatible lattice 
 *  classes. 
 * 
 * 
 * <br>
 * @author      James Silva <jbsilva @ bu.edu>                 
 * @since       2011-09    
 */
public interface Network extends SimSystem{

    
    
    /**
    *         setValue which updates the spin value in the lattice and the magnetization
    *       with update of time in image layer .
    * 
    *  @param i - i coordinate
    *  @param s - new spin value
    *  @param t - time 
    */ 
    public void setValue(int i,int s, int t);
    /**
    *         getValue which gets the spin value in the lattice. 3d case
    * 
    *  @param i - i coordinate
    */ 
    public int getValue(int i);
    /**
    *         getN gives the size of the lattice
    */ 
    public int getN();
    /**
    *         getMagnetization gives the magnetization of the lattice
    */ 
    public int getMagnetization();
    /**
    *         setFixedNetworkValues should set all fixed values in the lattice.
    */ 
    public void setFixedNetworkValues();
    /**
    *         initialize sets the lattice values to their initial values which is given as input.
    *
    *
    *  @param s - Initial lattice value; -2 for random values
    *
    */ 
    public void initialize(int s);
    /**
    *         isThisFixed should return true if the lattice coordinate is fixed.
    * 
    *  @param i - i coordinate
    */ 
    public boolean isThisFixed(int i);
    /**
    *         getNeighSum should return the sum of the spins within the interaction
    *   range centered by the coordinates given. 3d case
    * 
    *  @param i - i coordinate 
    */ 
    public int getNeighSum(int i);
    /**
    *         setInitialConfig should set all spins in the lattice to the values
    *   given by a file in the Config directory.
    *   
    *   @param t - time of config to set lattice to
    *   @param run - run of time to set lattice to
    *   @param  post - postfix or filename of lattice file
    */
    public void setInitialConfig(int t,int run,String post);
    /**
    *         getInstId returns the identification number of the instance of this
    *  program running. Useful for allowing each instance a directory to work in.
    */ 
    public int getInstId();
    /**
    *         getNFixed gives the amount of fixed spins in the lattice
    */ 
    public int getNFixed();
    
    public int getNetworkSeed();
    public int getTotalLinks();
    public void setNetworkSeed(int nseed);
    public HashMap<Integer,ArrayList<Integer>> getNetworkLinks();
    public HashMap<Integer,Integer> getNetworkValues();
    public AtomicIntegerArray getNetworkValuesArr();
    public HashMap<Integer,Boolean> getNetworkFixed();
    public void saveDegreeFrequency();

}


