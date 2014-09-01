package Backbone.System;

/**
* 
*   @(#) LatticeInt
*/  
import java.awt.image.BufferedImage;

/**
*      LatticeInt interface to be implemented with all compatible lattice 
*  classes. 
* 
* 
* <br>
* @author      James Silva <jbsilva @ bu.edu>                 
* @since       2011-09    
*/
public interface Lattice extends SimSystem {    
    /**
    *         getLength which gets the length of the lattice. 1D should be the same as N.
    * 
    */ 
    public int getLength();
    /**
    *         getDimension gives the dimensionality of the lattice
    */ 
    public int getDimension();
    /**
    *         getGeo gives the geometry of the lattice
    */ 
    public int getGeo();
    /**
    *         getN gives the size of the lattice
    */ 
    public int getN();
    /**
    *         getRange gives the interaction range of the lattice
    */ 
    public int getRange();
    /**
    *         getMagnetization gives the magnetization of the lattice
    */ 
    public int getMagnetization();
    /**
    *         getMagStaggered gives the staggered magnetization of the lattice
    */ 
    public int getMagStaggered();
    /**
    *         setFixedLatticeValues should set all fixed values in the lattice.
    */ 
    public void setFixedLatticeValues();
    /**
    *         isThisFixed should return true if the lattice coordinate is fixed.
    * 
    *  @param i - i coordinate
    *  @param j - j coordinate
    *  @param k - k coordinate 
    */ 
    public boolean isThisFixed(int i,int j,int k);
    /**
    *         getNeighSum should return the sum of the spins within the interaction
    *   range centered by the coordinates given. 3d case
    * 
    *  @param i - i coordinate
    *  @param j - j coordinate
    *  @param k - k coordinate 
    */ 
    /**
    *         getNinRnage should return the amount of spins within the interaction 
    *   range. Useful for metropolis algorithm in long range case.
    */ 
    public int getNinRange();
    
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
    *         getSystemImg should return an image of the lattice
    */ 
    public BufferedImage getSystemImg();
    /**
    *         makeVideo should set the lattice to begin making video which should
    *   require updating of the lattice image.
    */
    public void makeVideo();
    
    /**
    *         getFixSpinVal should return the average value of the fixed spins
    * 
    */ 
    public int getFixSpinVal();
    /**
    *         getNFixed gives the amount of fixed spins in the lattice
    */ 
    public int getNFixed();
}


