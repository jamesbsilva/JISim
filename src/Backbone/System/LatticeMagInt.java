package Backbone.System;

/**
* 
*   @(#) LatticeMagInt
*/  
import java.awt.image.BufferedImage;

/**
*      LatticeMagInt interface to be implemented with all compatible lattice 
*  classes. 
* 
* 
* <br>
* @author      James Silva <jbsilva @ bu.edu>                 
* @since       2011-09    
*/
public interface LatticeMagInt extends LatticeInt {    
    /**
    *         setValue which updates the spin value in the lattice and the magnetization
    *       with update of time in image layer .
    * 
    *  @param i - i coordinate
    *  @param j - j coordinate
    *  @param k - k coordinate
    *  @param s - new spin value
    *  @param t - time 
    */ 
    public void setValue(int i,int j, int k,int s, int t);
    /**
    *         setValue which updates the spin value in the lattice and the magnetization
    *       with update of time in image layer .
    * 
    *  @param i - i coordinate
    *  @param j - j coordinate
    *  @param k - k coordinate
    *  @param s - new spin value
    *  @param t - time 
    */ 
    public void setValue(int i,int s, int t);
    /**
    *         getValue which gets the spin value in the lattice. 3d case
    * 
    *  @param i - i coordinate
    *  @param j - j coordinate
    *  @param k - k coordinate
    */ 
    public int getValue(int i,int j, int k);
    /**
    *         getValue which gets the spin value in the lattice. 2d case
    * 
    *  @param i - i coordinate
    *  @param j - j coordinate
    */ 
    public int getValue(int i,int j);
    /**
    *         getValue which gets the spin value in the lattice.
    * 
    *  @param i - i coordinate
    */ 
    public int getValue(int i);
    /**
    *         getLength which gets the length of the lattice
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
    public int getNeighSum(int i,int j,int k);
    /**
    *         getNeighSum should return the sum of the spins within the interaction
    *   range centered by the coordinates given. 2d case
    * 
    *  @param i - i coordinate
    *  @param j - j coordinate
    */ 
    public int getNeighSum(int i,int j);
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
    *         highFixedAndStableMatched should return if the sublattice index
    *   for the highest density stable spins and fixed spins match.
    */ 
    public boolean highFixedAndStableMatched();
    /**
    *         getHighMagIndex returns the index of the highest density of 
    *   stable spins
    */ 
    public int getHighMagIndex();
    /**
    *         getHighFixIndex returns the index of the highest density of 
    *   fixed spins
    */ 
    public int getHighFixIndex();
    /**
    *         getLastAvgMagSubLat should return the average magnetization of the 
    *   sublattices.
    */ 
    public double getLastAvgMagSubLat();
    /**
    *         getSubLatSum should return the sum of the sublattice in the index
    *   given
    * 
    *  @param i - index of sublattice
    */ 
    public int getSubLatSum(int i);
    
    /**
    *         getFixedSubLatSum should return the sum of the fixed spin amounts 
    *   in sublattice in the index given
    * 
    *  @param i - index of sublattice
    */ 
    public int getFixedSubLatSum(int i);
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


