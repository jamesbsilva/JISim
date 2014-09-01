/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package Backbone.System;

/**
 * 
 *   @(#) Lattice
 */  
import java.awt.image.BufferedImage;

/**
 *      Lattice interface to be implemented with all compatible lattice 
 *  classes. 
 * 
 * 
 * <br>
 * @author      James Silva <jbsilva @ bu.edu>                 
 * @since       2011-09    
 */
public interface SimSystem  {

    
    public int getInstId();
    public int getNFixed();
    public int getN();
    /**
    *         getLatticeImg should return an image of the lattice
    */ 
    public BufferedImage getSystemImg();
    /**
    *         makeVideo should set the system to begin making video which should
    *   require updating of the system image.
    */
    public void makeVideo();
    /**
    *         getMagnetization gives the magnetization of the lattice
    */ 
    public int getMagnetization();


}


