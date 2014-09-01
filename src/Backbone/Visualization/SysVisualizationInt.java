package Backbone.Visualization;

/**
* 
*     @(#)   SysVisualizationInt
*/  

import java.awt.image.BufferedImage;
import java.util.ArrayList;

/** 
*   SysVisualizationInt is a visualization of a system
*  <br>
* 
*  @param simSys - the simSys for step
* 
* 
* @author      James Silva <jbsilva @ bu.edu>                 
* @since       2013-05    
*/
public interface SysVisualizationInt extends SysVisualization{
    
    public BufferedImage getImageOfVis();

    /**
    * initializeImg creates the image layer of the system
    *
    */
    public void initializeImg();

    /**
    * updateImg just updates the image layer instead of redrawing the layer.
    * Change time text of image
    *
    * @param i - i coordinate
    * @param j - j coordinate
    * @param spin - new spin value
    * @param t - new time in image
    *
    */
    public void updateImg(int i,int dat, int t);

    /**
    * spinDraw completely redraws the image layer based on the current
    * configuration of the system.
    */
    public void spinDraw(int[] data);
    
    /**
    * spinDraw completely redraws the image layer based on the current
    * configuration of the system.
    */
    public void spinDraw(ArrayList<Integer> data);

}