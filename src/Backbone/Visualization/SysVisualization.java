package Backbone.Visualization;

/**
* 
*     @(#)   SysVisualization
*/  

import java.awt.image.BufferedImage;
import java.util.ArrayList;

/** 
*   SysVisualization is a visualization of a system
*  <br>
* 
*  @param simSys - the simSys for step
* 
* 
* @author      James Silva <jbsilva @ bu.edu>                 
* @since       2013-05    
*/
public interface SysVisualization{
    
    public BufferedImage getImageOfVis();

    /**
    * initializeImg creates the image layer of the system
    *
    */
    public void initializeImg();
    
}