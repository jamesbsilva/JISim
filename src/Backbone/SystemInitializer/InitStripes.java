package Backbone.SystemInitializer;
/**
 * 
 *      @(#)  InitStripes
 */  

import Backbone.System.Lattice;
import Backbone.System.LatticeMagInt;
import Backbone.System.SimSystem;

/**  
 *   Class for initializing a system into stripes.
 *  
 *  <br>
 * 
 * @author      James Silva <jbsilva @ bu.edu>                 
 * @since       2013-06    
 */
public class InitStripes{
    public LatticeMagInt initializeSys(LatticeMagInt lat, int L, int R, double angle, int s0, int dim){    
        int i,j,k;
        int offset,domain;
        for(int u = 0; u < Math.pow(L, dim);u++){
            i = u%L;
            j = (int)((double)u/(double)L)%L;
            //k = (int)(u/((double)L*L))%L;
            offset = (int) (j*Math.cos(angle*Math.PI/180));
            domain = (int)((double)(i+offset)/(double)R);
            if(domain%2 == 0){
                lat.setValue(u, -1*s0, 0);
            }else{
                lat.setValue(u, 1*s0, 0);
            }
        }
        return lat;
    }  
}