package Backbone.Util.Interactions;

/**
*    @(#)   SqInteractionMap 
*/  

import Backbone.Util.Interactions.InteractionMap;
import java.util.ArrayList;

/**
*       SqInteractionMap is a sub step of the simulation used to recreate the simulation.
* 
* <br>
* 
* @author      James B. Silva <jbsilva @ bu.edu>                 
* @since       2013-10
*/
public class TriInteractionMap extends InteractionMap {
    public TriInteractionMap(int R, int dim){
        super();
        if(R > 0){
            initMap(R,dim);
        // nearest neighbor
        }else{
            initMapNN(dim);
        }
    }

    private void initMapNN(int dim){
        ArrayList<int[]> tempIn  = new ArrayList<int[]>();
        int[] temp = new int[dim];
        temp[0] = -1;
        tempIn.add(temp);
        
        temp = new int[dim];
        temp[0] = +1;
        tempIn.add(temp);
        
        if(dim > 1){
            temp = new int[dim];
            temp[1] = -1;
            tempIn.add(temp);
            temp = new int[dim];
            temp[1] = +1;
            tempIn.add(temp);
        }
        
        if(dim > 2){
            temp = new int[dim];
            temp[2] = -1;
            tempIn.add(temp);
            temp = new int[dim];
            temp[2] = +1;
            tempIn.add(temp);
        }
        map.add(tempIn);
    }
    
    private void initMap(int R , int dim){
        if (dim != 2){
            throw new IllegalArgumentException("TriInteractionMap | Only 2D allowed.");
        }
        ArrayList<int[]> tempIn  = new ArrayList<int[]>();
        int subL = 2*R + 1 ;
        int m, n, p;
        int yL = ((dim > 1) ? subL : 1);
        int zL = ((dim > 2) ? subL : 1);
        for(int z = 0; z < zL; z++){
        for(int y = 0; y < yL; y++){
        for(int x = 0; x < subL; x++){
            m = -R+x; n = -R+y; p = -R+z;
            int[] rt = new int[dim];
            rt[0] = m;
            if(dim > 1){rt[1] = n;}
            if(dim > 2){rt[2] = p;}
            if(dim == 1){
                if(m == 0){continue;}
            }else if(dim == 2){
                if(m == 0 && n == 0){continue;}
            }else if(dim == 3){
                if(m == 0 && n == 0 && p == 0){continue;}
            }
            tempIn.add(rt);
        }}}
        map.add(tempIn);
    }
    @Override
    public ArrayList<int[]> getLocation(int u){
        return map.get(0);
    }
}
