package Backbone.Util.Interactions;

/**
*    @(#)   InteractionMap 
*/  

import java.util.ArrayList;

/**
*       InteractionMap is a map which gives location of interacting sites.
* 
* <br>
* 
* @author      James B. Silva <jbsilva @ bu.edu>                 
* @since       2013-10
*/
public class InteractionMap {
    protected ArrayList<ArrayList<int[]>> map;
    
    public InteractionMap(){
        map = new ArrayList<ArrayList<int[]>>();
    }
    
    public ArrayList<int[]> getLocation(int u){
        return map.get(u);
    }

    public int getInteractingN(int u){
        return map.get(u).size();
    }
}

