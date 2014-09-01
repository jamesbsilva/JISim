package Triggers;

/**
* 
*    @(#) Trigger 
*/  
/** 
*      Basic Trigger interface to be implemented and used in simulations
*    to determine if a lattice needs to be triggered.
*  <br>
*  
* @author      James Silva <jbsilva @ bu.edu>                 
* @since       2011-09    
*/
public  interface  Trigger {
    public int getTcutoff();
    public boolean triggerNow();
    public boolean triggerNow(double m);
    public boolean triggerNow(double m,double e);
    public void reset();
    public void reset2();//optional reset 
    public void resetTime();
    public void setT(int t);
    public int getT();
    public int getTriggerTime();
}
