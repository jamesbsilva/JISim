package GPUBackend;
/**
 *   @(#) SciKitSimInterface 
 *
 */

import Backbone.Algo.IsingMC;
import Backbone.Algo.MCAlgo;
import java.util.Random;
import scikit.jobs.Control;
import scikit.jobs.SimulationCL;

/**
 *      SciKitSimInterface deals with handling interactions with SciKit Sim calls.
 *  
 *
 *
 *   <br>
 *
 * @author James B. Silva <jbsilva@bu.edu>
 * @since May 2012
 */
public final class SciKitSimInterface extends SimulationCL{
    MCAlgo simIsing =  null;
    
    public SciKitSimInterface(MCAlgo s){
        simIsing = s;
    }
    
    @Override
    public void load(Control c) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void run() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void animate() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void clear() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void oneCalcStep() {
        if(simIsing != null){simIsing.doOneStep();}
    }
}

