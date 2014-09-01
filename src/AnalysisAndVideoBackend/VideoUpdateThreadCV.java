package AnalysisAndVideoBackend;

/**
 * 
 *     @(#)   VideoUpdateThread
 */  
 

import Backbone.System.LatticeMagInt;
import Backbone.System.SimSystem;
import java.util.concurrent.Callable;

/** 
 *   Thread for a updating a video with a frame
 *  <br>
 * 
 *  @param simSys - the simSys for step
 *  @param vid - video maker class for current video
 * 
 * 
 * @author      James Silva <jbsilva @ bu.edu>                 
 * @since       2012-01    
 */


public class VideoUpdateThreadCV implements Callable<Boolean>{
    private volatile VideoMakerCV vid;
    private volatile SimSystem simSys;
    private boolean updateVid = true;
    
    /**  @param simSys - the simSys for step
    *  @param vid - video maker class for current video
    * 
    */ 
    public VideoUpdateThreadCV(SimSystem sys, VideoMakerCV vidm){
        simSys = sys;
        vid = vidm;
    }

    /**  @param simSys - the simSys for step
    *  @param vid - video maker class for current video
    * 
    */ 
    public VideoUpdateThreadCV(SimSystem sys, VideoMakerCV vidm, boolean up){
        simSys = sys;
        vid = vidm;
        updateVid = up;
    }
    
    /**
    *      call is the main logic that is going to be run by this callable object.
    *  For this object a frame is just written to the video and a boolean is 
    *  outputted based on if the video is complete.
    * 
    * @return true if video has not been written
    */
    @Override
    public Boolean call() {
        vid.addLatticeFrame(simSys.getSystemImg(),updateVid);
        return vid.isWritten();
    }    
}
