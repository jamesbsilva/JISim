package AnalysisAndVideoBackend;

/**
 * 
 *   @(#) VideoMaker
 */ 

import Backbone.System.AtomicLattice;
import Backbone.Util.DirAndFileStructure;
import Backbone.Util.ParameterBank;
import java.awt.image.BufferedImage;
import com.googlecode.javacv.CanvasFrame;
import com.googlecode.javacv.FFmpegFrameRecorder;
import com.googlecode.javacv.FrameRecorder;
import com.googlecode.javacv.cpp.avutil;
import com.googlecode.javacv.cpp.opencv_core.IplImage;
import java.util.logging.Level;
import java.util.logging.Logger;

/**  A Video maker class for lattices that uses openCV.
 *  <br>
 *  @param seed - seed for this particular run
 *  @param max - max amount of frames
 *  @param type - image file type that will go in avi
 *  @param Name - optional add in to file name
 * 
 * @author      James Silva <jbsilva @ bu.edu>                 
 * @since       2013-04    
 */
public class VideoMakerCV{
    private boolean displayVid;
    private String vidName; private String vidFormat="avi";
    private boolean vidInitialized = false; private boolean videoWritten = false;
    private FFmpegFrameRecorder recorder; private CanvasFrame canvasFrame ;
    private int frameRate = 3;
    private boolean knownLength=true;
    private int framesCaptured = 0;
    private ParameterBank param; private DirAndFileStructure dir;
    private int maxFrames;
    
    public VideoMakerCV(int seed,int max,int time){
        dir = new DirAndFileStructure(); param = new ParameterBank();
        if(max==0){knownLength=false;}
        maxFrames = max;
        String fname = "vid";
        vidName = dir.getVideoDirectory()+fname+"-h-"+param.hField+"-L-"+param.L+"-seed-"+seed+
        "-t-"+time+"."+vidFormat;
        dir = null;
        System.out.println("Writing " + vidName);
    }

    public VideoMakerCV(int seed,int max,String Name){
        dir = new DirAndFileStructure(); param = new ParameterBank();
        if(max==0){knownLength=false;}
        maxFrames = max;
        String fname = "vid";
        vidName = dir.getVideoDirectory()+fname+"--"+Name+"-seed-"+seed+".avi";
        dir = null;
        System.out.println("Writing " + vidName);
    }

    
    public VideoMakerCV(int seed,int max){
        this(seed,max,"");
    }
    
    /*
     *      Color map
     *      s = 1 (white) 255 255 255
     *      s = -1 ( red dark )  204 0 0  r g b
     *      s = -1 fixed blue 0 0 255 
     *      s = 1  light yellow 255 255 51
     *      s = 0  green  0 204 0
     *      s = 0 fixed orange 204 102 0
     * 
     */
   /* 
    Scalar getColor(int s, boolean fix){
        Scalar col = new Scalar(96,96,96);
        if(!fix){
            if(s == 0){
                col = new Scalar(0,204,0);
            }else if(s == -1){
                col = new Scalar(0,0,204);
            }else if(s == 1){
                col = new Scalar(255,255,255);
            }
        }else{
            if(s == 0){
                col = new Scalar(0,102,204);
            }else if(s == -1){
                col = new Scalar(255,0,0);
            }else if(s == 1){
                col = new Scalar(51,255,255);
            }
        }
        return col;
    }
    public void drawSquare(Mat img, int x, int y, int r){
        
        Core.rectangle(img, new Point(x-r,y-r), new Point(x+r,y+r), new Scalar(0, 255, 0));
    }
    
    
    public Mat convertBufferedImageToMat(BufferedImage imgBuff){
        byte[] pixels = ((DataBufferByte) imgBuff.getRaster().getDataBuffer()).getData();
        Mat m = new Mat(imgBuff.getHeight(), imgBuff.getWidth(), CvType.CV_8UC3, new Scalar(0));
        m.put(0, 0, pixels);
        return m;
    }

    
    Scalar getColor(float b, float g, float r){
        return new Scalar(b,g,r);
    }
    */
   
    public IplImage convertBuffImageToJCVimg(BufferedImage imgBuff){
        IplImage img = IplImage.createFrom(imgBuff); 
        return img;
    }
    
    public void setVidName(String vid){
        vidName = vid;
    }
    public void setDisplayVideo(boolean watchVid){
        displayVid = watchVid;
    }
    
    public void initVideo(BufferedImage buffImg){
        initCanvas(buffImg,false);
        // record video
        //System.out.println("framerate = " + frameRate);
        recorder = new FFmpegFrameRecorder(vidName,  buffImg.getWidth(),buffImg.getHeight());
        recorder.setVideoCodec(13);
        //recorder.setFormat("mp4");
        recorder.setFormat(vidFormat); recorder.setPixelFormat(avutil.AV_PIX_FMT_YUV420P);
        recorder.setFrameRate(frameRate); recorder.setVideoBitrate(10 * 1024 * 1024);        
        try {
            recorder.start();    
        } catch (FrameRecorder.Exception ex) {
            Logger.getLogger(VideoMakerCV.class.getName()).log(Level.SEVERE, null, ex);
        }
        vidInitialized = true;
    }
    
    private void initCanvas(BufferedImage buffImg , boolean canvasOnly){
        if(displayVid){
            canvasFrame = new CanvasFrame("Simulation");
            canvasFrame.setCanvasSize(buffImg.getWidth(),buffImg.getHeight());
            if(canvasOnly){vidInitialized = true;}
        }
    }
    
    /**
    *      writeVideo writes the video in avi format to the video directory.
    */
    public void writeVideo(){
        if(recorder != null){
            try {
                recorder.stop();
            } catch (FrameRecorder.Exception ex) {
                Logger.getLogger(VideoMakerCV.class.getName()).log(Level.SEVERE, null, ex);
            }
            System.out.println("Done with Video");    
        }
        if(displayVid){
            canvasFrame.dispose();
        }
        videoWritten=true;
    }
    /**
    *      setFramesPerCapture sets the amount of frames to push into the video
    *  when a frame is attempted to be added to the video.
    * 
    *  @param fr - frames per each lattice add
    */
    public void setFramesPerCapture(int fr){frameRate=fr;}
    
    public void addLatticeFrame(BufferedImage buffImg,boolean vidup){
        if(framesCaptured > maxFrames && knownLength){
            // Do Nothing
        }else if(framesCaptured==maxFrames && knownLength){
            System.out.print("Filled Max Frames - "+maxFrames+"   -");
            writeVideo();
            framesCaptured++;framesCaptured++;
        }else{
            if(displayVid && !vidup){
                if(!vidInitialized){
                    initCanvas(buffImg,true);               
                } 
                canvasFrame.showImage(buffImg);
            }
            if(vidup){
                if(!vidInitialized){
                    initVideo(buffImg);
                }
                if(displayVid){
                    canvasFrame.showImage(buffImg);
                }
                IplImage grabbedImage = convertBuffImageToJCVimg(buffImg);

                try {
                    recorder.record(grabbedImage);
                } catch (FrameRecorder.Exception ex) {
                    Logger.getLogger(VideoMakerCV.class.getName()).log(Level.SEVERE, null, ex);
                }
                framesCaptured++;
            }
        }
    }
    
    public void addLatticeFrame(BufferedImage buffImg){
        addLatticeFrame(buffImg,true);
    }

    public void addLatticeFrameToShownOnly(BufferedImage buffImg){
        addLatticeFrame(buffImg,false);
    }
    
    /**
    *       isWritten returns a boolean based on whether a video has been written 
    *   or not
    * 
    * @return  true if video has not been written
    */
    public boolean isWritten(){return !videoWritten;}
    
    public static void main(String[] args) {
        VideoMakerCV vid = new VideoMakerCV(0,100,"test");
        vid.setFramesPerCapture(1000);
        vid.setDisplayVideo(true);        
        AtomicLattice lat = new AtomicLattice(1);
        lat.makeVideo();
        vid.addLatticeFrame(lat.getSystemImg());
        lat.setValue(3, 4, 0, -1);
        vid.addLatticeFrame(lat.getSystemImg());
        lat.setValue(3, 5, 0, -1);
        vid.addLatticeFrame(lat.getSystemImg());
        lat.setValue(3, 6, 0, -1);
        vid.addLatticeFrame(lat.getSystemImg());
        lat.setValue(4, 4, 0, -1);
        vid.addLatticeFrame(lat.getSystemImg());
        //vid.writeVideo();
    }
}