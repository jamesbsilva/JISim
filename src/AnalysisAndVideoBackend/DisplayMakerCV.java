package AnalysisAndVideoBackend;

/**
 *   @(#) DisplayMakerCV
 */ 

import java.awt.image.BufferedImage;
import com.googlecode.javacv.CanvasFrame;
import com.googlecode.javacv.cpp.opencv_core.IplImage;

/**  DisplayMakerCV uses OpenCV to display 2D visualization images.
 *  <br>
 *  @param seed - seed for this particular run
 *  @param max - max amount of frames
 *  @param type - image file type that will go in avi
 *  @param Name - optional add in to file name
 * 
 * @author      James Silva <jbsilva @ bu.edu>                 
 * @since       2013-04    
 */
public class DisplayMakerCV{
    private boolean dispInitialized = false;
    private CanvasFrame canvasFrame ;
    private String windowStr = "Simulation Display";
    
    public DisplayMakerCV(String winName){
        windowStr = winName;
    }
 
    public DisplayMakerCV(){
    }
 
    public IplImage convertBuffImageToJCVimg(BufferedImage imgBuff){
        IplImage img = IplImage.createFrom(imgBuff); 
        return img;
    }
    
    private void initCanvas(BufferedImage buffImg ){
            canvasFrame = new CanvasFrame(windowStr);
            canvasFrame.setCanvasSize(buffImg.getWidth(),buffImg.getHeight());
            dispInitialized = true;
    }
    
    /**
    *      writeVideo writes the video in avi format to the video directory.
    */
    public void closeWindow(){
        canvasFrame.dispose();
    }
    
    public void addDisplayFrame(BufferedImage buffImg){
        if(!dispInitialized){
            initCanvas(buffImg);               
        } 
        canvasFrame.showImage(buffImg);    
    }
    
    public static void main(String[] args) {
    }
}