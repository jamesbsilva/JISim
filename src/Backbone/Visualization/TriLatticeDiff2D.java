package Backbone.Visualization;

/**
 * 
 *     @(#)   TriLatticeDiff2D
 */  

import Backbone.Util.ConfigDifference;
import com.googlecode.javacv.CanvasFrame;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.util.ArrayList;

        
/** 
 *   TriLatticeDiff2D is a visualization of the difference between configurations on a 2d triangular lattice
 *  <br>
 * 
 *  @param simSys - the simSys for step
 * 
 * 
 * @author      James Silva <jbsilva @ bu.edu>                 
 * @since       2013-05    
 */
public class TriLatticeDiff2D  {
    private BufferedImage visImg;
    private int ImageTextOffset = 90;
    private Graphics2D g;
    private int scale;
    private int L;
    private int minImgSize = 400;
    private double cirR = 1.7;
    private double max = 0;
    private double min = 0;
    private boolean drawAvg = true;
    private boolean binaryMode = false;
    private boolean logColor = false;
    private boolean expColor = true;
    private CanvasFrame canvasFrame ;
    
    public TriLatticeDiff2D(int len){
        L = len;
    }

    public BufferedImage getImageOfVis(){return visImg;}
    
    public void setBinaryDrawOffON(boolean bl){
        binaryMode = bl;
    }
    
    public void setMaxMinVal(double ma , double mi){
        max = ma;
        min = mi;
    }

    /**
    * initializeImg creates the image layer of the lattice
    *
    */
    public void initializeImg() {
        // Standard scale of images
        scale = 5;
        int sizeImage = scale * L;
        while (sizeImage < (minImgSize/2)) {
            scale = scale * 2;
            sizeImage = scale * L;
        }
        // Make the image
        ImageTextOffset = (int) (sizeImage/4);
        visImg = new BufferedImage(sizeImage*2+(int)((cirR-1)*scale),(int)(sizeImage*Math.sqrt(3))+
                ImageTextOffset+(int)((cirR-1)*scale), BufferedImage.TYPE_INT_RGB);
        g =  visImg.createGraphics();
    
        // intialize
        g.setColor(Color.WHITE);
        g.fillRect(0, 0, visImg.getWidth(), visImg.getHeight());
        //erase old text
        g.setColor(Color.DARK_GRAY);
        g.fillRect(0, 0, scale*2*L+(int)((cirR-1)*scale), ImageTextOffset);

    }

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
    public synchronized void updateImg(int i, int j, double spin, int t) {
        //erase old text
        g.setColor(Color.DARK_GRAY);
        g.fillRect(0, 0, scale*2*L+(int)((cirR-1)*scale), ImageTextOffset);
        g.setColor(Color.WHITE);
        if (t != 0) {
            Font font = new Font("Courier New", Font.BOLD, 72);
            g.setFont(font);
        }
        // 70 is the x coordinate of text
        g.drawString("t = " + t, ((int) L * scale / 10), ((int) ImageTextOffset * 8 / 10));

        // get color 
        if(binaryMode){
            g.setColor(getColorBinary(spin));
        }else{
            g.setColor(getColorRange(spin));
        }
        drawSpin(i,j);
    }

    /**
    * spinDraw completely redraws the image layer based on the current
    * configuration of the lattice.
    */
    public void spinDraw(ConfigDifference diff) {
        int spin;
        ArrayList<Double> spinDiff;
        double nMeasure = diff.getConfigN();
        if(drawAvg){
            max = diff.getMaxAvgSum();
            min = diff.getMinAvgSum();
            spinDiff = diff.getConfigDiffAvgSum();
            max /= nMeasure;
            min /= nMeasure;
            if(logColor){
                min = Math.log(min);
                max = Math.log(max);
            }else if(expColor){
                min = Math.exp(min);
                max = Math.exp(max);
            }
            //max = 1;min = 0;
            //max /= nMeasure;
            //min /= nMeasure;
        }else{        
            max = diff.getMaxDiff();
            min = diff.getMinDiff();
            spinDiff = diff.getConfigDiff();
        }
        
        //erase old text
        g.setColor(Color.DARK_GRAY);
        g.fillRect(0, 0, scale*2*L+(int)((cirR-1)*scale), ImageTextOffset);
        g.setColor(Color.WHITE);
        double currVal;
        for (int i = 0; i < L; i++) {
            for (int j = 0; j < L; j++) {
                // get color 
                currVal = spinDiff.get(i+j*L);
                //System.out.println(" sumVal | "+currVal+"      configs| "+nMeasure+"     max| "+max+"    min| "+min);
                if(drawAvg){currVal /= nMeasure;}
                if(binaryMode){
                    g.setColor(getColorBinary(currVal));
                }else{
                    g.setColor(getColorRange(currVal));
                }
                drawSpin(i,j);
            }
        }
    }
    
        
    /**
    * spinDraw completely redraws the image layer based on the current
    * configuration of the lattice.
    */
    public void spinDraw(int[] lattice) {
        int spin;
        g.setColor(Color.WHITE);
        g.fillRect(0, 0, visImg.getWidth(), visImg.getHeight());
        //erase old text
        g.setColor(Color.DARK_GRAY);
        g.fillRect(0, 0, visImg.getWidth(), ImageTextOffset);
        for (int i = 0; i < L; i++) {
            for (int j = 0; j < L; j++) {
                    spin = lattice[i+ j*L];
                    if(spin == 1){
                        g.setColor(Color.ORANGE);
                    }
                    if (spin == (-1)) {
                        g.setColor(Color.BLACK);
                    }
                    if (spin == (0)) {
                        g.setColor(Color.YELLOW);
                    }
                drawSpin(i,j);
            }
        }
    }
    
    public void spinDraw(ArrayList<Integer> lattice) {
        int spin;
        g.setColor(Color.WHITE);
        g.fillRect(0, 0, visImg.getWidth(), visImg.getHeight());
        //erase old text
        g.setColor(Color.DARK_GRAY);
        g.fillRect(0, 0, visImg.getWidth(), ImageTextOffset);
        for (int i = 0; i < L; i++) {
            for (int j = 0; j < L; j++) {
                    spin = lattice.get(i+ j*L);
                    if(spin == 1){
                        g.setColor(Color.ORANGE);
                    }
                    if (spin == (-1)) {
                        g.setColor(Color.BLACK);
                    }
                    if (spin == (0)) {
                        g.setColor(Color.YELLOW);
                    }
                drawSpin(i,j);
            }
        }
    }
    
    private Color getColorBinary(double spin){
        if(spin == 0){
            return Color.WHITE;
        }else{
            return Color.BLACK;
        }
    }

    private Color getColorRange(double spin){        
        spin = scaleColor(spin);
        //System.out.println("spin | "+spin+"     max| "+max+"    min| "+min);
        float col = (float)((spin-min)/(max-min));
        //System.out.println("Color in | "+col+"      spin | "+spin+"     max| "+max+"    min| "+min);
        return colorRange(col);
    }
    
    private Color colorRange(float col){
        float seg = 7.0f;
        int colScale = (int)(Math.floor(col*seg));
        col  = col*seg-(float)colScale;     
        switch (colScale) {
            // 0 1 1
            case 0: return new Color(1.0f-col,   1.0f,   1.0f);
            // 0 0 1
            case 1: return new Color(0.0f,  1.0f-col,   1.0f);
            // 1 0 1
            case 2: return new Color(col,  0.0f,  1.0f);
            // 1 0 0
            case 3: return new Color(1.0f,  0.0f,  1.0f-col);
            // 1 1 0
            case 4: return new Color(1.0f,   col,   0.0f);
            // 0 1 0
            case 5: return new Color(1.0f-col,  1.0f,  0.0f);
            // 0 0 0 
            default: return new Color(0.0f,   1.0f-col,   0.0f);      
        }        
    }
    
    private double scaleColor(double spin){
        if(logColor){
            //System.out.println("spin | "+spin+"     max| "+max+"    min| "+min);
            spin = Math.log(spin);
        }else if(expColor){
            spin = Math.exp(spin);
        }
        return spin;
    }
    
    private Color getColor(int spin){
        if(spin == 1){
            return Color.ORANGE;
        }else if (spin == (-1)) {
            return Color.BLACK;
        }else if (spin == (0)) {
            return Color.YELLOW;
        }else{
            return Color.WHITE;
        }
    }
    
    // Visualization increases spacing in lattice by a factor of to to have 
    // y coordinate spaced out by factor of sqrt(3) instead of sqrt(3)/2 
    private double getXhex(int i, int j){
        return (j%2 == 0) ? (double)2*i : (double)2*i + 1;
    }
    
    // Visualization increases spacing in lattice by a factor of to to have 
    // y coordinate spaced out by factor of sqrt(3) instead of sqrt(3)/2
    private double getYhex(int i, int j){
        return ((double)j)*Math.sqrt(3);
    }
    
    private synchronized void drawSpin(int u, int v){
        g.fillOval((int)(getXhex(u,v)*scale),((int)(getYhex(u,v)*scale))+
                ImageTextOffset+1, (int)(scale*cirR)-1, (int)(scale*cirR)-1);
    }
    
    public void initShowImage(){
        canvasFrame = new CanvasFrame("Lat");
        canvasFrame.setCanvasSize(getImageOfVis().getWidth(),getImageOfVis().getHeight());
    }
    public void updateImage(){
        canvasFrame.showImage(getImageOfVis());
    }
    
    public void endShowImage(){
        canvasFrame.dispose();
    }
    
    public void showColorGradient(){
        ArrayList<Double> grad = new ArrayList<Double>();
        max = 1.0;
        min = 0;
        //erase old text
        g.setColor(Color.DARK_GRAY);
        g.fillRect(0, 0, scale*2*L+(int)((cirR-1)*scale), ImageTextOffset);
        g.setColor(Color.WHITE);
        double currVal;
        if(logColor){
            currVal = ((1.0f/((double)(L*L)))-0.00001f);
            min = Math.log(currVal);
            max = Math.log(max);
        }else if(expColor){
            currVal = ((1.0f/((double)(L*L)))-0.00001f);
            min = Math.exp(currVal);
            max = Math.exp(max);
        }
        for (int i = 0; i < L; i++) {
            for (int j = 0; j < L; j++) {
                currVal = (float)(i+j*L)/(float)(L*L);
                if(currVal == 0){currVal = ((1.0f/((double)(L*L)))-0.00001f);}
                System.out.println("Val | "+currVal);
                g.setColor(getColorRange(currVal));
                drawSpin(i,j);
            }
        }    
    }
    
    // test the class
    public static void main(String[] args) {
        TriLatticeDiff2D vis = new TriLatticeDiff2D(50);
        vis.initializeImg();
        vis.initShowImage();
        vis.showColorGradient();
        vis.updateImage();
        
    }
}    
