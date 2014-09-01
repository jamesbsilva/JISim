package Backbone.Visualization;

/**
* 
*     @(#)   BarsDouble1D
*/  

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.util.ArrayList;

/** 
*   BarsDouble1D is a visualization of a 1d collection of doubles such in bars
*  <br>
* 
*  @param simSys - the simSys for step
* 
* 
* @author      James Silva <jbsilva @ bu.edu>                 
* @since       2013-07    
*/
public class BarsDouble1D implements SysVisualizationDbl{
    private BufferedImage visImg;
    private int ImageTextOffset = 0;
    private int currentTime = 0;
    private int N = 0;
    private Graphics2D g;
    private int delVal;
    private int heightDiv = 850;
    private ArrayList<Double> plotData;
    private int minImgSize = 1200;
    private int minCol = 0;
    private int maxRange = 1150; // Max possible size is 2300 before too small bars
    private int maxCol = 2*maxRange;
    private boolean logScale = false;
    private boolean reScale = false; // scale on basis of current max
    private boolean init = false;
    private boolean lastMaxInit = false;
    private double lastMax;
    private double widthPer = 0.99;
    
    public BarsDouble1D(int nen){
        N = nen;
        if(nen < maxRange){
            maxCol = nen;
        }
        plotData = new ArrayList<Double>();
    }
    
    /**
    *   getImageOfVis returns image of visualization
    */
    public BufferedImage getImageOfVis(){return visImg;}

    /**
    * initializeImg creates the image layer of the lattice
    */
    public void initializeImg() {
        // Standard scale of images
        int sizeImage = minImgSize;
        ImageTextOffset = (int)(heightDiv*0.1);
        // Make the image
        visImg = new BufferedImage(sizeImage, heightDiv + ImageTextOffset, BufferedImage.TYPE_INT_RGB);
        g = visImg.createGraphics();
        init = true;
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
    public synchronized void updateImg(int i,double dat, int t,double max) {
        if(max > lastMax){
            rescaleBars(max);
        }
        if(t != 0){
            //erase old text
            g.setColor(Color.DARK_GRAY);
            g.fillRect((int)( (double)visImg.getWidth() * 1.0 / 10.0), 0, 
                    (int)( (double)visImg.getWidth()), ImageTextOffset);

            g.setColor(Color.WHITE);
            if (t != 0) {
                Font font = new Font("Courier New", Font.BOLD, 32);
                g.setFont(font);
            }
            // 70 is the x coordinate of text
            g.drawString("t = " + t+"     range("+minCol+","+maxCol+")", (int)( 
                    (double)visImg.getWidth() * 1.0 / 10.0), 
                    (int)( ImageTextOffset * 8.0 / 10.0));
        }
        if(i > maxCol){
            rescaleBasedOnNCol(i);
        }else if(i < minCol){
            rescaleUnderMin(i);
        }
        drawBar(i,dat);
    }

    /**
    * spinDraw completely redraws the image layer based on the current
    * configuration of the lattice.
    */
    public void spinDraw(double[] data, double max) {
        double dat;
        N = data.length;
        g.setColor(Color.DARK_GRAY);
        g.fillRect(0, 0, visImg.getWidth(), visImg.getHeight());
        if(!reScale){
            if(!lastMaxInit){
                lastMax = max;
            }else{
                if(lastMax < max){
                    lastMax = max;
                }else{
                    max = lastMax;
                }
            }
        }
        if(logScale){
            max = logMin(data, max);
        }

        for (int i = 0; i < data.length; i++) {
            addData(i,data[i]);
        }
        int maxColDraw = (plotData.size() < maxCol) ? plotData.size() : maxCol ;
        for(int j = minCol; j < maxColDraw; j++) {
            drawBar(j,data[j]);
        }
    }
    
    /**
    * spinDraw completely redraws the image layer based on the current
    * configuration of the lattice.
    */
    public void spinDraw(ArrayList<Double> data, double max) {
        N = data.size();
        double dat;
        g.setColor(Color.DARK_GRAY);
        g.fillRect(0, 0, visImg.getWidth(), visImg.getHeight());
        if(!reScale){
            if(!lastMaxInit){
                lastMax = max;
            }else{
                if(lastMax < max){
                    lastMax = max;
                }else{
                    max = lastMax;
                }
            }
        }
        if(logScale){
            max = logMin(data, max);
        }
        for (int i = 0; i < data.size(); i++) {
            addData(i,data.get(i));
        }
        int maxColDraw = (plotData.size() < maxCol) ? plotData.size() : maxCol ;
        for(int j = minCol; j < maxColDraw; j++) {
            drawBar(j,data.get(j));
        }
    }

    /**
    * spinDraw completely redraws the image layer based on the current
    * configuration of the lattice.
    */
    public void spinDraw(float[] data, float max) {
        N = data.length;
        float dat;
        g.setColor(Color.DARK_GRAY);
        g.fillRect(0, 0, visImg.getWidth(), visImg.getHeight());
        if(!reScale){
            if(!lastMaxInit){
                lastMax = max;
            }else{
                if(lastMax < max){
                    lastMax = max;
                }else{
                    max = (float)lastMax;
                }
            }
        }
        if(logScale){
            max = logMin(data,max);
        }        
        for (int i = 0; i < data.length; i++) {
            addData(i,data[i]);
        }
        int maxColDraw = (plotData.size() < maxCol) ? plotData.size() : maxCol ;
        for(int j = minCol; j < maxColDraw; j++) {
            drawBar(j,data[j]);
        }
    }
    
    private double logMin(ArrayList<Double> data, double max){
        max = Math.log(max);
        //find min
        double min =10000000;
        for(int u =0; u < data.size();u++){
            if( Math.log(data.get(u)) < min){
                min = Math.log(data.get(u));
            }
        }
        for(int u =0; u < data.size();u++){
            data.set(u, Math.log(data.get(u))-min);
        }
        max = max-min;
        return max;
    }
    
    private float logMin(float[] data, float max){
        max = (float)Math.log(max);
        //find min
        double min =10000000;
        for(int u =0; u < data.length;u++){
            if( Math.log(data[u]) < min){
                min = Math.log(data[u]);
            }
        }
        for(int u =0; u < data.length;u++){
            data[u] = (float) (Math.log(data[u])-min);
        }
        max = (float) (max-min);
        return max;
    }
    
    private double logMin(double[] data, double max){
        max = Math.log(max);
        //find min
        double min =10000000;
        for(int u =0; u < data.length;u++){
            if( Math.log(data[u]) < min){
                min = Math.log(data[u]);
            }
        }
        for(int u =0; u < data.length;u++){
            data[u] = (Math.log(data[u])-min);
        }
        max = max-min;
        return max;
    }
    
    private void addData(int ind, double val){
        if(ind < plotData.size()){
            plotData.set(ind, val);
        }else{
            //System.out.println("Data  | "+ind);
            int sizeOld = plotData.size();
            for(int u = 0; u < (ind - sizeOld);u++){
                if((sizeOld+u) == ind ){
                    plotData.add(val);
                }else{
                    plotData.add(0.0);
                }
            }
        }
    }
    
    private synchronized void drawBar(int i, double val){
        // Do not draw if not in drawing range
        if(i > maxCol || i < minCol){
            return;
        }
        // draw 
        double n = (maxCol-minCol);
        double del = ((double)visImg.getWidth()/n);
        double width = (del*widthPer);
        int barAreaHeight =  visImg.getHeight()-ImageTextOffset-5 ;
        double bar = (double)(barAreaHeight)*val/lastMax;
        double ind = (i-minCol)*del;
        if(width < ((double)visImg.getWidth()*widthPer/(2*maxRange+1))){
            rescaleBasedOnColWidth(width);
        }
        //System.out.println("Drawing Bar in Pixel | "+(int)Math.floor(ind)+
        //        "    del | "+visImg.getHeight()+"     ind| "+ind+"    barSize | "+bar+"    val | "+val+"   max| "+lastMax);
        g.setColor(Color.WHITE);
        g.fillRect((int)Math.floor(ind) ,ImageTextOffset , (int)del, barAreaHeight);
        g.setColor(Color.BLUE);
        g.fillRect((int)Math.floor(ind) ,visImg.getHeight()-(int)bar,  (int)width, (int)bar );
    }
    
    
    private synchronized void rescaleUnderMin(int newVal){
        maxCol = newVal+(int)((double)maxRange/2.0);
        minCol = newVal-(int)((double)maxRange/2.0); 
        redrawBars();
    }
    
    private synchronized void rescaleBasedOnNCol(int newVal){
        if( newVal < maxCol && minCol == 0){
            double del = ((double)visImg.getWidth()/newVal);
            double width = (int)(del*widthPer);
            rescaleBasedOnColWidth(width);
        }else{
            maxCol = newVal+(int)((double)maxRange/2.0);
            minCol = newVal-(int)((double)maxRange/2.0);
            redrawBars();
        }
    }
    
    private synchronized void rescaleBasedOnColWidth(double newVal){
        int w = (int)((double)visImg.getWidth()/newVal);
        System.out.println("Current ColumnWidth | "+newVal+"    Width| "+visImg.getWidth()+"    w | "+w
                +"     N | "+N);
        visImg = new BufferedImage(w, visImg.getHeight(), BufferedImage.TYPE_INT_RGB);
        visImg.createGraphics();
        g = (Graphics2D) visImg.getGraphics();
        spinDraw(plotData,lastMax);
    }
    
    private synchronized void redrawBars(){
        N = plotData.size();
        //System.out.println("Redraw Bars | "+minCol+"     to | "+maxCol);
        int maxDrawCol = ( plotData.size() < maxCol ) ? plotData.size() : maxCol;
        for(int u = minCol; u < maxDrawCol; u++){
            updateImg(u,plotData.get(u), currentTime,lastMax);
        }        
    }
    
    private synchronized void rescaleBars(double newMax){
        if(newMax <= lastMax){return;}
        lastMax = newMax;
        N = plotData.size();
        //System.out.println("Rescale Bars | "+minCol+"     to | "+maxCol);
        int maxDrawCol = ( plotData.size() < maxCol ) ? plotData.size() : maxCol;
        for(int u = minCol; u < maxDrawCol; u++){
            updateImg(u,plotData.get(u), currentTime,lastMax);
        }        
    }
    
    // test the class
    public static void main(String[] args) {
            
    }
}