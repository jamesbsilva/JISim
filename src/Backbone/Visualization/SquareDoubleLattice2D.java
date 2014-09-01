package Backbone.Visualization;

/**
 * 
 *     @(#)   SquareLattice2D
 */  

import AnalysisAndVideoBackend.DisplayMakerCV;
import Backbone.Util.DataSaver;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.ArrayList;

/** 
 *   SquareLattice2D is a visualization of a 2d triangular lattice
 *  <br>
 * 
 *  @param simSys - the simSys for step
 * 
 * 
 * @author      James Silva <jbsilva @ bu.edu>                 
 * @since       2013-05    
 */
public class SquareDoubleLattice2D implements SysVisualizationDbl2D{
    private BufferedImage visImg;
    private int ImageTextOffset = 0;
    private Graphics2D g;
    private int scale;
    private int L;
    private int minImgSize = 600;
    private double cirR = 1.0;
    private boolean logScale = false;
    private boolean reScale = false; // scale on basis of current max
    private boolean lastMaxInit = false;
    private double lastMax;
    
    public SquareDoubleLattice2D(int len){
        L=len;
    }
    
    public BufferedImage getImageOfVis(){return visImg;}

    /**
    * initializeImg creates the image layer of the lattice
    *
    */
    public void initializeImg() {
        // Standard scale of images
        scale = 5;
        int sizeImage = scale * L;
        while (sizeImage < minImgSize) {
            scale = scale * 2;
            sizeImage = scale * L;
        }
        // Make the image
        visImg = new BufferedImage(sizeImage, sizeImage + ImageTextOffset, BufferedImage.TYPE_INT_RGB);
        g = visImg.createGraphics();
    }

    @Override
    /**
    * updateImg just updates the image layer instead of redrawing the layer.
    * Change time text of image
    *
    * @param i - i coordinate
    * @param spin - new spin value
    * @param t - new time in image
    * @param max - max value
    *
    */
    public void updateImg(int i, double dat, int t,double max) {
        updateImg(i%L, (int)((double)i/(double)L)%L, dat, t, max);
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
    public synchronized void updateImg(int i, int j, double dat, int t,double max) {
        float col = (float)((dat)/(max));
        g.setColor(colorRange(col));
        drawSpin(i,j);
    }

    public synchronized void updateImg(int i, int j,Color col) {
        g.setColor(col);
        drawSpin(i,j);
    }
    
    public synchronized void updateImg(int i, Color col) {
        g.setColor(col);
        int x = i%L;
        int y = ((int)((double)i/(double)L))%L;
        drawSpin(x,y);
    }
    
    /**
    * spinDraw completely redraws the image layer based on the current
    * configuration of the lattice.
    */
    public void spinDraw(double[] data, double max) {
        double dat;
        g.setColor(Color.LIGHT_GRAY);
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
            max = Math.log(max);
            //find min
            double min =10000000;
            for(int u =0; u < data.length;u++){
                if(Math.log(data[u]) < min){
                    min = Math.log(data[u]);
                }
            }
            for(int u =0; u < data.length;u++){
                data[u] = Math.log(data[u])-min;
            }
            max = max-min;
        }
        g.fillRect(0, 0, visImg.getWidth(), visImg.getHeight());
        for (int i = 0; i < L; i++) {for (int j = 0; j < L; j++) {
                dat = data[i+j*L];
                //System.out.println("Problem: "+b+"   dat: "+dat);
                float col = (float)((dat)/(max));
                g.setColor(colorRange(col));
                drawSpin(i,j);
            }
        }    
    }
    
    /**
    * spinDraw completely redraws the image layer based on the current
    * configuration of the lattice.
    */
    public void spinDraw(ArrayList<Double> data, double max) {
        double dat;
        g.setColor(Color.LIGHT_GRAY);
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
        }
        g.fillRect(0, 0, visImg.getWidth(), visImg.getHeight());
        float r; float gr;float b;
        float third = (1.0f/3.0f);
        float rbnd = third;
        float gbnd = 2.0f*third;
        for (int i = 0; i < L; i++) {for (int j = 0; j < L; j++) {
                dat = data.get(i+j*L);
                //System.out.println("Problem: "+b+"   dat: "+dat);
                float col = (float)((dat)/(max));
                g.setColor(colorRange(col));
                drawSpin(i,j);
            }
        }
    }


    /**
    * spinDraw completely redraws the image layer based on the current
    * configuration of the lattice.
    */
    public void spinDraw(float[] data, float max) {
        float dat;
        g.setColor(Color.LIGHT_GRAY);
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
            max = (float)Math.log(max);
            //find min
            float min =10000000;
            for(int u =0; u < data.length;u++){
                if(Math.log(data[u]) < min){
                    min = (float) Math.log(data[u]);
                }
            }
            for(int u =0; u < data.length;u++){
                data[u] = (float) Math.log(data[u])-min;
            }
            max = max-min;
        }        
        g.fillRect(0, 0, visImg.getWidth(), visImg.getHeight());
        int ind;
        for (int i = 0; i < L; i++) {for (int j = 0; j < L; j++) {
                ind = i+j*L;   
                dat = data[ind];
                float col = (float)((dat)/(max));
                g.setColor(colorRange(col));
                drawSpin(i,j);
            }
        }
    
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
    
    private synchronized void drawSpin(int i, int j){
        g.fillOval(i * scale,j * scale + ImageTextOffset, (int)(scale*cirR)-1, (int)(scale*cirR)-1);
    }
    
    // test the class
    public static void main(String[] args) {
        String dirAppnd = "Data/TestSFCL/Geo-6L80R7";
        String dir = "/home/j2/JISim/"+dirAppnd+"/";
        String dir2 = "/home/j2/JISim/";
        String filename = "sfAvg-geo-6-t-21499-seed--160176745-h-0.0";
        String fname = dir+filename;
        System.out.println(dir);
        File path = new File(dir);
        File[] sfFiles = path.listFiles();
        for(int k = 0; k < sfFiles.length;k++){
            if(!(sfFiles[k].getAbsolutePath().contains("sfAvg")) || sfFiles[k].getAbsolutePath().contains("png") ){continue;}
            System.out.println(sfFiles[k].getAbsolutePath());
         
            filename = sfFiles[k].getName();
            dirAppnd = sfFiles[k].getAbsolutePath().substring(dir2.length(),sfFiles[k].getAbsolutePath().length()-filename.length()-1);
            fname = dir+filename;
            
            DataSaver dSave = new DataSaver();
            ArrayList<Double> dat = dSave.getDoubleDataFromFileSingleLine(fname);
            int dim = 2;
            int L = (int)Math.pow(dat.size(),1/((double)dim));
            SquareDoubleLattice2D lat = new SquareDoubleLattice2D((int)(Math.sqrt(dat.size())));
            double max = 0;
            double min = 1000000000;
            System.out.println("L: "+ L);
            double[] data = new double[dat.size()];
            double val;
            for(int u = 0; u < dat.size();u++){
                //val = Math.log(dat.get(u));
                val = dat.get(u);
                if(val > max){max = val;}
                if(val < min){min = val;}
                data[u] = val;
            }
            System.out.println("max: "+max+"    min: "+min);
            for(int u = 0; u < data.length;u++){
                data[u] = data[u]-min;
            }
            max = max - min;

            lat.initializeImg();
            lat.spinDraw(data,max);
            dSave.saveImage(lat.getImageOfVis(),"png", filename+"-SF",dirAppnd);
            DisplayMakerCV dispSF = new DisplayMakerCV("Data");
            dispSF.addDisplayFrame(lat.getImageOfVis());
        }
        System.out.println("Done!");
            
    }
}