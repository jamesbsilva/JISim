package Backbone.Visualization;

/**
 * 
 *     @(#)   HoneycombLattice2D
 */  

import Backbone.System.LatticeMagInt;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.util.ArrayList;

/** 
 *   HoneycombLattice2D is a visualization of a 2d triangular lattice
 *  <br>
 * 
 *  @param simSys - the simSys for step
 * 
 * 
 * @author      James Silva <jbsilva @ bu.edu>                 
 * @since       2013-05    
 */
public class HoneycombLattice2D implements SysVisualizationInt2D{
    private BufferedImage visImg;
    private int ImageTextOffset = 90;
    private Graphics2D g;
    private int scale;
    private int L;
    private int minImgSize = 600;
    private boolean staggeredSpinDraw = false;
    private boolean useHeter = false;
    private double jInteraction = 0;
    private int R;
    private double cirR = 1.39;
    
    public HoneycombLattice2D(int len,int rin, double jin,boolean useHet, boolean stag){
        L = len;
        R = rin;
        cirR *= Math.pow(L/20,0.25);
        jInteraction = jin;
        staggeredSpinDraw = stag;
        useHeter = useHet;
        System.err.println("AtomicLatticeSumSpin | For Triangular Lattice parameters must be");
        System.err.println("AtomicLatticeSumSpin | 2d with L divisible by 4 for periodic conditions to work");
    }
    
    public BufferedImage getImageOfVis(){return visImg;}
    
    /**
    * initializeImg creates the image layer of the lattice
    *
    */
    public void initializeImg() {
        // Standard scale of images
        scale = 2;
        int sizeImage = scale  * L;
        while (sizeImage < (minImgSize/3)) {
            scale = scale * 2;
            sizeImage = scale * L;
        }
        ImageTextOffset = (int) (sizeImage/4);
        visImg = new BufferedImage(sizeImage*3+(int)(1.5*(cirR-1)*scale),(int)(sizeImage*Math.sqrt(3))+
                ImageTextOffset+(int)((cirR-1)*scale), BufferedImage.TYPE_INT_RGB);
        g = visImg.createGraphics();
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
    public synchronized void updateImg(int i, int j, int spin, int t) {
        //erase old text
        g.setColor(Color.DARK_GRAY);
        g.fillRect(0, 0, visImg.getWidth(), ImageTextOffset);
        g.setColor(Color.WHITE);
        if (t != 0) {
            Font font = new Font("Courier New", Font.BOLD, 72);
            g.setFont(font);
        }
        // 70 is the x coordinate of text
        g.drawString("t = " + t, ((int) L * scale * 6 / 10), ((int) ImageTextOffset * 7 / 10));
        // Draw the update
        if(staggeredSpinDraw && jInteraction < 0){
            if(R <= 1  && (((i)%2) == 0)){
                spin = -1*spin;
            }            
        }
        
        // get color for spin
        g.setColor(getColor(spin));

        drawSpin(i,j);
    }

    /**
    * spinDraw completely redraws the image layer based on the current
    * configuration of the lattice.
    */
    public void spinDraw(LatticeMagInt lattice) {
        int spin;
        g.setColor(Color.WHITE);
        g.fillRect(0, 0, visImg.getWidth(), visImg.getHeight());
        //erase old text
        g.setColor(Color.DARK_GRAY);
        g.fillRect(0, 0, visImg.getWidth(), ImageTextOffset);
        for (int i = 0; i < L; i++) {
            for (int j = 0; j < L; j++) {
                if (useHeter && lattice.isThisFixed(i, j, 0)) {
                    // Paint different color for fixed spins
                    // Cyan is like a blueish
                    // Magenta is purpleish
                    spin = lattice.getValue(i, j);
                    if(staggeredSpinDraw && (jInteraction < 0)){
                        if(R <= 1  && (((i)%2) == 0)){
                            spin = -1*spin;
                        }
                    }
                    // get color for spin
                    g.setColor(getColorFixed(spin));
                } else {
                    spin = lattice.getValue(i, j);
                    if(staggeredSpinDraw && (jInteraction < 0)){
                        if(R <= 1  && (((i)%2) == 0)){
                            spin = -1*spin;
                        }
                    }
                    
                    // get color for spin
                    g.setColor(getColor(spin));
    
                }
                drawSpin(i,j);
            }
        }
    }

    public void fixAspin(int i, int j, int spin){
        // get color for spin
        g.setColor(getColorFixed(spin));
        if(staggeredSpinDraw && jInteraction < 0){
            if(R <= 1  && (((i)%2) == 0)){
                spin = -1*spin;
            }            
        }
        drawSpin(i,j);
     
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
                if(staggeredSpinDraw && (jInteraction < 0)){
                    if(R <= 1  && (((i)%2) == 0)){
                        spin = -1*spin;
                    }
                }
                // get color for spin
                g.setColor(getColor(spin));
                drawSpin(i,j);
            }
        }
    }

    /**
    * spinDraw completely redraws the image layer based on the current
    * configuration of the lattice.
    */
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
                if(staggeredSpinDraw && (jInteraction < 0)){
                    if(R <= 1  && (((i)%2) == 0)){
                        spin = -1*spin;
                    }
                }
                // get color for spin
                g.setColor(getColor(spin));
                drawSpin(i,j);
            }
        }
    }

    public void showNeighbors(ArrayList<Integer> neighLat,int i, int j){
        int  u = 0; int v = 0;
        g.setColor(Color.CYAN);
        int offset = 0;
        for(int k = 0; k < neighLat.size();k++){
            int ind =neighLat.get(k);
            u = ind%L; v = ((int)((double)ind/(double)L))%L;
            if(ind != (i+j*L)){
                drawSpin(u,v);
            }
        }
        g.setColor(Color.MAGENTA);
        drawSpin(i,j);
        
    }
    
    private Color getColorFixed(int spin){
        if(spin == 1){
            return Color.BLUE;
        }else if (spin == (-1)) {
            return Color.YELLOW;
        }else if (spin == (0)) {
            return Color.RED;
        }else{
            return Color.PINK;
        }
    }
    
    private Color getColor(int spin){
        if(spin == 1){
            return Color.WHITE;
        }else if (spin == (-1)) {
            return Color.BLACK;
        }else if (spin == (0)) {
            return Color.GREEN;
        }else{
            return Color.PINK;
        }
    }

    
    private synchronized void drawSpin(int u, int v){
        g.fillOval((int)(getXcoordHoney(u,v)*scale),((int)(getYcoordHoney(u,v)*scale))+
                ImageTextOffset, (int)(scale*cirR)-1, (int)(scale*cirR)-1);
    }
    
       
    public double getXcoordHoney(int i, int j){
        double scaleX = 2.0;
        double offset = (j%2 == 0) ? Math.floor(i/2)*scaleX : (Math.floor(i/2)+1.5)*scaleX;
        return (((double)i)*scaleX)+offset;
    }
    
    public double getYcoordHoney(int i, int j){
        double scaleY = 2.0;
        return ((double)j)*Math.sqrt(3)*scaleY/2.0;
    }

    @Override
    public void updateImg(int i, int dat, int t) {
        updateImg(i%L, (int)((double)i/(double)L)%L, dat, t);
    }
}