package Backbone.Visualization;

/**
 *     @(#)   SquareLattice2D
 */  

import Backbone.System.LatticeInt;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
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
public class SquareLattice2D implements SysVisualizationInt2D{
    private BufferedImage visImg;
    private int minImgSize = 600; private int ImageTextOffset = 90;
    private Graphics2D g; private int L; private int R;
    private double jInteraction = 0;
    private int scale;
    private boolean staggeredSpinDraw = false;
    private boolean useHeter = true;
    private double cirR = 1.0;
    
    public SquareLattice2D(int len,int rin, double jin,boolean useHet, boolean stag){
        L = len; R = rin; jInteraction = jin;
        staggeredSpinDraw = stag;
        useHeter = useHet;
    }
    
    public BufferedImage getImageOfVis(){return visImg;}

    /**
    * initializeImg creates the image layer of the lattice
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
        visImg.createGraphics();
        g = (Graphics2D) visImg.getGraphics();
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
    public synchronized void updateImg(int i, int spin, int t) {
        updateImg(i%L,(int)((double)i/(double)L)%L,spin,t);
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
        g.fillRect(((int) L * scale * 3 / 10), 0, ((int) L * scale * 6 / 10), ImageTextOffset);
        
        g.setColor(Color.WHITE);
        if (t != 0) {
            Font font = new Font("Courier New", Font.BOLD, 72);
            g.setFont(font);
        }
        // 70 is the x coordinate of text
        g.drawString("t = " + t, ((int) L * scale * 3 / 10), ((int) ImageTextOffset * 6 / 10));
        // Draw the update
        if(staggeredSpinDraw && jInteraction < 0){
            if(R < 1  && (((i+j)%2) == 0)){
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
    public void spinDraw(LatticeInt lattice) {
        int spin;
        g.setColor(Color.LIGHT_GRAY);        
        g.fillRect(0, 0, visImg.getWidth(), visImg.getHeight());
        g.setColor(Color.DARK_GRAY);
        g.fillRect(0, 0, L * scale, ImageTextOffset);
        
        for (int i = 0; i < L; i++) {
            for (int j = 0; j < L; j++) {
                if (useHeter && lattice.isThisFixed(i, j, 0)){
                    // Paint different color for fixed spins
                    spin = lattice.getValue(i, j);
                    if(staggeredSpinDraw && (jInteraction < 0)){
                        if(R < 1  && (((i+j)%2) == 0)){
                            spin = -1*spin;
                        }
                    }
                    // get color for spin
                    g.setColor(getColorFixed(spin));
                
                } else {
                    spin = lattice.getValue(i, j);
                    if(staggeredSpinDraw && (jInteraction < 0)){
                        if(R < 1  && (((i+j)%2) == 0)){
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
    
    public Color getColorFixed(int spin){
        if(spin == 1){
            return Color.BLUE;
        }else if (spin == (-1)) {
            return Color.MAGENTA;
        }else if (spin == (0)) {
            return Color.RED;
        }else{
            return Color.PINK;
        }
    }
    
    public Color getColor(int spin){
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

    
    private synchronized void drawSpin(int i, int j){
        g.fillOval(i * scale,j * scale + ImageTextOffset, (int)(scale*cirR)-1, (int)(scale*cirR)-1);
    }
    
    public void showNeighbors(ArrayList<Integer> neighLat,int i, int j){
        int  u = 0; int v = 0;
        System.out.println("SquareLattice2D | Showing Neighbors");
        g.setColor(Color.CYAN);
        for(int k = 0; k < neighLat.size();k++){
            int ind = neighLat.get(k);
            u = ind%L; v = ((int)((double)ind/(double)L))%L;            
            if((u+v*L) != (i+j*L)){
                drawSpin(u,v);
            } 
        }
        g.setColor(Color.MAGENTA);
        drawSpin(i,j);
    }    
}