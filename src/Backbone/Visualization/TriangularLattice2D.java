package Backbone.Visualization;

/**
 * 
 *     @(#)   TriangularLattice2D
 */  

import Backbone.System.LatticeMagInt;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.util.ArrayList;

        
/** 
 *   TriangularLattice2D is a visualization of a 2d triangular lattice
 *  <br>
 * 
 *  @param simSys - the simSys for step
 * 
 * 
 * @author      James Silva <jbsilva @ bu.edu>                 
 * @since       2013-05    
 */
public class TriangularLattice2D implements SysVisualizationInt2D{

    private BufferedImage visImg;
    private int ImageTextOffset = 90;
    private Graphics2D g;
    private int scale;
    private int L;
    private int minImgSize = 400;
    private boolean staggeredSpinDraw = false;
    private boolean useHeter = true;
    private double jInteraction = 0;
    private int R;
    private double cirR = 1.7;
    
    public TriangularLattice2D(int len,int rin, double jin,boolean useHet, boolean stag){
        L = len;
        R = rin;
        jInteraction = jin;
        staggeredSpinDraw = stag;
        useHeter = useHet;
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
        while (sizeImage < (minImgSize/2)) {
            scale = scale * 2;
            sizeImage = scale * L;
        }
        // Make the image
        ImageTextOffset = (int) (sizeImage/4);
        visImg = new BufferedImage(sizeImage*2+(int)((cirR-1)*scale),(int)(sizeImage*Math.sqrt(3))+
                ImageTextOffset+(int)((cirR-1)*scale), BufferedImage.TYPE_INT_RGB);
        g =  visImg.createGraphics();
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
        g.fillRect(0, 0, scale*2*L+(int)((cirR-1)*scale), ImageTextOffset);
        g.setColor(Color.WHITE);
        if (t != 0) {
            Font font = new Font("Courier New", Font.BOLD, 72);
            g.setFont(font);
        }
        // 70 is the x coordinate of text
        g.drawString("t = " + t, ((int) L * scale / 10), ((int) ImageTextOffset * 8 / 10));
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
    public void spinDraw(LatticeMagInt lattice) {
        int spin;
        g.setColor(Color.WHITE);
        g.fillRect(0, 0, visImg.getWidth(), visImg.getHeight());
        //erase old text
        g.setColor(Color.DARK_GRAY);
        g.fillRect(0, 0, scale*2*L+(int)((cirR-1)*scale), ImageTextOffset);
        for (int i = 0; i < L; i++) {
            for (int j = 0; j < L; j++) {
                if (useHeter && lattice.isThisFixed(i, j, 0)) {
                    // Paint different color for fixed spins
                    // Cyan is like a blueish
                    // Magenta is purpleish
                    spin = lattice.getValue(i, j);
                    if(staggeredSpinDraw && (jInteraction < 0)){
                        if(R < 1  && (((i+j)%2) == 0)){
                            spin = -1*spin;
                        }
                    }                    
                    // get color for fixed spin
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
    
    public void showNeighbors(ArrayList<Integer> neighLat,int i, int j){
        int  u = 0; int v = 0;
        System.out.println("TriangularLattice2D | Showing Neighbors");
        g.setColor(Color.CYAN);
        for(int k = 0; k < neighLat.size();k++){
            int ind = neighLat.get(k);
            u = ind%L; v = ((int)((double)ind/(double)L))%L;            
            if((u+v*L) != (i+j*L)){
                drawSpin(u,v);
            } 
        }
        g.setColor(Color.MAGENTA);
        g.fillOval((int)(getXhex(i,j)*scale),((int)(getYhex(i,j)*scale))+
                ImageTextOffset+1, (int)(scale*cirR)-1, (int)(scale*cirR)-1); 
    }
}    
