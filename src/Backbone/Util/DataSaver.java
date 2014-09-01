package Backbone.Util;

/**
*    @(#) DataSaver
*/

import Backbone.System.LatticeMagInt;
import Backbone.System.Network;
import Backbone.System.SimSystem;
import Backbone.System.SimpleLattice;
import Backbone.Visualization.SquareLattice2D;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Scanner;
import java.util.concurrent.atomic.AtomicIntegerArray;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.imageio.ImageIO;

/** 
*   A class that deals with saving the data files that are necessary for
*  any simulations.
*  <br>
* 
* 
* @author      James Silva <jbsilva @ bu.edu>                 
* @since       2011-09    
*/
public class DataSaver {
    private ParameterBank param; private DirAndFileStructure dir;
    private int instId=1; // Instance Id for multiple instances
    private int maxInstances=8;
    private int minImageSize = 600;
    private int minScale = 5;
    private int L;
    private boolean systemIsNetwork =false;
    private SquareLattice2D sq = new SquareLattice2D(1,1,0.0,true,false);
        
    public DataSaver(){
        dir = new DirAndFileStructure(); param = new ParameterBank(); L = param.L;
        updateInstances(1);
    }

    /**
    * 
    * @param id - id of instance of this program
    */
    public DataSaver(int id){
        this(id,"");
    }

    /**
    * 
    * @param id - id of instance of this program
    * @param paramPost - postfix of the parameter file
    */
    public DataSaver(int id,String paramPost){
        dir = new DirAndFileStructure();instId=id; param = new ParameterBank(paramPost);L=param.L;
    }

    /**
    * 
    * @param id - id of instance of this program
    * @param param2 - parameters
    */
    public DataSaver(int id,ParameterBank param2){
        dir = new DirAndFileStructure();instId=id; param = param2;L=param.L;
    }
    
    /**
    * 
    * @param paramPost - postfix of the parameter file
    */
    public DataSaver(String paramPost){
        dir = new DirAndFileStructure();updateInstances(1); 
        param = new ParameterBank(paramPost);L=param.L;
    }
    
    /**
    * 
    * @param param2 - parameters
    */
    public DataSaver(ParameterBank param2){
        dir = new DirAndFileStructure();updateInstances(1); param = param2;L=param.L;
    }
    
    /**
    *   setL sets the system size.
    * 
    * @param length - new System Size
    */    
    public void setL(int length){ L = length; }
    
    /**
    *      updateClosedInstance updates the instances file to reflect closing of 
    *  this program
    */
    public void updateClosedInstance(){updateInstances(-1);}
    
    /**
    *      getInstId returns the instance id of this program instance
    * 
    * @return id- instance id of program
    */
    public int getInstId(){ return instId; }

    /**
    *      saveArrayFileAsImg saves a configuration file for the given time and 
    *  does this in the form of an image in png format.
    * 
    * @param fixtype - type of fixed spin configuration
    * @param t - current time
    * @param run - current run
    * @param seed - current seed
    * @param tOff - any offset to the time for this image
    */
    public void saveArrayFileAsImg(String fixtype,int t, int run,int seed,int tOff){
        saveArrayFileAsImg(fixtype,t,run,minScale,"png",seed,tOff,"");
    }
    
    /**
    *      saveArrayFileAsImg saves a configuration file for the given time and 
    *  does this in the form of an image in png format.
    * 
    * @param fixtype - type of fixed spin configuration
    * @param t - current time
    * @param run - current run
    * @param seed - current seed
    * @param tOff - any offset to the time for this image
    * @param sub - any subdirectory of images to go in
    */
    public void saveArrayFileAsImg(String fixtype,int t, int run,int seed,int tOff,String sub){
        saveArrayFileAsImg(fixtype,t,run,minScale,"png",seed,tOff,sub);
    }
    
    /**
    *      saveArrayFileAsImg saves a configuration file for the given time and 
    *  does this in the form of an image in png format.
    * 
    * @param t - current time
    * @param run - current run
    * @param seed - current seed
    * @param tOff - any offset to the time for this image
    */
    public void saveArrayFileAsImg(int t, int run,int seed,int tOff){
        saveArrayFileAsImg("",t,run,minScale,"png",seed,tOff,"");
    }
    
    /**
    *      saveArrayFileAsImg saves a configuration file for the given time and 
    *  does this in the form of an image in png format.
    * 
    * @param t - current time
    * @param run - current run
    * @param seed - current seed
    */
    public void saveArrayFileAsImg(int t, int run,int seed){
        saveArrayFileAsImg(t,run,minScale,"png",seed);
    }
    
    /**
    *      saveArrayFileAsImg saves a configuration file for the given time and 
    *  does this in the form of an image in png format.
    * 
    * @param t - current time
    * @param run - current run
    */
    public void saveArrayFileAsImg(int t, int run){
        saveArrayFileAsImg(t,run,minScale,"png",0);
    }

    /**
    *      updateInstances opens the instances file and changes the value
    *  by the input amount
    * 
    * @param changeInst - change instance file by this amount
    */
    @SuppressWarnings({"CallToThreadDumpStack", "ConvertToTryWithResources"})
    private void updateInstances(int changeInst){
        String fname = dir.getSettingsDirectory()+"instances.txt"; int id=1;
        //read file first    
        try{
            Scanner scanner = new Scanner(new File(fname));
            id = scanner.nextInt();
        }catch (FileNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        //update 
        id = id+changeInst;
        
        // reset if clearly too many instances
        if(id>maxInstances){id=1;}
        try {
            PrintStream out = new PrintStream(new FileOutputStream(fname));
            out.println(id); out.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        instId = id;
    }

    /**
    *       clearInstances opens the instances file and sets the value to zero.
    */
    @SuppressWarnings("ConvertToTryWithResources")
    public void clearInstances(){
        String fname = dir.getSettingsDirectory()+"instances.txt";
        try {
            PrintStream out = new PrintStream(new FileOutputStream(fname));
            out.println(0); out.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    /**
    *       doesConfigExist returns true if the given configuration with given
    *   run and time actually exist.
    * 
    * @param run - run for the configuration file
    * @param t - time of lattice state in the given run
    * @return true if configuration file is present
    */
    public boolean doesConfigExist(int run, int t){
    	boolean fstate; String fname = "";
    	fname = dir.getConfigDirectory(instId)+getLatticeTempFilename(run, t, fname);
    	fstate = (new File(fname)).exists();
    	return fstate;
    }
    
    
    public void saveImage(BufferedImage img, String type,String filename){
        saveImage(img,type,filename,"");
    }
    
    /**
    *       saveImage saves the given image in the image directory
    * 
    * @param img - input image
    * @param type - image format typer
    * @param filename -filename of image
    */
    public void saveImage(BufferedImage img, String type,String filename, String newDir){
        if(filename.equalsIgnoreCase("")){filename = "Config";} String fname;
        if( newDir == "" ){
            fname = dir.getImageDirectory()+filename+"."+type;
        }else{
            fname = dir.getRootDirectory()+newDir+"/"+filename+"."+type;
        }
        createDirIfUnInit(fname);
        File output = new File(fname);
        try {
            ImageIO.write(img, type, output);
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }

    public void saveLatticeAsImg(int t ,int run,int seed,LatticeMagInt lattice){
        saveLatticeAsImg(t,run,seed,lattice,"");
    }
    
    
    /**
    *      saveLatticeAsImg saves the given lattice as an image
    *  in the image directory.
    * 
    * 
    * @param t - time of lattice state
    * @param run - run to save the lattice tim
    * @param seed - seed of given run
    * @param lattice - lattice to save as image
    */
    public void saveLatticeAsImg(int t ,int run,int seed,LatticeMagInt lattice,String newDir){
        String fname;
        if(newDir == ""){
            fname = dir.getImageDirectory()+getLatticeImgFilename(run, t, seed);
        }else{
            fname = dir.getRootDirectory()+newDir+"/"+getLatticeImgFilename(run, t, seed);
        }
        try {
            // retrieve image
            BufferedImage bi = lattice.getSystemImg();
            File outputfile = new File(fname);
            createDirIfUnInit(fname);
            ImageIO.write(bi, "png", outputfile);
        } catch (IOException e) {
        }
    }
    
    /**
    *      saveArrayFileAsImg saves a configuration file for the given time and 
    *  does this in the form of an image in png format.
    * 
    * @param t - current time
    * @param run - current run
    * @param scale - scale of image
    * @param type - format type of image
    * @param seed - current seed
    */
    public void saveArrayFileAsImg(int t, int run,int scale, String type,int seed){
        saveArrayFileAsImg("",t,run,scale,type,seed,0,"");
    }
    
    /**
    *      saveArrayFileAsImg saves a configuration file for the given time and 
    *  does this in the form of an image in png format.
    * 
    * @param fixType - type of fixed spin configuration
    * @param t - current time
    * @param run - current run
    * @param scale - scale of image
    * @param type - format type of image
    * @param seed - current seed
    * @param tOffset - any offset to the time for this image
    * @param subDir - any directory to go deeper in Images Dir
    */
    public void saveArrayFileAsImg(String fixType, int t, int run,int scale,
            String type,int seed,int tOffset,String subDir){
        if(!subDir.contains("/")){subDir += "/";}
        
        String fname;
        if(systemIsNetwork){
            fname = dir.getConfigDirectory(instId)+getNetworkTempFilename(run, t, "");
        }else{    
            fname = dir.getConfigDirectory(instId)+getLatticeTempFilename(run, t, "");
        }

        int sizeImage = scale*L ;
        while(sizeImage < minImageSize){scale = scale*2;sizeImage = scale*L;}

        BufferedImage img = new BufferedImage(sizeImage, sizeImage, BufferedImage.TYPE_INT_RGB);  
        Graphics2D g = img.createGraphics(); 
        boolean fileFound = true;
        g.setColor(Color.LIGHT_GRAY);        
        g.fillRect(0, 0, img.getWidth(), img.getHeight());
        try{
            Scanner scanner = new Scanner(new File(fname)); Scanner scanner2;
            while(scanner.hasNextLine() && scanner.hasNextInt()) {
                scanner2 = new Scanner(scanner.nextLine());
                   
                int x,y;
                if(systemIsNetwork){
                    x = scanner2.nextInt();
                    y = ((int)((double)x/((double)param.L)))%param.L;
                    x = x%param.L;
                }else{
                    x = scanner2.nextInt();
                    y = scanner2.nextInt();
                }
                int fSpin = scanner2.nextInt();
                g.setColor(sq.getColor(fSpin));
                drawSpin(x, y,scale,0,1.0,g);
                //g.fillRect(x*scale, y*scale, scale, scale);    			
            }
        } catch (FileNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
            fileFound = false;   
        }
        
        Scanner scan2;
        try {
            scan2 = new Scanner(new File(dir.getFixedDirectory()+"fixed"+fixType+".txt"));
            Scanner scanner2;
            while(scan2.hasNextLine() && scan2.hasNextInt()) {
                scanner2 = new Scanner(scan2.nextLine());
                int x,y;
                if(systemIsNetwork){
                    x = scanner2.nextInt();
                    y = ((int)((double)x/((double)param.L)))%param.L;
                    x = x%param.L;
                }else{
                    x = scanner2.nextInt();
                    y = scanner2.nextInt();
                }
                int fSpin = scanner2.nextInt();
                g.setColor(sq.getColorFixed(fSpin));
                drawSpin(x, y,scale,0,1.0,g);
                //g.fillRect(x*scale, y*scale, scale, scale);    		
            }
        } catch (FileNotFoundException ex) {
            Logger.getLogger(DataSaver.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        t = t-tOffset;
        
        if(systemIsNetwork){
            fname = dir.getImageDirectory()+subDir+getNetworkImgFilename(run, t, seed);
        }else{
            fname = dir.getImageDirectory()+subDir+getLatticeImgFilename(run, t, seed);
        }
        
        createDirIfUnInit(fname);
        File output = new File(fname);
        if(fileFound){
            try {
                ImageIO.write(img, type, output);
            } catch (IOException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }		
        }
    }

    private synchronized void drawSpin(int i, int j,int scale,int ImageTextOffset,double cirR, Graphics2D g){
        g.fillOval(i * scale,j * scale + ImageTextOffset, (int)(scale*cirR)-1, (int)(scale*cirR)-1);
    }
    
    /**
    *      saveIntArrayAsImg saves an integer array  and 
    *  does this in the form of an image in png format.
    * 
    * @param lat - Integer array of lattice
    * @param fname - output file name
    */
    public void saveIntArrayAsImg(int[][] lat, String fname){
        fname = dir.getConfigDirectory(instId)+fname;
        int scale =minScale;
        int sizeImage = scale*L;
        while(sizeImage < minImageSize){scale = scale*2;sizeImage = scale*L;}

        BufferedImage img = new BufferedImage(sizeImage, sizeImage, BufferedImage.TYPE_INT_RGB);  
        Graphics2D g = img.createGraphics();          
        g.setColor(Color.LIGHT_GRAY);        
        g.fillRect(0, 0, img.getWidth(), img.getHeight());
        int fSpin;
        for(int x=0;x<L;x++){ for(int y=0;y<L;y++){
            fSpin = lat[x][y];
            g.setColor(sq.getColor(fSpin));
                drawSpin(x, y,scale,0,1.0,g);
        }}

        createDirIfUnInit(fname);
        File output = new File(fname);
        try {
            ImageIO.write(img, "png", output);
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }		
    }

    /**
    *      saveIntArrayAsImg saves an integer array  and 
    *  does this in the form of an image in png format.
    * 
    * @param lat - Integer array of lattice
    * @param fname - output file name
    */
    public void saveFloatArrayAsImgInImages(float[][] lat, String fname){
    	saveFloatArrayAsImgInImages(lat, fname,"");
    }
    
    /**
    *      saveIntArrayAsImg saves an integer array  and 
    *  does this in the form of an image in png format.
    * 
    * @param lat - Integer array of lattice
    * @param fname - output file name
    * @param outString - string to post on image
    */
    public void saveFloatArrayAsImgInImages(float[][] lat, String fname,String outString){
        fname = dir.getImageDirectory()+fname;
        int scale = minScale;
        int sizeImage = scale*L ;
        while(sizeImage < minImageSize){scale = scale*2;sizeImage = scale*L;}
        BufferedImage img = new BufferedImage(sizeImage, sizeImage, BufferedImage.TYPE_INT_RGB);  
        Graphics2D g = img.createGraphics();
        g.setColor(Color.LIGHT_GRAY);        
        g.fillRect(0, 0, img.getWidth(), img.getHeight());
        
        float fSpin;float rcol;float gcol;float bcol;
        for(int x=0;x<L;x++){ for(int y=0;y<L;y++){
            fSpin = lat[x][y];
            Color col = colorRange(fSpin);
            g.setColor(col);
            drawSpin(x, y,scale,0,1.0,g);
        }}
        
        g.setColor(Color.BLACK);
        Font font = new Font("Courier New", Font.BOLD, 32);
        g.setFont(font);
        // 70 is the x coordinate of text
        g.drawString(outString, ((int) L*scale*4/10), ((int) 90*7/10));


        createDirIfUnInit(fname);
        File output = new File(fname);
        try {
            ImageIO.write(img, "png", output);
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
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
    
    /**
    *      getImageFromArray saves an integer array and 
    *  does this in the form of an image in png format.
    * 
    * @param lat - Integer array of lattice
    */
    public BufferedImage getImageFromArray(float[][] lat){
        int scale =minScale;
        int sizeImage = scale*L ;
        while(sizeImage < minImageSize){scale = scale*2;sizeImage = scale*L;}

        BufferedImage img = new BufferedImage(sizeImage, sizeImage, BufferedImage.TYPE_INT_RGB);  
        Graphics2D g = img.createGraphics();
        g.setColor(Color.LIGHT_GRAY);        
        g.fillRect(0, 0, img.getWidth(), img.getHeight());
        
        float fSpin;float rcol;float gcol;float bcol;
        for(int x=0;x<L;x++){ for(int y=0;y<L;y++){
            fSpin = lat[x][y];
            Color col = colorRange(fSpin);
            g.setColor(col);
            drawSpin(x, y,scale,0,1.0,g);
        }}            
        return img;
    }
    
    /**
    *      saveIntArrayAsImg saves an integer array  and 
    *  does this in the form of an image in png format.
    * 
    * @param lat - Integer array of lattice
    * @param fix - fixed array of lattice
    * @param fname - ouput file name
    */
    public void saveIntArrayAsImg(int[][] lat, int fixed[][], String fname){
        fname = dir.getConfigDirectory(instId)+fname;

        int scale = minScale;
        int sizeImage = scale*L ;
        while(sizeImage < minImageSize){scale = scale*2;sizeImage = scale*L;}

        BufferedImage img = new BufferedImage(sizeImage, sizeImage, BufferedImage.TYPE_INT_RGB);  
        Graphics2D g = img.createGraphics(); 
        g.setColor(Color.LIGHT_GRAY);        
        g.fillRect(0, 0, img.getWidth(), img.getHeight());
        int fSpin;int fix;
        for(int x=0;x<L;x++){ for(int y=0;y<L;y++){
            fSpin = lat[x][y];
            fix = fixed[x][y];
            if(fix == 1){
                g.setColor(sq.getColorFixed(fSpin));
            }else{
                g.setColor(sq.getColor(fSpin));
            }
            drawSpin(x, y,scale,0,1.0,g);
        }}

        createDirIfUnInit(fname);
        File output = new File(fname);
        try {
            ImageIO.write(img, "png", output);
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }		
    }

    
    /**
    *      deleteConfigNotNucleationRange delete all configurations except
    *  for the given time and the time one step before and after this given time.
    * 
    * @param run - run to delete configurations for
    * @param time - time to exclude from deletions
    */
    public void deleteConfigNotNucleationRange(int run, int time){
        File path = new File(dir.getConfigDirectory(instId));
        File[] files = path.listFiles();
        String fname = new String("-Run-");String fname2 = new String("-Run-");
        String fname3 = new String("-Run-");
        if(run>0){
            fname = dir.getConfigDirectory(instId)+"Config-t-"+time+fname+run+".txt";
            time  = time+1;
            fname2 = dir.getConfigDirectory(instId)+"Config-t-"+time+fname+run+".txt";
            time  = time-2;
            fname3 = dir.getConfigDirectory(instId)+"Config-t-"+time+fname+run+".txt";
        }else{
            fname = dir.getConfigDirectory(instId)+"Config-t-"+time+".txt";
            time  = time+1;
            fname2 = dir.getConfigDirectory(instId)+"Config-t-"+time+".txt";
            time  = time-2;
            fname3 = dir.getConfigDirectory(instId)+"Config-t-"+time+".txt";	
        }
        
        // replaces forward slash with backslash
        fname = fname.replace('/', '\\');

        for(int i=0; i<files.length; i++) {
            if((files[i].getPath()).equalsIgnoreCase(fname)
                    ||(files[i].getPath()).equalsIgnoreCase(fname2)
                    ||(files[i].getPath()).equalsIgnoreCase(fname3)){
            }else{files[i].delete();}
        }
    }

    /**
    *       deleteConfigArrays deletes all the configurations files for the current
    *   instance of the program.
    */
    public void deleteConfigArrays(){
        File path = new File(dir.getConfigDirectory(instId));
        File[] files = path.listFiles();

        for(int i=0; i<files.length; i++) {
            String ext= files[i].getAbsolutePath();
            ext = ext.substring(ext.length()-3, ext.length());
            if(ext.equalsIgnoreCase("txt")){files[i].delete();}
        }
    }


    /**
    *        saveConfig2D3D saves the given lattice into the configuration directory
    *   as a text file
    * 
    * @param lat - input lattice
    * @param run - run for this given lattice
    * @param time - time for current state of lattice
    */
    public void saveConfig(LatticeMagInt lat,int run, int time){
        saveConfig2D3D(lat,run,time,"");
    }
  
    /**
    *       saveConfigPermanent saves the given lattice as an image file into
    *   the saved config directory in the configuration directory
    * 
    * @param lat - input lattice
    * @param run - run for this given lattice
    * @param time - time for current state of lattice
    * @param post - postfix to add to the filename
    */
    public void saveConfigPermanent(LatticeMagInt lat,int run, int time,String post)
    {
        systemIsNetwork =false;
        String fname = new String(dir.getSavedConfigDirectory(instId));
        fname += getLatticePermFilename(run, time, post);
        try {
            // retrieve image
            BufferedImage bi = lat.getSystemImg();
            createDirIfUnInit(fname);
            File outputfile = new File(fname);
            ImageIO.write(bi, "png", outputfile);
        } catch (IOException e){
        }
    }
    
    
    /**
    *       saveConfigPermanent saves the given lattice as an image file into
    *   the saved config directory in the configuration directory
    * 
    * @param lat - input lattice
    * @param run - run for this given lattice
    * @param time - time for current state of lattice
    * @param post - postfix to add to the filename
    */
    public void saveConfigPermanent(SimSystem sys,int run, int time,String post,int nfixed,String configSavePost){
        if(sys.getClass().getName().contains("etwork")){
            saveNetworkPermanent((Network)sys,run,time,post,nfixed,configSavePost);
        }else if(sys.getClass().getName().contains("attice")){
            saveConfigPermanent((LatticeMagInt)sys,run,time,post,nfixed,configSavePost);
        }
    }
    
    public void saveConfigPermanent(SimSystem lat,int run, int time,String post,int nfixed){
        saveConfigPermanent(lat,run, time,post,nfixed,"");
    }

    public void saveNetworkPermanent(Network net,int run, int time,String post,int nfixed,String confPost){
        saveNetworkConfig(post, run, time, net.getNetworkValuesArr(), 
                net.getNetworkLinks(), net.getNetworkFixed(),nfixed,confPost);
    }
    /**
    *       saveConfigPermanent saves the given lattice as an image file into
    *   the saved config directory in the configuration directory
    * 
    * @param lat - input lattice
    * @param run - run for this given lattice
    * @param time - time for current state of lattice
    * @param post - postfix to add to the filename
    */
    public void saveConfigPermanent(LatticeMagInt lat,int run, int time,String post,int nfixed){
        saveConfigPermanent(lat,run, time,post,nfixed,"");
    }
    
    /**
    *       saveConfigPermanent saves the given lattice as an image file into
    *   the saved config directory in the configuration directory
    * 
    * @param lat - input lattice
    * @param run - run for this given lattice
    * @param time - time for current state of lattice
    * @param post - postfix to add to the filename
    */
    public void saveConfigPermanent(LatticeMagInt lat,int run, int time,String post,int nfixed,String configSavePost)
    {
        String fname = new String(dir.getSavedConfigDirectory(instId)+"Fixed_"+nfixed+"_"+configSavePost+"/");

        fname += getLatticePermFilename(run, time, post);
        
        // Check if directory exists if not create
        File f = new File(fname);
        createDirIfUnInit(fname);
        
        // push into directory

        try {// retrieve image
            BufferedImage bi = lat.getSystemImg();
            File outputfile = new File(fname);
            ImageIO.write(bi, "png", outputfile);
        } catch (IOException e) {
            
        }
    }
    
    /**
    *        saveConfig2D3D saves the given lattice into the configuration directory
    *   as a text file
    * 
    * @param lat - input lattice
    * @param run - run for this given lattice
    * @param time - time for current state of lattice
    * @param post - postfix to add to the filename
    */
    public void saveConfig2D3D(LatticeMagInt lat,int run, int time,String post){
        systemIsNetwork =false;
        String fname = dir.getConfigDirectory(instId)+getLatticeTempFilename(run, time, post);
        try {
            PrintStream out = new PrintStream(new FileOutputStream(
                fname));
            for (int i = 0; i < lat.getLength(); i++){for (int j = 0; j < lat.getLength(); j++){
                if(param.D == 3){
                    for (int k = 0; k < lat.getLength(); k++){out.println(i + "   " + j + "  "+k + "  "+   lat.getValue(i, j,k));}
                }else{out.println(i + "   " + j + "  "+   lat.getValue(i, j));}
            }}
            out.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    
    public void saveNetworkVals(Network net,int run, int time,String post){
        systemIsNetwork = true;
        String fname = dir.getConfigDirectory(instId)+getNetworkTempFilename(run, time, post);
        try {
            PrintStream out = new PrintStream(new FileOutputStream(
                fname));
            for (int i = 0; i < net.getN(); i++){
                out.println(i + "   " +   net.getValue(i)+"    "+net.getNetworkSeed());
            }
            out.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }
        
    /**
    *        saveConfig2D3D saves the given lattice into the configuration directory
    *   as a text file
    * 
    * @param sys - input system
    * @param run - run for this given lattice
    * @param time - time for current state of lattice
    * @param post - postfix to add to the filename
    */
    public void saveConfig(SimSystem sys,int run, int time,String post){
        if(sys.getClass().getName().contains("etwork")){
            saveNetworkVals((Network)sys,run,time,post);
        }else if(sys.getClass().getName().contains("attice")){
            saveConfig2D3D((LatticeMagInt)sys,run,time,post);
        }
    }
    
    public void saveConfig(SimSystem sys,int run, int time){
        saveConfig(sys,run,time,"");
    }
    
    /**
    *       saveSeed saves the given seed as a text file in the JISim root directory
    * 
    * @param seed - input seed number
    */
    public void saveSeed(int seed){
        String fname = "seed.txt";
        fname = dir.getRootDirectory()+fname;
        try {
            PrintStream out = new PrintStream(new FileOutputStream(fname));
            out.println(seed);
            out.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    /**
    *       saveStringDataKey saves a given string array as a text file with 
    *   the given filename. This is supposed to be used to save a key for the
    *   output data.
    * 
    * @param data - string array to be used as key for data output
    * @param name - filename 
    */
    public void saveStringDataKey(String[] data, String name){
        String fName = dir.getDataDirectory()+name;
        int size = data.length;
        // Data is being appended since data is part of a set
        try {
            PrintStream out = new PrintStream(new FileOutputStream(fName));
            for (int i = 0; i < size; i++){
                out.print(data[i]+ "     ");
            }
            out.println();
            out.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

    }
    
    /**
    * 
    *       saveDoubleArrayData saves the given double array as a text file in
    *   the data directory.
    * 
    * @param data - double array that is input data
    * @param name - filename of output file
    */
    public void saveDoubleArrayData(double[] data, String name){
        String fName = dir.getDataDirectory()+name;
        int size = data.length;
        createDirIfUnInit(fName);
        // Data is being appended since data is part of a set
        try {
            PrintStream out = new PrintStream(new FileOutputStream(fName,true));
            for (int i = 0; i < size; i++){
                out.print(data[i]+ "     ");
            }
            out.println();
            out.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    /**
    * 
    *       saveDoubleArrayListData saves the given double arraylist as a text file in
    *   the data directory.
    * 
    * @param data - double array that is input data
    * @param name - filename of output file
    */
    public void saveDoubleArrayListData(ArrayList<Double> data, String name){
        String fName = dir.getDataDirectory()+name;
        createDirIfUnInit(fName);
        // Data is being appended since data is part of a set
        try {
            PrintStream out = new PrintStream(new FileOutputStream(fName,true));
            for (int i = 0; i < data.size(); i++){
                out.print(data.get(i)+ "     ");
            }
            out.println();
            out.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }
    
    /**
    * 
    *       saveIntArrayListData saves the given int arraylist as a text file in
    *   the data directory.
    * 
    * @param data - double array that is input data
    * @param name - filename of output file
    */
    public void saveIntArrayListData(ArrayList<Integer> data, String name){
        String fName = dir.getDataDirectory()+name;
        createDirIfUnInit(fName);
        // Data is being appended since data is part of a set
        try {
            PrintStream out = new PrintStream(new FileOutputStream(fName,true));
            for (int i = 0; i < data.size(); i++){
                out.print(data.get(i)+ "     ");
            }
            out.println();
            out.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }
    
    /**
    * 
    *       saveIntLLData saves the given integer linked list as a text file in
    *   the data directory.
    * 
    * @param data - integer linked list that is input data
    * @param name - filename of output file
    */
    public void saveIntLLData(LinkedList<Integer> data, String name){
        String fName = dir.getDataDirectory()+name;
        int size = data.size();
        try {
            PrintStream out = new PrintStream(new FileOutputStream(fName,true));
            for (int i = 0; i < size; i++){ out.println(data.remove()); }
            out.println();
            out.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    
    /**
    * 
    *       saveIntLLData saves the given double linked list as a text file in
    *   the data directory.
    * 
    * @param data - double linked list that is input data
    * @param name - filename of output file
    */
    public void saveDoubleLLData(LinkedList<Double> data, String name){
        String fName = dir.getDataDirectory()+name;
        createDirIfUnInit(fName);
        int size = data.size();
        try {
            PrintStream out = new PrintStream(new FileOutputStream(fName,true));
            for (int i = 0; i < size; i++){ out.println(data.remove()); }
            out.println();
            out.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }
    
    public void saveNetworkConfig(String post,int run, int time, AtomicIntegerArray nodeIsingVals,
        HashMap<Integer,ArrayList<Integer>> networkBonds,HashMap<Integer,Boolean> networkFixed,
        int nfixed, String configSavePost){    
        systemIsNetwork = true;
        String fname = new String(dir.getSavedConfigDirectory(instId)
                +"Fixed_"+nfixed+"_"+configSavePost+"/");
        fname += getNetworkPermFilename(run, time, post);
        createDirIfUnInit(fname);
        
        try {
            PrintStream out = new PrintStream(new FileOutputStream(fname));
            for(int u = 0;u < nodeIsingVals.length();u++){
                String outString = ""+u+"   "+nodeIsingVals.get(u);
                if(networkFixed.get(u) == null){
                    outString = outString +"   "+0;
                }else{
                    if(networkFixed.get(u) == true){
                        outString = outString +"   "+1;}else{
                        outString = outString +"   "+0;
                    }
                }
                ArrayList<Integer> bonds = networkBonds.get(u);
                for(int k = 0; k < bonds.size();k++){
                    outString = outString +"   "+bonds.get(k);
                }
                out.println(outString);    
            }                    
            out.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }
    
    public void createDirIfUnInit(String filename){
        if(!filename.contains("/")){return;}
        String[] temp = filename.split("/");
        String dir = filename.substring(0,filename.length()-temp[temp.length-1].length());
        if((new File(dir)).exists()){}else{
            temp = dir.split("/");   
            String currDir = "/";
            for(int u = 0; u < temp.length;u++){
                currDir += temp[u]+"/";
                File f = new File(currDir);    
                f.mkdir();  
            }
        }
    }
    
    public String getLatticeTempFilename(int run, int time, String post){
        String fname = new String("-Run-");
        if(run>0){
            fname = "Config-t-"+time+fname+run+post+".txt";
        }else{
            fname = "Config-t-"+time+post+".txt";
        }
        return fname;
    }
    
    public String getLatticePermFilename(int run, int time, String post){
        String fname="";
        if(run>0){
            fname = fname+
                "Config-t-"+time+"-Run-"+run+post+"-temp-"+param.temperature
                +"-field-"+param.hField+"-R-"+param.R+".png";
        }else{
            fname = fname+
                "Config-t-"+time+post+"-temp-"+param.temperature
                +"-field-"+param.hField+"-R-"+param.R+".png";        	
        }
        return fname;
    }
    public String getLatticeImgFilename(int run, int t, int seed){
        String fname = "";
        if(seed!=0){
            fname = "Config-t-"+t+"-Run-"+run+"-field-"+param.hField+"-temp-"
                    +param.temperature+"-R-"+param.R+"-Seed-"+seed+".png";
        }else if(seed==0 && run<=0){
            fname = dir.getImageDirectory()+"Config-t-"
                +t+"-field-"+param.hField+".png";            
        }else{
            fname = "Config-t-"
                    +t+"-Run-"+run+"-field-"+param.hField+".png";            
        }
        return fname;
    }
    
    public String getNetworkImgFilename(int run, int t, int seed){
        String fname = "";
        if(seed!=0){
        fname = "NetworkConfig-t-"
                +t+"-Run-"+run+"-Seed-"+seed+"-field-"+param.hField+".png";
            }else if(seed==0 && run<=0){
                    fname = dir.getImageDirectory()+"Config-t-"
                        +t+"-field-"+param.hField+".png";            
            }else{
        fname = "NetworkConfig-t-"
                +t+"-Run-"+run+"-field-"+param.hField+".png";            
        }
        return fname;
    }
    public String getNetworkPermFilename(int run, int time, String post){
        String fname="";
        if(run>0){
            fname = fname+
                "NetworkConfig-t-"+time+"-Run-"+run+post+"-temp-"+param.temperature
                +"-field-"+param.hField+"-R-"+param.R+".txt";
        }else{
            fname = fname+
                "NetworkConfig-t-"+time+post+"-temp-"+param.temperature
                +"-field-"+param.hField+"-R-"+param.R+".txt";        	
        }
        return fname;
    }
    public String getNetworkTempFilename(int run, int time, String post){
        String fname = new String("-Run-");
        if(run>0){
            fname = "NetworkConfig-t-"+time+fname+run+post+".txt";
        }else{
            fname = "NetworkConfig-t-"+time+post+".txt";
        }
        return fname;
    }   
    
    /**
    *      getDoubleDataFromFile gets data from a data file consisting of a single
    *   line.
    * 
    * @param fname - file name
    */
    public ArrayList<Double> getDoubleDataFromFileSingleLine(String fname){
        ArrayList<Double> dat = new ArrayList<Double>();
        //read file first    
        try{
            Scanner scanner = new Scanner(new File(fname));
            while(scanner.hasNextDouble()){
                dat.add(scanner.nextDouble());
            }
        }catch (FileNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        return dat;
    }
    
    // test the class
    public static void main(String[] args) {
        DataSaver dat = new DataSaver();
        SimpleLattice lat = new SimpleLattice(4,1,2,2,1);
        double[] data = new double[4];
        data[1]=1; data[0]=0;
        data[2]=2; data[3]=3;
        //dat.deleteConfigArrays();        
        System.out.println("DataSaver | Done!");
    }
}
