package Backbone.Util;

/**
*    @(#) DirAndFileStructure 
*/  

import JISim.MCSimulation;
import JISim.SetupJISim;
import java.io.*;
import java.nio.channels.FileChannel;
import java.util.Scanner;
import java.util.logging.Level;
import java.util.logging.Logger;

/** 
* 	This class keeps all information on the project directory and 
*   structure. This is only useful after setting the structure up by using 
*   SetupJISim class.
*  <br>
*   
*  @see JISim.SetupJISim
*   
* @author      James B Silva <jbsilva @ bu.edu>                 
* @since       2011-10   
*/
public class DirAndFileStructure {
    private String directory;	
    public DirAndFileStructure(){	
        try {
            setDirectory();
        } catch (IOException e) {
            // TODO Auto-generated catch block
           e.printStackTrace();
        }
    }

    /**
    *         setDirectory reads the directory which was setup using setRootDirectory class.
    *         
    *         Alternatively one could setup a txt file called directory.txt and include it in the
    *         java project directory
    *
    *  
    *  @see JISim.SetupJISim
    *
    */
    public void setDirectory() throws IOException {
        FileReader fin = new FileReader("directory.txt");
        Scanner src = new Scanner(fin);
        directory = src.next();
    }

    /**
    *         getRunDirectory gives the directory for extra Run data not covered in the other directory.
    *  
    *
    */
    public String getRunDirectory(){
        return directory+"Run/";
    }

    /**
    *         getFixedDirectory gives the directory for fixed lattice value files.
    *  
    *
    */
    public String getFixedDirectory(){
        return directory+"Fixed/";
    }

    /**
    *         getSettingsDirectory gives the directory for program settings files.
    *  
    *
    */
    public String getSettingsDirectory(){
        return directory+"Settings/";
    }

    /**
    *         assertFixedFileExist checks for fixed file existence and makes
    *   fixed file if it doesnt exist.
    *  
    *
    */
    public void assertFixedFileExist(String postfix){
        // Does file for this instance exists
        if(!((new File(directory+"Fixed/"+"fixed"+postfix+".txt")).exists())){
            System.out.println("DirAndFileStructure Fixed file does not exist. Making : "+
                    directory+"Fixed/"+"fixed"+postfix+".txt");
            String fname = directory+"Fixed/"+"fixed"+postfix+".txt";
            try {
                PrintStream out = new PrintStream(new FileOutputStream(
                        fname));
                out.println();
                out.println();
                out.close();
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }
        }
    }
    
    /**
    *         getConfigDirectory gives the directory for all configuration data.
    *  
    *
    */
    public String getConfigDirectory(int id){
        // Does file for this instance exists
        if(!((new File(directory+"Config/"+"Instance "+id)).exists())){
            File f = new File(directory+"Config/"+"Instance "+id+"/");
            f.mkdir();  
        }
        return directory+"Config/"+"Instance "+id+"/";
    }
    
    /**
    *         getLatticeConfig gives the lattice data for given input.
    *  
    *	@param id - instance id
    *	@param t - time of lattice
    *	@param run - run for this lattice
    *   @param post - any postscript for file
    */
    public String getLatticeConfig(int id, int t,int run, String post){
        // Does file for this instance exists
        String fname = "-Run-";
        if(run>0){
            fname = getConfigDirectory(id)+"Config-t-"+t+fname+run+post+".txt";}
        else{
            fname = getConfigDirectory(id)+"Config-t-"+t+post+".txt";	
        }
        return fname;
    }
    
    /**
    *         getSavedConfigDirectory gives the directory to save configuration data.
    *  
    *
    */
    public String getSavedConfigDirectory(int id){
        String postfix="";
        if(id > 0){ postfix = postfix+id; }
            
        // Does file for this instance exists
        if(!((new File(directory+"Config/"+"Saved Configs"+postfix)).exists())){
            File f = new File(directory+"Config/"+"Saved Configs"+postfix+"/");
            f.mkdir();  
        }
        return directory+"Config/"+"Saved Configs"+postfix+"/";
    }

    /**
    *         getDataDirectory gives the directory for actual data created by the simulation.
    *  
    *
    */
    public String getDataDirectory(){ return directory+"Data/"; }

    /**
    *         getDataDirectory gives the directory for actual data created by the simulation.
    *  
    *
    */
    public String getDataDirectory(String SubDir){
       // Does file for this instance exists
        if(!((new File(getDataDirectory()+SubDir+"/")).exists())){
            File f = new File(getDataDirectory()+SubDir+"/");
            f.mkdir();  
        }
        return getDataDirectory()+SubDir+"/";
    }

    /**
    *         getVideoDirectory gives the directory for actual video created by the simulation.
    *  
    *
    */
    public String getVideoDirectory(){
        return directory+"Videos/";
    }

    /**
    *         getVideoDirectory gives the directory for actual images (usually not configuratiosn) created by the simulation.
    *  
    *
    */
    public String getImageDirectory(){
        return directory+"Images/";
    }

    /**
    *         getRootDirectory gives the directory for the the simulation project ie the directory that contains
    *         the directories for data/config/run/fixed.
    *  
    *
    */
    public String getRootDirectory(){return directory;}
    
    public void createTempSettingsFile(String inParam,String tempParam)  {
        String prefix = getSettingsDirectory()+"parameters";
        String fileType =".txt";
        //System.out.println("in: "+prefix+inParam+fileType);
        //System.out.println("out: "+prefix+tempParam+fileType);
        File sourceFile = new File(prefix+inParam+fileType);
        File destFile = new File(prefix+tempParam+fileType);
        try {
            if(!destFile.exists()) { destFile.createNewFile(); }
            FileChannel source = null;
            FileChannel destination = null;
            try {
                source = new FileInputStream(sourceFile).getChannel();
                destination = new FileOutputStream(destFile).getChannel();
                destination.transferFrom(source, 0, source.size());
            }
            finally {
                if(source != null) { source.close(); }
                if(destination != null) { destination.close(); }
            }
        } catch (IOException ex) {
              Logger.getLogger(DirAndFileStructure.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    public void deleteTempSettingsFile(String tFile){
        String prefix = getSettingsDirectory()+"parameters";
        String fileType =".txt";
        tFile = prefix+tFile+fileType;
        //System.out.println("DELETING FILE"+tFile);
        File f = new File(tFile);
        if(f.exists()){ f.delete(); }
    }
    
    /*
    *       checkForSettingsFiles checks that all necessary settings are present
    *    and makes them in default form if they are not.
    * 
    *   @param mc - mcsimulation to check for existance of settings files
    */
    public void checkForSettingsFiles(MCSimulation mc){
        // Check for Files 
        DirAndFileStructure dir = new DirAndFileStructure();
        SetupJISim set = new SetupJISim(dir.getRootDirectory());
        // Check for param file
        File f = new File(dir.getSettingsDirectory()+"parameters"+mc.getParamPostFix());     
        if(!f.exists()){ set.makeBasicParameters(mc.getParamPostFix()); }
        // Sim Process
        f = new File(dir.getSettingsDirectory()+"SimProcess"+mc.getSimPostFix());     
        if(!f.exists()){ set.makeBasicSimProcess(mc.getSimPostFix()); }
        // Fixed Config File
        String fname =dir.getFixedDirectory()+"fixed"+mc.getFixedPostFix()+".txt";
        f = new File(fname);
        if(!f.exists()){
            try {
                PrintStream out = new PrintStream(new FileOutputStream(
                    fname));
                out.println(); out.close();
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }
        }
    }
}
