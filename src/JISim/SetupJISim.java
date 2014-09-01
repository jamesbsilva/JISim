package JISim;

/**
 * 
 *   @(#) setupJISim
 */  
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.Scanner;


 /** 
 *   A class for writing a configuration file for the directory for simulation. 
 *  It also writes the necessary directories. It simply sets up everything.
 *  <br>
 * 
 * 
 * @author      James Silva <jbsilva @ bu.edu>                 
 * @since       2011-09    
 */



public class SetupJISim {

    // linux for myself /home/jbsilva/Research/JISim/
    // windows D:/Research/Workspace/JISim/
    public String directory =  "";

    /**
    * @param d - directory in which to run JISim
    */
    public SetupJISim(){
    }
    
    /**
     * @param d - directory in which to run JISim
     */
    public SetupJISim(String d){
        directory =d;
    }

    /**
    *      setDirectory creates the setup file for establishing the directory
    *  to save all information.
    */
    public void setDirectory(){
        String fname = new String("directory.txt");

        try {
            PrintStream out = new PrintStream(new FileOutputStream(
                fname));

                        out.println(directory);

            out.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    /**
    *      setDirectory creates the setup file for establishing the directory
    *  to save all information.
    * 
    *   @param dir - directory to save all information to.
    */
    public void setDirectory(String dir){
        String fname = "directory.txt";
        directory = dir;
        try {
            PrintStream out = new PrintStream(new FileOutputStream(
                fname));

                        out.println(directory);

            out.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }
	
    /**
    *      makeDirs creates the directory structure that needs to be made in the root directory
    */
    public void makeDirs(){
        File f = new File(directory+"/Config/");
        f.mkdir();
        f = new File(directory+"/Config/Saved Configs/");
        f.mkdir();
        f = new File(directory+"/Data/");
        f.mkdir();  
        f = new File(directory+"/Fixed/");
        f.mkdir();  
        f = new File(directory+"/Run/");
        f.mkdir();  
        f = new File(directory+"/Videos/");
        f.mkdir();  
        f = new File(directory+"/Images/");
        f.mkdir();  
        f = new File(directory+"/Settings/");
        f.mkdir();  
    
        String fname = directory+"/Fixed/fixed.txt";
            try {
                    PrintStream out = new PrintStream(new FileOutputStream(
                        fname));
                    out.println();
                    out.close();
                } catch (FileNotFoundException e) {
                    e.printStackTrace();
            }
            
        createInstancesFile();
      
    }
	
    /**
    *      createInstancesFile makes the instances file which keeps track
    *   of how many instances of the program are running
    */
    private void createInstancesFile(){
        
        String fname = directory+"/Settings/instances.txt";
        try {
                PrintStream out = new PrintStream(new FileOutputStream(
                    fname));
                out.println(0);
                out.close();
            } catch (FileNotFoundException e) {
                e.printStackTrace();
        }

    }
	
        
    /**
    *       makeBasicSimProcess makes a basic simulation process file and saves
    *   it in the settings directory.
    * 
    * @param type - post fix for the filename of outgoing simulation process file
    */
    public void  makeBasicSimProcess(String type){
        String fname = directory+"/Settings/SimProcess"+type+".txt";
        if((new File(fname)).exists()){}else{
        try {
            PrintStream out = new PrintStream(new FileOutputStream(
                fname));
            out.println("0 changeT 10.0");
            out.println("50 ChangeBackT");
            out.println("150 FlipField");

            out.close();
            } catch (FileNotFoundException e) {
                e.printStackTrace();
        }}

    }
        
    /**
    *       makeBasicParameters makes a basic parameters file and saves
    *   it in the settings directory.
    * 
    * @param type - post fix for the filename of outgoing parameters file
    */
    public void  makeBasicParameters(String type){
        String fname = directory+"/Settings/parameters"+type+".txt";
        if((new File(fname)).exists()){}else{
        try {
                PrintStream out = new PrintStream(new FileOutputStream(
                    fname));
            out.println(":Length: 200      ");
            out.println(":initialSpin: 1     ");
            out.println(":dilution: false   // useDilution is false");
            out.println(":heterogenous: true      // use Heterogenous System, true for yes");
            out.println(":range: 5         // use Long Range Interactions if R > 0");
            out.println(":geometry: 2        // Geometry type  2 for square lattice, 3 for hex");
            out.println(":dimension: 2             // dimensions of system ");
            out.println(":temperature: 1.9         //  temperature of system "); 
            out.println(":hfield: 0.45              // magnetic field of system  ");
            out.println(":jInteraction: 1.0         // interaction parameter of system ");
            out.println(":mcalgo: metropolis            // algorithm of simulation");
            out.println(":trigger: simplevalfit 85 55     // trigger for simulation any "
                    + "ensuing double are parameters for trigger");
            out.println(":triggerparam1: 80                  // parameter for trigger 1 "
                    + "(overides value given in trigger line ");
            out.println(":triggerparam2: 50                   // parameter for trigger 2 "
                    + "(overides value given in trigger line ");
            

            out.close();
            } catch (FileNotFoundException e) {
                e.printStackTrace();
        }}

    }
        
        
    public static void main(String[] args) {
        boolean ResponseNotDone = true; 
        boolean ConfirmNotDone = true; 
        String response="";
        String confirm;
        Scanner scan = new Scanner(System.in);
        
        while(ResponseNotDone){
            System.out.println("Input the full directory to create/setup JISim directories and files?");	
            response = scan.next();
            while(ConfirmNotDone){
                System.out.println("Is this the right directory? (y or n)");		
                System.out.println(response);
                confirm = scan.next();
            if(confirm.equalsIgnoreCase("y") ){
                ConfirmNotDone=false;
                ResponseNotDone=false;
            }else if(confirm.equalsIgnoreCase("n")){
                ConfirmNotDone=false;
            }else{
                System.out.println("y or n");}
            }
            ConfirmNotDone=true;
        }
        
        SetupJISim set = new SetupJISim();

        set.setDirectory(response);
        set.makeDirs();
        set.makeBasicSimProcess("");
        set.makeBasicParameters("");

        System.out.print("Done!");
    }
}
