package JISim;
/**
 * 
 *    @(#)   JISimCommandLine
 */  
import Backbone.Util.DataSaver;
import Backbone.Util.DirAndFileStructure;
import Backbone.Util.ParameterBank;
import Examples.MCFieldExplorer;
import Examples.MCSimulationMovieRuns;
import java.io.File;
import java.util.Scanner;

/**
*      JISimCommandLine is a console line program for running basic example
*   simulations and also perform tasks like clearing instances file and making
*   fixed configurations
* 
* <br>
* 
* @author      James Silva <jbsilva @ bu.edu>                 
* @since       2012-02
*/
public class JISimCommandLine {
    Scanner scan = new Scanner(System.in);
    	
    /**
    *          printWelcomeScreen outputs a welcome screen to the console
    */
    public void printWelcomeScreen(){
        System.out.println("*********************************************");
        System.out.println("*********************************************");
        System.out.println("********W*e*l*c*o*m*e************************");
        System.out.println("*********************************************");
        System.out.println("************t*o******************************");
        System.out.println("*********************************************");
        System.out.println("**********J*I*S*i*m**************************");
        System.out.println("*********************************************");
        System.out.println("*********************************************");
    }
	
    /**
    *      assertSetup makes sure the setup of the program is complete and 
    *  will send user to setup prompt if not setup.
    */
    public void assertSetup(){
        File f = new File("directory.txt");
        if(f.exists()){}else{
        SetupJISim.main(null);    
        }
    }

    /**
    *          doAction performs the action previously selected by user then ask 
    *      user if a new action must be performed
    * 
    * @param act - action number selected 
    * @return - true if need to perform another action
    */
    public boolean doAction(int act){
            DirAndFileStructure d = new DirAndFileStructure();
            SetupJISim set= new SetupJISim();
            System.out.println("JISim set in "+d.getRootDirectory());
            set.setDirectory(d.getRootDirectory());
            
        switch (act) {
        // add new simulation actions here
        case 1: MCSimulationIntervention.main(null);break;
        case 2: {if(MakeFixedConfig.movieSetExists()){}else{MakeFixedConfig.makeMovieSet();}
                        set.makeBasicParameters("Movie");set.makeBasicSimProcess("Movie");
                        MCSimulationMovieRuns.main(null);break;}
        case 3: {if(MakeFixedConfig.movieSetExists()){}else{MakeFixedConfig.makeMovieSet();}
                        set.makeBasicParameters("Explorer");set.makeBasicSimProcess("Explorer");
                        MCFieldExplorer.main(null);break;}
        case 4: MakeFixedConfig.main(null);break;
        case 5: {System.out.println(" What file post fix parameter file?");String post = scan.next();
                                ParameterBank p = new ParameterBank(post);post=p.postFix;
                                System.out.println("what is the new field?");
                                p.changeField(scan.nextDouble(), post);
                        }break;
        case 6: {System.out.println(" What file post fix parameter file?");String post = scan.next();
                ParameterBank p = new ParameterBank(post);post=p.postFix;
                System.out.println("what is the new temp?");
                p.changeTemp(scan.nextDouble(), post);
                }break;
        case 7:  {DataSaver dsave = new DataSaver();
        dsave.clearInstances();
        System.out.println("cleared");}break;    

        default: System.out.println("exiting.");return false;
        }
        boolean ConfirmNotDone = true; 
        String confirm;
        boolean continueJISim=false;
        while(ConfirmNotDone){
                        System.out.println("Continue with another simulation?(y or n)");		
                        confirm = scan.next();
        if(confirm.equalsIgnoreCase("y") ){
                ConfirmNotDone=false;continueJISim=true;
        }else if(confirm.equalsIgnoreCase("n")){
                        ConfirmNotDone=false;
        }else{
        System.out.println("y or n");}
                }
        return continueJISim;
    }
	
    /**
    *      getAction queries the user as to what action is to be performed and returns
    *  the action number of this action
    * 
    * @return action number of action selected by user
    */
    public int getAction(){
        boolean ResponseNotDone = true; 
        boolean ConfirmNotDone = true; 
        String confirm;
        int act=0;
        while(ResponseNotDone){
            // Add new simulations here
            System.out.println("Choose simulation to run");
            System.out.println("0) Exit JISim.");
            System.out.println("1) Intervention Method on the Ising Model.");
            System.out.println("2) Create videos of Ising simulations on different fixed" +
                            " spin arrangements.");
            System.out.println("3) Explore the field range that is best for a " +
                            "particular arrangement.");
            System.out.println("4) Make fixed spin configurations.");
            System.out.println("5) Edit your parameters file. Change h field.");
            System.out.println("6) Edit your parameters file. Change temperature.");
            System.out.println("7) Clear your instances file.");

            act = scan.nextInt();
            while(ConfirmNotDone){
                    System.out.println("Is this the right action? (y or n)");		
                    System.out.println("action: "+act);
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

            return act;
    }
	
    public static void main(String[] args) {
        JISimCommandLine jcmd = new JISimCommandLine();
        jcmd.assertSetup();
        boolean continueCmd=true;
        while(continueCmd){
        continueCmd = jcmd.doAction(jcmd.getAction());}
        System.out.println("Bye.");
    }
}
