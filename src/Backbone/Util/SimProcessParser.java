package Backbone.Util;

/**
 * 
 *    @(#) SimProcessParser
 */  
import Backbone.Algo.IsingMC;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.LinkedList;
import java.util.Queue;
import java.util.Random;
import java.util.Scanner;


 /** 
 *   A parser for defining the process for a given simulation.
 *  <br>
 * 
 *  @param fname - postfix for SimProcess file.
 * 
 * @author      James Silva <jbsilva @ bu.edu>                 
 * @since       2011-09    
 */

public final class SimProcessParser {

    private Queue<double[]> processTimeActions;
    private Queue<double[]> persistProcessTimeActions;
    private double[] currentProcess;
    private DataSaver dSave;
    private String type = "";
    private String paramPostfix = "";
    private ParameterBank param=null;

    public SimProcessParser(){
        this("");
    }


    /**  @param fname - postfix for SimProcess file.
    */
    public SimProcessParser(String fname){
        this(null,fname,"");
    }

    /**  @param fname - postfix for SimProcess file.
    */
    public SimProcessParser(SimProcessParser sim,String fname,String parampost){
        if(sim!=null){            
            type =fname;
            paramPostfix=parampost;
            persistProcessTimeActions = new LinkedList<double[]>(sim.getListOfActions());
            processTimeActions = new LinkedList<double[]>(sim.getListOfActions());
        }else{
            type =fname;
            paramPostfix=parampost;
            processTimeActions = new LinkedList<double[]>();
            readSimProcess(fname);
            persistProcessTimeActions = new LinkedList<double[]>(processTimeActions);
        }
    }

    /**  @param fname - postfix for SimProcess file.
    */
    public SimProcessParser(SimProcessParser sim,String fname,String parampost , ParameterBank pbank){
        if(sim!=null){            
            type =fname;
            paramPostfix=parampost;
            param = pbank;
            persistProcessTimeActions = new LinkedList<double[]>(sim.getListOfActions());
            processTimeActions = new LinkedList<double[]>(sim.getListOfActions());
        }else{
            type =fname;
            paramPostfix=parampost;
            param = pbank;
            processTimeActions = new LinkedList<double[]>();
            readSimProcess(fname);
            persistProcessTimeActions = new LinkedList<double[]>(processTimeActions);
        }
    }
    
    public LinkedList<double[]> getListOfActions(){return (new LinkedList<double[]>(persistProcessTimeActions));}

    /**
    * 
    *       nextSimProcessTime returns the next time that a simulation
    *   process will take place and removes this time from the stack
    * 
    * @return time of next simulation process step
    */
    public int nextSimProcessTime(){
        currentProcess = processTimeActions.remove();
        return (int)currentProcess[0];
    }

    /**
    *       processesLeft returns the amount of simulation process time points left.
    * 
    * @return simulation processes left 
    */
    public int processesLeft(){return processTimeActions.size();}

    /**
    *      simProcess performs the next process on the given monte carlo simulation.
    * 
    * @param m - monte carlo simulation algorithm to perform the sim process in
    * @return the updated monte carlo simulation after the simulation process action
    */
    public IsingMC simProcess(IsingMC m){
        int procId = (int) currentProcess[1];
        int time = (int) currentProcess[0];

        switch (procId) {
            case 1: m.flipField();break;
            case 2: dSave.saveConfig(m.getSimSystem(),m.getRun(),time);break;
            case 3: m.flipField();dSave.saveConfig(m.getSimSystem(),m.getRun(),time);break;
            case 4: m.flipField();dSave.saveSeed(m.getSeed());break;
            case 5: m.setSeed(getRandomInt());m.flipField();break;
            case 6: m.setSeed(getRandomInt());dSave.saveConfig(m.getSimSystem(),m.getRun(),time);break;
            case 7: m.flipField(); m.setSeed(getRandomInt());dSave.saveConfig(m.getSimSystem(),m.getRun(),time);break;
            case 8: dSave.saveConfig(m.getSimSystem(),m.getRun(),time);dSave.saveSeed(m.getSeed());break;
            case 9: dSave.saveSeed(m.getSeed());break;
            case 10: m.setSeed(getRandomInt());break;
            case 11: m.changeT(currentProcess[2]);break;
            case 12: m.changeTFlipField(currentProcess[2]);break;
            case 13: m.changeTandH(currentProcess[2],currentProcess[2]);break;
            case 14: m.changeTandH(getTemp(), -1*getField());break;
            case 15: m.changeH(currentProcess[2]);break;
            case 16: m.changeT(getTemp());break;
            case 17: m.changeH(getField());break;
            case 18: m.changeTFlipField(getTemp());break;
            default: break;
        }

        return m;
    }

    private double getTemp(){
        ParameterBank pbank;
        if(param == null){
            pbank = new ParameterBank(paramPostfix);
        }else{
            pbank = param;
        }
        return param.temperature;
    }

    private double getField(){
        ParameterBank pbank;
        if(param == null){
            pbank = new ParameterBank(paramPostfix);
        }else{pbank = param;}
        return param.hField;
    }
    
    /**
    *   getRandomInt gets a random integer  
    * 
    * @return random integer
    */
    public int getRandomInt(){
        Random Ran = new Random();
        return Ran.nextInt();
    }

    /**
    *       readSimProcess reads the simulation process file with the given
    *   filename
    * 
    * @param fname - postfix of simulation process filename to read
    */
    public void readSimProcess(String fname){
        DirAndFileStructure directory = new DirAndFileStructure();
        String dir = directory.getSettingsDirectory(); 
        String fileName = dir + "SimProcess"+fname+".txt";

        Scanner scanner;
        
        try {
            scanner = new Scanner(new File(fileName));

            // Going to process line by line
            scanner.useDelimiter("\n");

            while(scanner.hasNext()) {
                String y = scanner.next();
                // process the String
                String[] fullY = y.split("\\s");
                if(!(fullY.length==0)){
                // Second value is supposed to be a string
                y = fullY[1];

                // First value is supposed to be a time of action
                double x =  Double.valueOf((fullY[0].trim())).doubleValue();

                // If Action requires a parameter it gets read here and its the third part of line
                double yp=0;
                if(fullY.length>2){
                yp = Double.valueOf((fullY[2].trim())).doubleValue();}

                double[] temp = new double[3];
                temp[0]= x;

                if(y.equalsIgnoreCase("FlipField")){temp[1]=1;}
                else if(y.equalsIgnoreCase("SaveConfig")){temp[1]=2;}
                else if(y.equalsIgnoreCase("FlipFieldSaveConfig")){temp[1]=3;}
                else if(y.equalsIgnoreCase("SaveConfigFlipField")){temp[1]=3;}
                else if(y.equalsIgnoreCase("FlipFieldSaveSeed")){temp[1]=4;}
                else if(y.equalsIgnoreCase("SaveSeedFlipField")){temp[1]=4;}
                else if(y.equalsIgnoreCase("FlipFieldChangeSeed")){temp[1]=5;}
                else if(y.equalsIgnoreCase("FlipFieldChangeSeed")){temp[1]=5;}

                else if(y.equalsIgnoreCase("ChangeSeedSaveConfig")){temp[1]=6;}
                else if(y.equalsIgnoreCase("SaveConfigChangeSeed")){temp[1]=6;}

                else if(y.equalsIgnoreCase("FlipFieldChangeSeedSaveConfig")){temp[1]=7;}
                else if(y.equalsIgnoreCase("SaveConfigFlipFieldChangeSeed")){temp[1]=7;}
                else if(y.equalsIgnoreCase("ChangeSeedFlipFieldSaveConfig")){temp[1]=7;}
                else if(y.equalsIgnoreCase("ChangeSeedSaveConfigFlipField")){temp[1]=7;}
                else if(y.equalsIgnoreCase("FlipFieldSaveConfigChangeSeed")){temp[1]=7;}
                else if(y.equalsIgnoreCase("SaveConfigChangeSeedFlipField")){temp[1]=7;}


                else if(y.equalsIgnoreCase("DoNothing")){temp[1]=0;}
                else if(y.equalsIgnoreCase("SaveConfigSaveSeed")){temp[1]=8;}
                else if(y.equalsIgnoreCase("SaveSeedSaveConfig")){temp[1]=8;}

                else if(y.equalsIgnoreCase("SaveSeed")){temp[1]=9;}
                else if(y.equalsIgnoreCase("ChangeSeed")){temp[1]=10;}
                else if(y.equalsIgnoreCase("ChangeT")){temp[1]=11;temp[2]=yp;}
                else if(y.equalsIgnoreCase("sameTandH")){temp[1]=13;temp[2]=yp;}
                else if(y.equalsIgnoreCase("sameHandT")){temp[1]=13;temp[2]=yp;}
                else if(y.equalsIgnoreCase("changeH")){temp[1]=15;temp[2]=yp;}
                else if(y.equalsIgnoreCase("FlipFieldChangeBackTandH")){temp[1]=14;} 
                else if(y.equalsIgnoreCase("ChangeTFlipField")){temp[1]=12;temp[2]=yp;}
                else if(y.equalsIgnoreCase("FlipFieldChangeT")){temp[1]=12;temp[2]=yp;}
            else if(y.equalsIgnoreCase("FlipFieldChangeBackT")){temp[1]=18;
            }else if(y.equalsIgnoreCase("ChangeBackT") || y.equalsIgnoreCase("ChangeBackTemp")||
                    y.equalsIgnoreCase("ChangeBackTemperature")){
                temp[1]=16;           
            }else if(y.equalsIgnoreCase("ChangeBackH") || y.equalsIgnoreCase("ChangeBackField") ||
                    y.equalsIgnoreCase("ChangeBackHField")){
                temp[1]=17;            
            }
            processTimeActions.add(temp);

            }}

            //push a end line like last process which cannot occur.
            double[] temp = new double[2];
            temp[0]= 0;
            temp[1]=0;
            processTimeActions.add(temp);
        } catch (FileNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }


    /**
    *       nextPrintSimActionNumber prints out the action number of the next 
    *   simulation process.
    */
    public void printSimActions(){
        LinkedList<double []> acts = getListOfActions();
        for (int u = 0; u< acts.size();u++){
            printProcess(acts.get(u));
        }
    }

    
    /**
    *       nextPrintSimActionNumber prints out the action number of the next 
    *   simulation process.
    */
    public void nextPrintSimActionNumber(){
        System.out.println((int)currentProcess[1]);
    }

    /**
    *       nextSimActionNumber returns the action number of the next 
    *   simulation process.
    */
    public int nextSimActionNumber(){
        return (int)currentProcess[1];
    }
    
    /**
    *       nextSimAction prints out the action of the next 
    *   simulation process.
    */
    public void nextSimAction(){
        printProcess(currentProcess);
    }

    private void printProcess(double[] action){
        int procId = (int)action[1];
        System.out.println("SimProcessParser | At t: "+action[0]+"   action number: "+procId);
        switch (procId) {
            case 1: System.out.println("SimProcessParser | Flip Field");break;
            case 2: System.out.println("SimProcessParser | Save Config");break;
            case 3: System.out.println("SimProcessParser | Flip Field  and Save Config");break;
            case 4: System.out.println("SimProcessParser | FlipField and Save Seed");break;
            case 5: System.out.println("SimProcessParser | Flip Field and Change Seed");break;
            case 6: System.out.println("SimProcessParser | Change Seed and Save Config");break;
            case 7: System.out.println("SimProcessParser | Flip Field ,Change Seed , and Save Config");break;
            case 8: System.out.println("SimProcessParser | Save Config and Seed");break;
            case 9: System.out.println("SimProcessParser | Save Seed");break;
            case 10: System.out.println("SimProcessParser | Change Seed");break;
            case 11: System.out.println("SimProcessParser | Change Temperature to "+ action[2]);break;
            case 12: System.out.println("SimProcessParser | Flip Field and Change Temperature to "+ action[2]);break;
            case 13: System.out.println("SimProcessParser | Change Temperature and hField to "+ action[2]);break;
            case 14: System.out.println("SimProcessParser | Change Temperature and hField back");break;
            case 15: System.out.println("SimProcessParser | Change hField to "+ action[2]);break;
            case 16: System.out.println("SimProcessParser | Change Temp back");break;
            case 17: System.out.println("SimProcessParser | Change Field back");break;
            case 18: System.out.println("SimProcessParser | Flip Field Change Temp back");break;
            default: System.out.println("SimProcessParser | Do Nothing");break;
        }
    }
    
    /**
    *       printAvailableProcessesKeywords prints available processes keywords as defined
    *  in the simprocessparser class.
    */
    public void printAvailableProcessesKeywords(){
        System.out.println("SimProcessParser | FlipField");
        System.out.println("SimProcessParser | SaveConfig");
        System.out.println("SimProcessParser | FlipFieldSaveConfig");
        System.out.println("SimProcessParser | SaveConfigFlipField");
        System.out.println("SimProcessParser | FlipFieldSaveSeed");
        System.out.println("SimProcessParser | SaveSeedFlipField");
        System.out.println("SimProcessParser | FlipFieldChangeSeed");
        System.out.println("SimProcessParser | FlipFieldChangeSeed");

        System.out.println("SimProcessParser | ChangeSeedSaveConfig");
        System.out.println("SimProcessParser | SaveConfigChangeSeed");

        System.out.println("SimProcessParser | FlipFieldChangeSeedSaveConfig");
        System.out.println("SimProcessParser | SaveConfigFlipFieldChangeSeed");
        System.out.println("SimProcessParser | ChangeSeedFlipFieldSaveConfig");
        System.out.println("SimProcessParser | ChangeSeedSaveConfigFlipField");
        System.out.println("SimProcessParser | FlipFieldSaveConfigChangeSeed");
        System.out.println("SimProcessParser | SaveConfigChangeSeedFlipField");

        System.out.println("SimProcessParser | DoNothing");
        System.out.println("SimProcessParser | SaveConfigSaveSeed");
        System.out.println("SimProcessParser | SaveSeedSaveConfig");

        System.out.println("SimProcessParser | SaveSeed");
        System.out.println("SimProcessParser | ChangeSeed");
        System.out.println("SimProcessParser | ChangeT");
        System.out.println("SimProcessParser | sameTandH");
        System.out.println("SimProcessParser | sameHandT");
        System.out.println("SimProcessParser | changeH");
        System.out.println("SimProcessParser | FlipFieldChangeBackTandH"); 
        System.out.println("SimProcessParser | ChangeTFlipField");
        System.out.println("SimProcessParser | FlipFieldChangeT");
        System.out.println("SimProcessParser | FlipFieldChangeBackT");
        System.out.println("SimProcessParser | ChangeBackTemp");
        System.out.println("SimProcessParser | ChangeBackField");
    }

    /**
    *       timeToFlipField returns the amount of time until the field was first
    *   flipped.
    */
    public int timeToFlipField(){
        int t;int tflip=0; int nextAct;
        LinkedList<double []> acts = getListOfActions();
        for (int u = 0; u< acts.size();u++){
            double[] actin = acts.get(u);
            t = (int)actin[0];
            nextAct = (int)actin[1];
            
            if(tflip==0){
                switch (nextAct) {
                case 1: tflip = t;break;
                case 2: break;
                case 3: tflip = t;break;
                case 4: tflip = t;break;
                case 5: tflip = t;break;
                case 6: break;
                case 7: tflip = t;break;
                case 8: break;
                case 9: break;
                case 10: break;
                case 11: break;
                case 12: tflip = t;break;
                case 14: tflip = t;break;
                default: break;}
            }
        }
        return tflip;
    }
    
    // test the class
    public static void main(String[] args) {
        SimProcessParser sim = new SimProcessParser("Explorer");
        System.out.println("SimProcessParser | t: "+sim.nextSimProcessTime());
        sim.nextSimActionNumber();
        sim.nextSimAction();
        System.out.println("SimProcessParser | t: "+sim.nextSimProcessTime());
        sim.nextPrintSimActionNumber();
        sim.nextSimAction();
        System.out.println("SimProcessParser | t: "+sim.nextSimProcessTime());
        sim.nextPrintSimActionNumber();
        sim.nextSimAction();

        System.out.println("SimProcessParser | tflip: "+sim.timeToFlipField());
        System.out.println("SimProcessParser | Done!");
    }
}
