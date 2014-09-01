package Backbone.Util;

/**
* 
*    @(#) ParameterBank
*/
import Backbone.System.SimpleLattice;
import Triggers.*;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;

/** 
*   A parser for defining the parameters for the simulation.
*  <br>
* 
*  @param fname - postfix for ParameterBank file.
* 
* @author      James Silva <jbsilva @ bu.edu>                 
* @since       2011-09    
*/
public class ParameterBank {
	public double temperature;
	public double jInteraction;
	public double hField;
	public boolean useLongRange;
	public boolean useDilution;
	public boolean useHeter;
	public int Q;
	public int L;
	public int N;
	public int D;
	public int run;
	public int Geo;
        public int s; //initial spin (-2) is random
        public String mcalgo;
        public String trigger;
        public double triggerParam1;
        public double triggerParam2;
        public int R = 1;
        public int Granularity = 0;
        public boolean CompleteInfoGranular = false;
        public String postFix;
        private ArrayList<String> file=null;
        private HashMap<String,ArrayList<String>> aliases;
        private ParsingHelper parser = new ParsingHelper();
        private DirAndFileStructure dir;
        
    public ParameterBank(){
    	this("");
    }
 
    /**
    *  creates map of aliases for parameter keywords
    */
    private void createAliasesMap(){
        aliases = new HashMap<String,ArrayList<String>>();
        ArrayList lengthArr = new ArrayList<String>();
        lengthArr.add("Length");
        lengthArr.add("length");
        lengthArr.add("L");
        lengthArr.add("l");
        ArrayList sArr = new ArrayList<String>();
        sArr.add("s0");
        sArr.add("initialSpin");
        sArr.add("spinInitial");
        sArr.add("spin");
        ArrayList algoArr = new ArrayList<String>();
        algoArr.add("algo");
        algoArr.add("algorithm");
        algoArr.add("mcalgorithm");
        algoArr.add("mcalgo");
        algoArr.add("mc");
        algoArr.add("MC");
        ArrayList interactionArr = new ArrayList<String>();
        interactionArr.add("interaction");
        interactionArr.add("jInteraction");
        interactionArr.add("j");
        interactionArr.add("J");
        ArrayList rangeArr = new ArrayList<String>();
        rangeArr.add("R");
        rangeArr.add("r");
        rangeArr.add("range");
        rangeArr.add("Range");    
        ArrayList dimArr = new ArrayList<String>();
        dimArr.add("dim");
        dimArr.add("dimension");
        dimArr.add("Dimension");
        dimArr.add("D");
        ArrayList dilutionArr = new ArrayList<String>();
        dilutionArr.add("dilution");
        dilutionArr.add("diluted");
        dilutionArr.add("useDilution");
        dilutionArr.add("Dilution");
        dilutionArr.add("Diluted");
        ArrayList heterArr = new ArrayList<String>();
        heterArr.add("heterogenous");
        heterArr.add("Heterogenous");
        heterArr.add("homogenous");
        heterArr.add("Homogenous");
        heterArr.add("hetero");
        ArrayList granularArr = new ArrayList<String>();
        granularArr.add("granularity");
        granularArr.add("Granularity");
        ArrayList fieldArr = new ArrayList<String>();
        fieldArr.add("hField");
        fieldArr.add("Field");
        fieldArr.add("field");
        fieldArr.add("hfield");
        fieldArr.add("h");
        ArrayList tempArr = new ArrayList<String>();
        tempArr.add("temperature");
        tempArr.add("Temperature");
        tempArr.add("temp");
        tempArr.add("Temp");
        tempArr.add("T");
        ArrayList geoArr = new ArrayList<String>();
        geoArr.add("geometry");
        geoArr.add("Geometry");
        geoArr.add("geo");
        geoArr.add("Geo");
        ArrayList trigArr = new ArrayList<String>();
        trigArr.add("trigger");
        trigArr.add("Trigger");
        ArrayList trigparamArr = new ArrayList<String>();
        trigparamArr.add("triggerparam");
        trigparamArr.add("Triggerparam");
        trigparamArr.add("triggerParam");
        trigparamArr.add("TriggerParam");
        trigparamArr.add("TriggerParameter");
        trigparamArr.add("trigger");
        trigparamArr.add("Trigger");
        
        aliases.put("L",lengthArr);
        aliases.put("Geo",geoArr);
        aliases.put("D",dimArr);
        aliases.put("temperature",tempArr);
        aliases.put("trigger",trigArr);
        aliases.put("triggerParam",trigparamArr);
        aliases.put("hField", fieldArr);
        aliases.put("jInteraction", interactionArr);
        aliases.put("s", sArr);
        aliases.put("R", rangeArr);
        aliases.put("mcalgo", algoArr);
        aliases.put("useHeter", heterArr);
        aliases.put("useDilution", dilutionArr);
        aliases.put("Granularity", granularArr);
    }
    
    
    private int getTriggerNParams(String trig){
        if(trig.equalsIgnoreCase("simplevaluefit")
                ||trig.equalsIgnoreCase("simplevalfit")){
            return 2;
        }else{
            return 1;
        }
    }
    
    
    /*
    *       getProperAlias gets the proper alias to search for in parameter file
    *   by testing all possible aliases.
    * 
    * @param param - parameter to search for
    * @param filename - filename to be parsed
    * 
    */
    private String getProperAlias(String param, String filename){
        String alias;
        int maxAttempts = 15;
        int attempt=0;
        alias = getProperAlias(param, filename,0);
        while(alias.equalsIgnoreCase("noGoodString") && attempt <= maxAttempts){
            attempt++;alias = getProperAlias(param,filename,0);
        }   
        
        if(alias.equalsIgnoreCase("noGoodString")){
            System.err.println("ERROR PARAMETERS: "+param+"   NOT FOUND");
        }
        
        return alias;
    }
    
    /*
    *       getProperAlias gets the proper alias to search for in parameter file
    *   by testing all possible aliases.
    * 
    * @param param - parameter to search for
    * @param filename - filename to be parsed
    * 
    */
    private String getProperAlias(String param, ArrayList<String> filename){
        String alias;
        int maxAttempts = 15;
        int attempt=0;
        alias = getProperAlias(param, filename,0);
        while(alias.equalsIgnoreCase("noGoodString") && attempt <= maxAttempts){
            attempt++;alias = getProperAlias(param,filename,0);
        }   
        return alias;
    }
    
    /*
    *       getProperAlias gets the proper alias to search for in parameter file
    *   by testing all possible aliases.
    * 
    * @param param - parameter to search for
    * @param filename - filename to be parsed
    * @param postNum - number to append to parameter (used to disntinguish many params)
    */    
    private String getProperAlias(String param, String filename, int postNum){
        ArrayList<String> alias = aliases.get(param);
        String goodString="noGoodString";
        String post = "";
        if(postNum > 0){post = post+postNum;}
        
        for(int u = 0; u < alias.size();u++){
            if(goodString.equals("noGoodString") &&
                parser.checkParameterExist(filename, alias.get(u)+post)){goodString = alias.get(u)+post;}
        }
        //System.out.println("ParameterBank | Good String: "+goodString);
        return goodString;
    }

    /*
    *       getProperAlias gets the proper alias to search for in parameter file
    *   by testing all possible aliases.
    * 
    * @param param - parameter to search for
    * @param filename - filename to be parsed
    * @param postNum - number to append to parameter (used to disntinguish many params)
    */    
    private String getProperAlias(String param, ArrayList<String> filename, int postNum){
        ArrayList<String> alias = aliases.get(param);
        String goodString="noGoodString";
        String post = "";
        if(postNum > 0){post = post+postNum;}
        
        for(int u = 0; u < alias.size();u++){
            if(goodString.equals("noGoodString") &&
                parser.checkParameterExist(filename, alias.get(u)+post)){goodString = alias.get(u)+post;}
        }
        //System.out.println("ParameterBank | Good String: "+goodString);
        return goodString;
    }
    
    /**
    * 
    * @param name - the post of the parameter file name 
    */
    public ParameterBank(String name){
        this(null,name);
    }
    /**
    * 
    * @param name - the post of the parameter file name 
    */
    public ParameterBank(ParameterBank param,String name){
        if(param == null){
            setupNewParamBank(name);
        }else{
            copyParamBank(param);
        }
    }
    
    private void copyParamBank(ParameterBank param){        
        dir  = new DirAndFileStructure();
        createAliasesMap();
        
        L = param.L;
        Granularity = param.Granularity;
        postFix = param.postFix;
        R = param.R;
        Geo = param.Geo;
        N = param.N;
        Q = param.Q;
        hField = param.hField;
        jInteraction = param.jInteraction;
        mcalgo = param.mcalgo;
        s = param.s;
        temperature = param.temperature;
        trigger = param.trigger;
        triggerParam1 = param.triggerParam1;
        triggerParam2 = param.triggerParam2;
        useDilution = param.useDilution;
        useHeter = param.useHeter;
        useLongRange = param.useLongRange;
        CompleteInfoGranular = param.CompleteInfoGranular;
        run = param.run;
    }
    
    private void setupNewParamBank(String name){
        dir  = new DirAndFileStructure();
        createAliasesMap();
        postFix =name;
        String fname = "parameters"+name+".txt";
        fname = dir.getSettingsDirectory()+fname;
        File f = new File(fname);
        if(f.exists()){}else{
            System.out.println(fname+"  does not exist. Using default.");
            fname = dir.getSettingsDirectory()+"parameters.txt";
        }
        file = parser.parseFileIntoStringList(fname);

        L = parser.parseFileInt(file, getProperAlias("L", file), 0);
        Granularity = L;
        if(!getProperAlias("Granularity",file).equals("noGoodString")){
            Granularity = parser.parseFileInt(file, getProperAlias("Granularity", file), 0);
        }
        s = parser.parseFileInt(file, getProperAlias("s", file), 0);
        Geo = parser.parseFileInt(file, getProperAlias("Geo", file), 0);
        D = parser.parseFileInt(file, getProperAlias("D", file), 0);

        useDilution = parser.parseFileBoolean(file, getProperAlias("useDilution", file), 0);
        useHeter = parser.parseFileBoolean(file, getProperAlias("useHeter", file), 0);
        if(getProperAlias("useHeter", file).contains("homo") || 
                getProperAlias("useHeter", file).contains("Homo")){
            useHeter = !useHeter;
        }

        useLongRange = false;R=1;
        if(!getProperAlias("R",file).equals("noGoodString")){
            R = parser.parseFileInt(file, getProperAlias("R", file), 0);
            if(R > 0){useLongRange = true;}
        }

        hField = parser.parseFileDouble(file, getProperAlias("hField", file), 0);
        jInteraction = parser.parseFileDouble(file, getProperAlias("jInteraction", file), 0);
        temperature = parser.parseFileDouble(file, getProperAlias("temperature", file), 0);

        mcalgo =  parser.parseFileFirstString(file, getProperAlias("mcalgo", file));
        trigger =  parser.parseFileFirstString(file, getProperAlias("trigger", file));            
        ArrayList<Double> params = parser.parseIntoDblArray(file, getProperAlias("trigger", file));
        triggerParam1=-1;
        
        if(params.size() > 0){
            int u = params.size();
            int i = 0;
            while(u > 0){
                if(i == 0){triggerParam1 = params.get(i);}
                if(i == 1){triggerParam2 = params.get(i);}
                u--;
            }
        }

        String alias;
        if(!getProperAlias("triggerParam", file,1).equals("")){
            alias = getProperAlias("triggerParam", file,1);
            while(alias.equalsIgnoreCase("noGoodString")){
                System.out.println("attempting to get alias again");
                alias = getProperAlias("triggerParam", file,1);}
            triggerParam1 =  parser.parseFileDouble(file,alias,0);            
        }
        
        if(trigger.equalsIgnoreCase("simplevaluefit")){
            if(!getProperAlias("triggerParam", file,2).equals("") ){
                System.out.println("attempting to get alias again");
                alias = getProperAlias("triggerParam", file,2);
                while(alias.equalsIgnoreCase("noGoodString")){alias = getProperAlias("triggerParam", file,2);}
                triggerParam2 =  parser.parseFileDouble(file, alias,0);            
            }
        }

        //SimpleLattice lat = new SimpleLattice(L,R,D,Geo,s);
        N = (int)Math.pow(L, D);
        Q = (int)Math.pow(2*R+1, D)-1;
    }
    
    /**
    *       setN sets the size of Lattice
    * 
    * @param num - size of lattice
    */
    public void setN(int num){N = num;}

    /**
    *       setQ sets the amount of spins in range for this parameter bank
    * 
    * @param qnum - number of spins in range
    */
    public void setQ(int qnum){Q = qnum;}	

    /**
    *      changeTemp writes the new temperature in the parameter file.
    * 
    * 
    * @param tempnow - new temperature
    * @param name - postfix of filename of parameter file
    */
    public void changeTemp(double tempnow,String name) {
        dir  = new DirAndFileStructure();
        postFix =name;
        String fname = "parameters"+name+".txt";
        fname = dir.getSettingsDirectory()+fname;

        String edit = ""+tempnow;
        parser.editFile(fname, getProperAlias("temperature", fname), 0, edit, 1);
        temperature = tempnow;
    }

    /**
    *      changeJinteraction changes the interaction in the parameter file.
    * 
    * 
    * @param tempnow - new temperature
    * @param name - postfix of filename of parameter file
    */
    public void changeJinteraction(double jnow,String name) {
        dir  = new DirAndFileStructure();
        postFix =name;
        String fname = "parameters"+name+".txt";
        fname = dir.getSettingsDirectory()+fname;

        String edit = ""+jnow;
        parser.editFile(fname, getProperAlias("jInteraction", fname), 0, edit, 1);
        jInteraction = jnow;
    }
    
    /**
    *      changeLength writes the new length in the parameter file.
    * 
    * 
    * @param lNow - new length
    * @param name - postfix of filename of parameter file
    */
    public void changeLength(int lNow,String name) {
        dir  = new DirAndFileStructure();
        postFix =name;
        String fname = "parameters"+name+".txt";
        fname = dir.getSettingsDirectory()+fname;

        String edit = ""+lNow;
        parser.editFile(fname, getProperAlias("L", fname), 0, edit, 1);
        L = lNow;
    }

    /**
    *      changeRange writes the new range in the parameter file.
    * 
    * 
    * @param rNow - new range
    * @param name - postfix of filename of parameter file
    */
    public void changeRange(int rNow,String name) {
        dir  = new DirAndFileStructure();
        postFix =name;
        String fname = "parameters"+name+".txt";
        fname = dir.getSettingsDirectory()+fname;

        String edit = ""+rNow;
        parser.editFile(fname, getProperAlias("R", fname), 0, edit, 1);
        R = rNow;
    }
    
    /**
    *      changeField writes the new magnetic field in the parameter file.
    * 
    * 
    * @param hnow - new magnetic field
    * @param name - postfix of filename of parameter file
    */
    public void changeField(double hnow,String name) {
            dir  = new DirAndFileStructure();
            postFix =name;
            String fname = "parameters"+name+".txt";
            fname = dir.getSettingsDirectory()+fname;
            String edit = ""+hnow;
            parser.editFile(fname, getProperAlias("hField", fname), 0, edit, 1);     
            hField = hnow;
    }

    
    /**
     *      printParameters prints out the relevant parameter into the console
     */  
    public void printParameters(){
        System.out.println("ParameterBank | ________________________________________________________");
        System.out.println("ParameterBank | L: "+L);
        System.out.println("ParameterBank | jInteraction: "+jInteraction);
        System.out.println("ParameterBank | granularity: "+Granularity);
        System.out.println("ParameterBank | s: "+s);
        System.out.println("ParameterBank | temp: "+temperature);
        System.out.println("ParameterBank | hField: "+hField);
        System.out.println("ParameterBank | dim: "+D);
        System.out.println("ParameterBank | Geo: "+Geo);
        System.out.println("ParameterBank | Range: "+R);
        System.out.println("ParameterBank | mc algo :" +mcalgo);
        System.out.println("ParameterBank | trigger :" +trigger);
        System.out.println("ParameterBank | using   Heter: "+useHeter+"    Dilution:"+useDilution+"     Long Range:"+useLongRange);
        System.out.println("ParameterBank | ________________________________________________________");		
    }

    /**
    *       printAvailableParamOptions prints available parameter options as defined
    *  in the ParameterBank class.
    */
    public void printAvailableParamOptions(){
        //iterating over keys only
        for (String key : aliases.keySet()) {
            System.out.println("ParameterBank | Parameter Option : " + key);
        }
    
    }
    
    public Trigger setProperTrigger(ParameterBank param, Trigger trig, String SimPost,boolean output){
        Trigger trigger = null;
        if(param.trigger.equalsIgnoreCase("deviation")){
            trigger = new
                DeviationTrigger(param.s,
                param.triggerParam1,param.triggerParam2,SimPost);
        }else if(param.trigger.equalsIgnoreCase("simplevalue") || param.trigger.equalsIgnoreCase("simple")){
            trigger = new 
                SimpleValueTrigger(param.s,
                param.triggerParam1,param.N,SimPost,param.postFix);}
        else if(param.trigger.equalsIgnoreCase("simplevaluefit") ||
            param.trigger.equalsIgnoreCase("simplefit")|| 
            param.trigger.equalsIgnoreCase("simplevalfit")){
            trigger = new
                SimpleValFitTrigger(param.s,param.triggerParam1,
                    param.N,SimPost,param.triggerParam2,output);}
        else if (param.trigger.equalsIgnoreCase("simpleEvalue") || param.trigger.equalsIgnoreCase("simpleE")
                ||param.trigger.equalsIgnoreCase("simpleEval")){
            trigger = new 
                SimpleEValueTrigger(param.s,
                param.triggerParam1,param.N,SimPost,param.postFix);
        }else{
            System.err.println("ERROR: NO APPLICABLE TRIGGER FOUND FOR : "+param.trigger);
        }
        return trigger;
    }
    
    public static void main(String[] args) {
        ParameterBank param = new ParameterBank("");
        //param.changeTemp(-8.4, "");
        //param.changeRange(0, "");
        //param.changeField(8.4, "");
        
        System.out.println("ParameterBank | L: "+param.L);
        System.out.println("ParameterBank | N: "+param.N);
        System.out.println("ParameterBank | s: "+param.s);
        System.out.println("ParameterBank | temp: "+param.temperature);
        System.out.println("ParameterBank | hField: "+param.hField);
        System.out.println("ParameterBank | dim: "+param.D);
        System.out.println("ParameterBank | Geo: "+param.Geo);
        System.out.println("ParameterBank | Q: "+param.Q);
        System.out.println("ParameterBank | R: "+param.R);
        System.out.println("ParameterBank | mc algo :" +param.mcalgo);
        System.out.println("ParameterBank | trigger :" +param.trigger);
        System.out.println("ParameterBank | using   Heter: "+param.useHeter+"    Dilution:"+param.useDilution+"     Long Range:"+param.useLongRange);
        System.out.println("ParameterBank | Done!");
    }	
}
