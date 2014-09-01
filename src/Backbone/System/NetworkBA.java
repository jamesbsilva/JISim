/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package Backbone.System;
/**
 * 
 *   @(#) Network
 */  
import Backbone.Util.DataSaver;
import Backbone.Util.DirAndFileStructure;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;
import java.util.Scanner;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicIntegerArray;

/**
 *      Network interface to be implemented with all compatible lattice 
 *  classes. 
 * 
 * 
 * <br>
 * @author      James Silva <jbsilva @ bu.edu>                 
 * @since       2012-11    
 */
public class NetworkBA implements Network{
    private AtomicIntegerArray nodeIsingVals;   // lattice want to implement using just one index
    private HashMap<Integer,ArrayList<Integer>> networkBonds;   // lattice want to implement using just one index
    private HashMap<Integer,Boolean> networkFixed;   // lattice want to implement using just one index
    private BufferedImage networkImg;
    private Graphics2D g;
    private int scale = 5;
    private int ImageTextOffset = 90;
    private int networkEdges = 0;
    private int initM = 2;
    private int addedM = 1;
    private int N;
    private int instId;
    private int s0=1;
    private int nfixed = 0;
    private double bondProb;
    private AtomicInteger magnetization;
    private boolean initialized = false;
    private boolean useHeter = false;
    private boolean makingVideo = false;
    private String fixedFname ="";
    private int seed = 0;
    private int seedBonds = 0;
    private int seedVals = 0;
    private String degreeOutFname ="Degree/degreeDataBA-";
    private Random Ran = new Random();
   
    public NetworkBA(){
        seed = (int)(Math.random()*Integer.MAX_VALUE);
        setNetworkSeed(seed);
        initialize();
        degreeOutFname +=".txt"; 
    }
    
    public NetworkBA(int nodes,int m0,int edgesM, String fixPost, int seedin, boolean heter){
        seed = seedin;
        setNetworkSeed(seed);
        fixedFname = fixPost;
        initialize(nodes,m0,edgesM,heter);
        degreeOutFname += "N-"+nodes+"-m0"+m0+"-mAdd-"+edgesM+"-seed-"+seed+".txt";
    }

    public NetworkBA(int nodes,int m0,int edgesM, String fixPost, boolean heter){
        seed = (int)(Math.random()*Integer.MAX_VALUE);
        setNetworkSeed(seed);
        degreeOutFname += "N-"+nodes+"-m0"+m0+"-mAdd-"+edgesM+"-seed-"+seed+".txt";
        fixedFname = fixPost;
        initialize(nodes,m0,edgesM,heter);
    }
    
    public void setNetworkSeed(int nseed){
        seed= nseed;
        Ran.setSeed(nseed);
        seedBonds = Ran.nextInt();
        seedVals = Ran.nextInt();
    }
    public int getNetworkSeed(){return seed;}
    public HashMap<Integer,ArrayList<Integer>> getNetworkLinks(){return networkBonds;}
    public HashMap<Integer,Integer> getNetworkValues(){
        HashMap<Integer,Integer> vals = new HashMap<Integer,Integer>();
        for(int u = 0; u < N;u++){
            vals.put(u, getValue(u));
        }
        return vals;
    }
    public AtomicIntegerArray getNetworkValuesArr(){return nodeIsingVals;}
    public HashMap<Integer,Boolean> getNetworkFixed(){return networkFixed;}
   
    private void initialize(){
        if(networkBonds ==  null){networkBonds = new HashMap<Integer,ArrayList<Integer>>();}
        if(nodeIsingVals ==  null){nodeIsingVals = new AtomicIntegerArray(N);}
        connectNetwork();
        sortNetworkByDegree();
        initialize(s0);
    }

    private void initialize(int nodes,int m0,int edgesM, boolean heter){
        useHeter = heter;
        initM =m0;
        addedM = edgesM;
        N = nodes;
        if(networkBonds ==  null){networkBonds = new HashMap<Integer,ArrayList<Integer>>();}
        if(nodeIsingVals ==  null){nodeIsingVals = new AtomicIntegerArray(N);}
        connectNetwork();
        sortNetworkByDegree();
        initialize(s0);
    }
    
    /**
    *         setValue which updates the node spin value in the lattice and the magnetization
    *       with update of time in image layer .
    * 
    *  @param i - i coordinate
    *  @param s - new spin value
    *  @param t - time 
    */ 
    public void setValue(int i,int s, int t){
        if (this.getValue(i) == s) {
            return;
        }
        // Magnetization changes differently in diluted site
        if (this.getValue(i) != 0) {
            magnetization.set(magnetization.get() + 2 * s);
        } else {
            magnetization.set(magnetization.get() + s);
        }
        if(makingVideo){
            updateImg(i, s, t);
        }
        // Change the spin at site    
        nodeIsingVals.set(i, s);
    }
         
    /**
    *         getValue which gets the spin value in the lattice. 2d case
    * 
    *  @param i - i coordinate
    *  @param j - j coordinate
    */ 
    public int getValue(int i){
        return nodeIsingVals.get(i);
    }
    
    /**
    *         getN gives the size of the lattice
    */ 
    public int getN(){return N;}
    /**
    *         getMagnetization gives the magnetization of the lattice
    */ 
    public int getMagnetization(){return magnetization.get();}
    public int getTotalLinks(){return networkEdges;}
     /**
    *         setFixedNetworkValues should set all fixed values in the lattice.
    */ 
    public void setFixedNetworkValues(){
        networkFixed.clear();
        
        DirAndFileStructure dir = new DirAndFileStructure();
        String fileName  = "fixed" + fixedFname + ".txt";
        fileName = dir.getFixedDirectory() + fileName;        
        //Scanner scanner = new Scanner(fileName);
        Scanner scanner;
        int sum = 0;
        nfixed = 0;
	try {
            scanner = new Scanner(new File(fileName));
	    while(scanner.hasNextLine() && scanner.hasNextInt()) {
                Scanner scanner2 = new Scanner(scanner.nextLine());
                ArrayList<Integer> bonds = new ArrayList<Integer>();
                int i = scanner2.nextInt();
                int val = scanner2.nextInt();
                if(networkFixed.get(i) != null){networkFixed.remove(i);}
                networkFixed.put(i, true);
                if(nodeIsingVals.get(i)!=val){
                    sum += (val-nodeIsingVals.get(i));
                }
                nodeIsingVals.set(i, val);
                nfixed++;
	    }
	} catch (FileNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
	}
        magnetization.set(sum+magnetization.get());
        System.out.println("NetworkBA | Nodes fixed in network : "+nfixed+"    M: "+magnetization.get());
    }
    
    /**
    *         initialize sets the lattice values to their initial values which is given as input.
    *
    *
    *  @param s - Initial lattice value; -2 for random values
    *
    */ 
    public void initialize(int s){
        if(nodeIsingVals ==  null){nodeIsingVals = new AtomicIntegerArray(N);}
        magnetization = new AtomicInteger();
        magnetization.set(0);
        if(networkBonds ==  null){networkBonds = new HashMap<Integer,ArrayList<Integer>>();}
        networkFixed = new HashMap<Integer,Boolean>();
        int spin = s;
        int sum = 0;
        Ran.setSeed(seedVals);
        for(int u = 0;u < N;u++){
            if (s == (-2)) {
                spin = getRandom(0, 1);
                if (spin == 0) {
    
                    spin = -1;
                }
            }
            sum += spin;
            nodeIsingVals.set(u, spin);
            networkFixed.put(u, false);
        }    
        magnetization.set(sum);
        if(useHeter){setFixedNetworkValues();}
        initializeImg();
        spinDraw();
    }
 
    /**
    *         isThisFixed should return true if the lattice coordinate is fixed.
    * 
    *  @param i - i coordinate
    */ 
    public boolean isThisFixed(int i){
        if(networkFixed.get(i) == null){return false;}
        return networkFixed.get(i);
    }
    
    /**
    *         getNeighSum should return the sum of the spins within the interaction
    *   range centered by the coordinates given. 2d case
    * 
    *  @param i - i coordinate
    */ 
    public int getNeighSum(int i){
        int sum = 0;
        ArrayList<Integer> bonds = networkBonds.get(i);
        for(int u = 0;u < bonds.size();u++){
            sum += nodeIsingVals.get(bonds.get(u));
        }
        return sum;
    }
    
    /**
    *         setInitialConfig should set all spins in the lattice to the values
    *   given by a file in the Config directory.
    *   
    *   @param t - time of config to set lattice to
    *   @param run - run of time to set lattice to
    *   -*@param  post - postfix or filename of lattice file
    */
    public void setInitialConfig(int t,int run,String post){
        DirAndFileStructure dir = new DirAndFileStructure();
        DataSaver dsave = new DataSaver(instId);
        String fileName = fixedFname;
        fileName = dir.getConfigDirectory(instId) + dsave.getNetworkTempFilename(run, t, post);
        //Scanner scanner = new Scanner(fileName);
        Scanner scanner;
        int sum = 0;
        nfixed = 0;
	try {
		scanner = new Scanner(new File(fileName));
	    while(scanner.hasNextLine()) {
                Scanner scanner2 = new Scanner(scanner.nextLine());
                ArrayList<Integer> bonds = new ArrayList<Integer>();
                int i = scanner2.nextInt();
                int val = scanner2.nextInt();
                sum += val;
                int fix = scanner2.nextInt();
                while(scanner2.hasNextInt()){
                    bonds.add(scanner2.nextInt());
                }
                if(fix == 1){networkFixed.put(i, true);nfixed++;
                }else{networkFixed.put(i, false);}
                networkBonds.remove(i);
                networkBonds.put(i,bonds);
                nodeIsingVals.set(i, val);
	    }
	} catch (FileNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
	}
        magnetization.set(sum);
    }
    
    public void connectNetwork(){
        resetBonds();
        Ran.setSeed(seedBonds);
        for(int u=1; u < N;u++){
            if(u < initM){
                int connect = (int)(Ran.nextDouble()*u);
                makeBondIJ(u,connect);
                makeBondIJ(connect, u);
                networkEdges++;
            }else{
                for(int v = 0;v<addedM;v++){
                    boolean bondMade = false;
                    while(!bondMade){
                        int connect = (int)(Ran.nextDouble()*u);
                        calcBondProb(connect);
                        if(Ran.nextDouble() < bondProb){
                            bondMade = true;
                            makeBondIJ(u,connect);
                            makeBondIJ(connect,u);
                            networkEdges++;}
                    }
            }}
        }
    }
    
    private void resetBonds(){
        networkEdges = 0;
        networkBonds.clear();
        for(int u=0; u < N;u++){
            ArrayList<Integer> bonds = new ArrayList<Integer>();
            networkBonds.put(u,bonds);
        }
    }
        
    private void calcBondProb(int u){
        int uLinks = (networkBonds.get(u)).size();
        bondProb = (double)(uLinks)/(double)(networkEdges);
    }
    
    private void makeBondIJ(int i, int j){
        ArrayList<Integer> bondsI = networkBonds.get(i);
        networkBonds.remove(i);
        bondsI.add(j);
        networkBonds.put(i, bondsI);
    }
    
    /**
    *         getInstId returns the identification number of the instance of this
    *  program running. Useful for allowing each instance a directory to work in.
    */ 
    public int getInstId(){return instId;}
    /**
    *         getNFixed gives the amount of fixed spins in the lattice
    */ 
    public int getNFixed(){return nfixed;}

    public void printLinks(int u){
        ArrayList<Integer> bonds = networkBonds.get(u);
        for(int k = 0;k < bonds.size();k++){
            System.out.println("Node: "+u+"   Bonded: "+bonds.get(k));
        }
    }
    
    public void printDegree(int u){
        ArrayList<Integer> bonds = networkBonds.get(u);
        System.out.println("Node: "+u+"   Degree: "+bonds.size());
        //System.out.println(+bonds.size());
    }
    
    public void printDegrees(){
        for(int u=0;u<networkBonds.size();u++){
            printDegree(u);
        }
    }
    
    public void printDegreeFrequency(){
        ArrayList<Integer> freq = new ArrayList<Integer>();
        for(int u = 1; u < networkBonds.get(0).size();u++){
            int count = 0;
            for(int v = 1; v < networkBonds.size();v++){
                if(networkBonds.get(v).size() == u){count++;}
            }
            freq.add(count);
        }
        for(int u = 0; u < networkBonds.get(0).size()-1;u++){
            System.out.println("node: "+u+"     freq: "+freq.get(u));
        }
    }

    public void saveDegreeFrequency(){
        ArrayList<Integer> freq = new ArrayList<Integer>();
        for(int u = 1; u < networkBonds.get(0).size();u++){
            int count = 0;
            for(int v = 1; v < networkBonds.size();v++){
                if(networkBonds.get(v).size() == u){count++;}
            }
            freq.add(count);
        }
        DataSaver dsave = new DataSaver(instId);
        for(int u = 0; u < networkBonds.get(0).size()-1;u++){
            double[] data = new double[2];
            data[0] = u+1;data[1]=freq.get(u);
            dsave.saveDoubleArrayData(data, degreeOutFname);
        }
    }
    
    public int getRandom(int min, int max) {
        return (min + (int) (Ran.nextDouble() * ((max - min) + 1)));
    }
    
    private void sortNetworkByDegree(){
        boolean sorted=false;
        while(!sorted){
            sorted = true;
            for(int u = 0;u < networkBonds.size()-1;u++){
                if(networkBonds.get(u+1).size() > networkBonds.get(u).size()){
                    exchangeNodes(u+1, u);
                    sorted=false;
                }
            }
        }
    }
    
    private void exchangeNodes(int i, int j){
        ArrayList<Integer> bondsI = networkBonds.get(i);
        ArrayList<Integer> bondsJ = networkBonds.get(j);
        int nodeValI = nodeIsingVals.get(i);
        int nodeValJ = nodeIsingVals.get(j);
        nodeIsingVals.set(i, nodeValJ);
        nodeIsingVals.set(j, nodeValI);
        networkBonds.put(j, bondsI);
        networkBonds.put(i, bondsJ);
    }
    
    @Override
    public BufferedImage getSystemImg() {
        return networkImg;
    }

    
    private void updateImg(int i, int spin,int t){
        //erase old text
        int L = (int)(Math.ceil(Math.sqrt(N)));
        g.setColor(Color.BLACK);
        g.fillRect(((int)L*scale*6/10),0,((int)L*scale*4/10),ImageTextOffset);
        g.setColor(Color.WHITE);
        if(t!=0){
            Font font = new Font("Courier New", Font.BOLD, 72);
            g.setFont(font);
            // 70 is the x coordinate of text
            g.drawString("t = "+ t, ((int) L*scale*6/10), ((int) ImageTextOffset*7/10));
        }
        spin = getValue(i);
        if(spin==1){g.setColor(Color.MAGENTA);}else 
            if(spin==(-1)){g.setColor(Color.CYAN);}else
                    if(spin==(0)){g.setColor(Color.WHITE);}
        int x = i%L;
        int y = ((int)((double)i/(double)L))%L;
        g.fillRect(x*scale, y*scale+ImageTextOffset, scale, scale);
    }

    private void spinDraw(){
        int spin;
        int L = (int)(Math.ceil(Math.sqrt(N)));
        for(int i = 0;i < N;i++){ 
            if(isThisFixed(i) ){
                    // Paint different color for fixed spins
                    spin = getValue(i);
                    if(spin==1){g.setColor(Color.GRAY);} 
                    if(spin==(-1)){g.setColor(Color.RED);}
                    if(spin==(0)){g.setColor(Color.MAGENTA);}
            }else{
            spin = getValue(i);
            if(spin==1){g.setColor(Color.YELLOW);}
                if(spin==(-1)){g.setColor(Color.CYAN);}
                        if(spin==(0)){g.setColor(Color.WHITE);}}
            int x = i%L;
            int y = ((int)((double)i/(double)L))%L;
            g.fillRect(x*scale, y*scale+ImageTextOffset, scale, scale);		 
        }
    }
    private void initializeImg(){
        int L = (int)(Math.ceil(Math.sqrt(N)));
        int sizeImage = scale*L;
        while(sizeImage < 600){scale = scale*2;sizeImage = scale*L;}
        networkImg = new BufferedImage(sizeImage, sizeImage+ImageTextOffset, BufferedImage.TYPE_INT_RGB);  
        g = networkImg.createGraphics();
    }
    
    
    @Override
    public void makeVideo() {
        makingVideo = true;
    }
    
    // test the class
    public static void main(String[] args) {
        NetworkBA network = new NetworkBA(200,10,1,"",false);
        System.out.println("Total Links: "+network.networkEdges);
        //network.saveDegreeFrequency();
        //network.printDegreeFrequency();
        //network.printDegrees();
    }
}


