package JISim;

/**
* 
*    @(#) MakeFixedConfig2D
*/  
 
import Backbone.Util.DirAndFileStructure;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.Random;
import java.util.Scanner;
import java.util.Vector;

/** 
 *      Helper class for making a fixed spin configuration.
 *  <br>
 * 
 *  @param type - which type of config - box,square,random,circle,crack,cross
 *  @param s - spin for the fixed spins
 *  @param x - x coordinate of center of configuration
 *  @param y - y coordinate of center of configuration
 *  @param L - size of lattice
 *  @param L1 - first dimension of configuration , radius of circle
 *  @param L2 - second dimension of configuration , length of crack or amount 
 *              of fixed spins
 *  
 * @author      James Silva <jbsilva @ bu.edu>                 
 * @since       2011-12    
 */

public class MakeFixedConfig {    
/**
*       clearFixedConfig2D clears the fixed spin configuration file that 
* has the given post fix.
* 
* @param type - post fix for the filename of this fixed spin configuration
*/    
public static void clearFixedConfig2D(String type){    	
    DirAndFileStructure dir = new DirAndFileStructure();
    String fName = dir.getFixedDirectory()+"fixed"+type+".txt";
 
    try {
        PrintStream out = new PrintStream(new FileOutputStream(
            fName,false));
        out.println("   ");

        out.close();
    } catch (FileNotFoundException e) {
        e.printStackTrace();
    }

    System.out.println("Done with Clearing Fixed Config!");
}
 
/**
* 
*       MakeFixedConfig2D makes the fixed spin configuration and saves it to a
*   file in the fixed spin configuration directory.
* 
* @param type - type of fixed spin configuration
* @param s - spin of these fixed spins
* @param x - x coordinate of the center of the configuration
* @param y -  y coordinate of the center of the configuration
* @param L - Length of the lattice.
* @param L1 - First Parameter of the fixed configuration. For circle this is R
* @param L2 - Second parameter of the fixed configuration. For crack this is the length.
*/
public static int[][] makeFixedConfig2D(String type,int s, int x, int y,int L, int L1, int L2){
	return makeFixedConfig2D(type,"",s,x,y,L,L1,L2,(int)(Math.random()*1026475));
}

/**
* 
*       MakeFixedConfig2D makes the fixed spin configuration and saves it to a
*   file in the fixed spin configuration directory.
* 
* @param type - type of fixed spin configuration
* @param s - spin of these fixed spins
* @param x - x coordinate of the center of the configuration
* @param y -  y coordinate of the center of the configuration
* @param L - Length of the lattice.
* @param L1 - First Parameter of the fixed configuration. For circle this is R
* @param L2 - Second parameter of the fixed configuration. For crack this is the length.
*/
public static int[][] makeFixedConfig2D(String type,int s, int x, int y,int L, int L1, int L2, int seed){
	return makeFixedConfig2D(type,"",s,x,y,L,L1,L2,seed);
}

/**
* 
*       MakeFixedConfig2D makes the fixed spin configuration and saves it to a
*   file in the fixed spin configuration directory.
* 
* @param type - type of fixed spin configuration
* @param postFix - postfix for fixed spin configuration file
* @param s - spin of these fixed spins
* @param x - x coordinate of the center of the configuration
* @param y -  y coordinate of the center of the configuration
* @param L - Length of the lattice.
* @param L1 - First Parameter of the fixed configuration. For circle this is R
* @param L2 - Second parameter of the fixed configuration. For crack this is the length.
*/
public static int[][] makeFixedConfig2D(String type,String post,int s, int x, int y,int L, int L1, int L2){
	return makeFixedConfig2D(type,post,s,x,y,L,L1,L2,(int)(Math.random()*1026475));
}

/**
* 
*       makeFixedConfig2D makes the fixed spin configuration and saves it to a
*   file in the fixed spin configuration directory.
* 
* @param type - type of fixed spin configuration
* @param postFix - postfix for fixed spin configuration file
* @param s - spin of these fixed spins
* @param x - x coordinate of the center of the configuration
* @param y -  y coordinate of the center of the configuration
* @param L - Length of the lattice.
* @param L1 - First Parameter of the fixed configuration. For circle this is R
* @param L2 - Second parameter of the fixed configuration. For crack this is the length.
* @param seed - random number seed.
*/
public static int[][] makeFixedConfig2D(String type,String postFix,int s, int x, int y,int L, int L1, int L2, int seed){    
    int[][] fixed;
    Random Ran = new Random();
    Ran.setSeed(seed);
    System.out.println("MakeFixedConfig | Making fixed with config : "+type);
        
    if(type.equalsIgnoreCase("randomspin")||type.equalsIgnoreCase("crack")){
        fixed = new int[2*L2+10][3];
    }else if(type.equalsIgnoreCase("box")){
        fixed = new int[4*(L1+L2)+3][3];
    }else if (type.equalsIgnoreCase("circle")){
        L1 = (int) (L1/2.0);
        fixed = new int[(L1*L1)+3][3];
    }else if (type.equalsIgnoreCase("rectangle")){
        fixed = new int[L1*L2+10][3];
    }else if(type.equalsIgnoreCase("square")){
        L2=L1;
        fixed = new int[L1*L2+3][3];       
    }else if(type.equalsIgnoreCase("balancedrandom")){
        L2 = 2*L2;
        L1=L2;
        fixed = new int[2*L2+10][3];       
    }else{
        fixed = new int[2*(L1+L2)+10][3];
    }
    
    int nfixed=0;
    if(L1 == 0){L1=1;}
    if(L2 == 0){L2=1;}
    int u=0;
    int v=0;
    int sum = 0;
    
    // need to add something to make sure not in same
    HashMap<Integer,Integer> currFixed = new HashMap<Integer,Integer>();
    for(int i = 0;i<(L*L);i++){currFixed.put(i, 0);}
        
    if(type.equalsIgnoreCase("randomspin")){
        boolean goodSpin;
        for(int j = 0;j<L2;j++){
            goodSpin = false;
            while(!goodSpin){
                u = (int) (Ran.nextDouble()*L) ;
                v = (int) (Ran.nextDouble()*L) ;
                // if spin hasnt been fixed then continue
                if(currFixed.get(u+v*L)== 0){goodSpin=true;}    
            }
            s = (int) (Ran.nextDouble()*2);
            if(s==0){s=-1;}
            nfixed++;
            fixed[j][0]= u;
            fixed[j][1]= v;
            currFixed.put((u+v*L), 1);
            fixed[j][2]= s;
            sum += s;
        }

    }
        
        
    // Random fixing of spins where number fixed is predetermined and imposed
    if((type.contains("random") || type.contains("Random"))&&(
            type.contains("conserved")|| type.contains("Conserved"))){
        boolean goodSpin = false;
        for(int j = 0;j<L2;j++){
            goodSpin = false;
            while(!goodSpin){
                u = (int) (Ran.nextDouble()*L) ;
                v = (int) (Ran.nextDouble()*L) ;
                // if spin hasnt been fixed then continue
                if(currFixed.get(u+v*L)== 0){goodSpin=true;}    
            }
            
            nfixed++;
            fixed[j][0]= u;
            fixed[j][1]= v;
            currFixed.put((u+v*L), 1);
                if((type.contains("balanced") || type.contains("Balanced")) && ((j%2)==1)){
                    fixed[j][2]= (-1)*s;sum += (-1)*s;
                }else{                    
                    fixed[j][2]= s;sum += s;
                }
        }
    }
    
    // Random fixing of spins where number fixed isnt predetermined
    if(type.contains("random")||type.contains("Random")){
        double probAccept = ((double)L2)/((double)(L*L));
        for(int i = 0;i<L;i++){
            for(int j = 0;j<L;j++){
            if(Ran.nextDouble() < probAccept){
                nfixed++;
                fixed[nfixed][0]= i;
                fixed[nfixed][1]= j;
                currFixed.put((i+j*L), 1);
                    if((type.contains("balanced") || type.contains("Balanced"))){
                        if(Ran.nextDouble() > 0.5){
                            fixed[nfixed][2]= (-1)*s;sum += (-1)*s;
                        }
                    }else{                    
                        fixed[nfixed][2]= s;sum += s;
                    }
            }
        }}
    }
    
    if(type.contains("interactionquilt") || type.contains("InteractionQuilt")
            || type.contains("interactionQuilt") || type.contains("Interactionquilt")){
        int R = 2*L2+1;
        double p = (double)L1/(double)(L*L);
        int subNfixed = (int) (p*R*R);
        Vector<int[]> patches = new Vector<int[]>();
        int nPatches = x;
        for(int e=0;e<nPatches;e++){
        int[] fixedSub = new int[subNfixed];
            for(int k = 0;k < subNfixed;k++){
                boolean prevFix =true;
                int t = (int) (Ran.nextDouble()*(R*R));
                while(prevFix){
                    prevFix = false;
                    for(int m = 0;m < subNfixed;m++){
                        if(fixedSub[m] == t){prevFix=true;}
                    }
                    if(prevFix){t = (int) (Ran.nextDouble()*(R*R));}
                }

                fixedSub[k] = t;
            }
            patches.add(fixedSub);
        }
        int currentPatch=1;
        int patchInd =0;
        int patchMod = y;
        
        System.out.println("Making Quilt "+L1+" fixed  in patches of "+subNfixed
                + " with "+patches.size()+" possible patches");
        
        int maxSub  = (int)((double)L/R);
        int xx=0;int yy=0;
        // Paste all over lattice
        // choose patch
        for(int mx = 0;mx < maxSub;mx++){for(int my = 0;my < maxSub;my++){
            currentPatch = (patchInd)%patches.size();
            if(((mx+maxSub*my)%y) ==0){patchInd++;}
            if(type.contains("sections") || type.contains("Sections")){
                if(type.contains("Mix") || type.contains("mix")){
                    if(mx > ((int)(maxSub/2)-1)){currentPatch = 0;
                    }else{
                        if(mx %2 ==0){
                        currentPatch=1;
                        }else{
                            currentPatch=2;
                        }
                    }
                }else{    
                    if(mx > (int)(maxSub/2)){currentPatch = 1;}else{currentPatch=0;}
                } 
            }
            // pasting of patches
            int patchNfix = 0;
            for(int k = 0;k < subNfixed;k++){
                int t = patches.get(currentPatch)[k];
                int tx = t%R;
                int ty = ((int)((double)t/R))%R;        
                xx = mx*R+tx;
                yy = my*R+ty;
                if(type.contains("balanced") || type.contains("Balanced")){
                    if(patchNfix%2==0){fixed[nfixed][2] = s;}else{fixed[nfixed][2] = -1*s;sum += s;}
                }else{fixed[nfixed][2] = s;sum += s;}  
                fixed[nfixed][0] = xx;
                fixed[nfixed][1] = yy;
                nfixed++;
                patchNfix++;
            }
        }}
    }
    
    
    if(type.equalsIgnoreCase("square")||type.equalsIgnoreCase("rectangle")){
        for(int i = 0;i<L1;i++){
            for(int j = 0;j<L2;j++){
                u= (x-(L1/2)+i+L)%L;
                v= (y-(L2/2)+j+L)%L;

                nfixed++;
                fixed[nfixed][0]= u;
                fixed[nfixed][1]= v;
                fixed[nfixed][2]= s;
                sum += s;
            }
        }
    }
	
    if(type.equalsIgnoreCase("box")){
        for(int i = 0;i<L1;i++){
            for(int j = 0;j<L2;j++){
                if(i==0 || i==(L1-1)||j==(L2-1)||j==0){
                u= (x-(L1/2)+i+L)%L;
                v= (y-(L2/2)+j+L)%L;

                nfixed++;
                fixed[nfixed][0]= u;
                fixed[nfixed][1]= v;
                fixed[nfixed][2]= s;}
            }
        }
    }
	
    if(type.equalsIgnoreCase("cross")){
        for(int i = 0;i<L1;i++){
            for(int j = 0;j<L2;j++){
                if(i==(L1/2)||j==(L2/2)){
                u= (x-(L1/2)+i+L)%L;
                v= (y-(L2/2)+j+L)%L;

                nfixed++;
                fixed[nfixed][0]= u;
                fixed[nfixed][1]= v;
                fixed[nfixed][2]= s;sum += s;}
            }
        }
    }
	
    if(type.equalsIgnoreCase("crack")){

        for(int j = 0;j<L2;j++){

            v= (y-(L2/2)+j+L)%L;

            nfixed++;
            fixed[j][0]= v;
            fixed[j][1]= x;
            fixed[j][2]= s;sum += s;}
    }
	
    if(type.equalsIgnoreCase("circle")){
        double r;double r2=(L1/2.0)*(L1/2.0);
        for(int i = 0;i<L1;i++){
            for(int j = 0;j<L1;j++){
                u= x-(L1/2)+i;
                v= y-(L1/2)+j;
                r = (u-x)*(u-x)+(v-y)*(v-y);

                if(r<r2){
                        nfixed++;
                fixed[nfixed][0]= u;
                fixed[nfixed][1]= v;
                fixed[nfixed][2]= s;sum += s;}
            }
        }
    }	

	
    DirAndFileStructure dir = new DirAndFileStructure();
    
    String fName = dir.getFixedDirectory()+"fixed"+type+".txt";
    
    // Regularly just name after type but for postfix case use the postfix
    if(!postFix.equalsIgnoreCase("")){
        fName = dir.getFixedDirectory()+"fixed"+postFix+".txt";
    }

    // Make file just in case
    try {
        PrintStream out = new PrintStream(new FileOutputStream(
            fName));
        out.println();
        out.close();
    } catch (FileNotFoundException e) {
        e.printStackTrace();
    }

    try {
        PrintStream out = new PrintStream(new FileOutputStream(
            fName,false));
        int sumCheck = 0;
        for (int i = 0; i < fixed.length; i++){
                if(fixed[i][0]==0 && fixed[i][1]==0 && fixed[i][2]==0){}else{
                            if(type.equalsIgnoreCase("balancedrandom") || 
                                    type.equalsIgnoreCase("randombalanced")){sumCheck = sumCheck+fixed[i][2];}
                out.println(fixed[i][0]+ "     "+fixed[i][1]+ "     "+fixed[i][2]+ "     ");}
        }
        if(type.equalsIgnoreCase("balancedrandom") || type.equalsIgnoreCase("randombalanced")){System.out.println("Sum of Fixed: "+sumCheck);}
        System.out.println("MakeFixedConfig | Done with creating fixed");
        out.println();

        out.close();
    } catch (FileNotFoundException e) {
        e.printStackTrace();
    }

    System.out.println("Done with fixing "+nfixed+" of config  "+type+"   sum: "+sum+"  in file : "+fName+"  !");
    return fixed;
}


public static void makeFixedRandomNetwork(String type,String postFix,int nfix,int s, int N, int seed){
    int[][] fixed;
    Random Ran = new Random();
    Ran.setSeed(seed);
    int nfixed = 0;
    int u=0;
    fixed = new int[2*nfix+1][2];
    
    // need to add something to make sure not in same
    HashMap<Integer,Integer> currFixed = new HashMap<Integer,Integer>();
    for(int i = 0;i<(N);i++){currFixed.put(i, 0);}
    
    
    if(type.equalsIgnoreCase("highestdegree")){
        for(int j = 0;j<nfix;j++){
            fixed[j][0]= j;
            nfixed++;
            if(type.contains("balanced") && ((j%2)==1)){
                fixed[j][1]= (-1)*s;
            }else{                    
                fixed[j][1]= s;
            }
        }        
    }else if(type.equalsIgnoreCase("lowestdegree")){
        for(int j = 0;j<nfix;j++){
            fixed[j][0]= N-1-j;
            nfixed++;
            if(type.contains("balanced") && ((j%2)==1)){
                fixed[j][1]= (-1)*s;
            }else{                    
                fixed[j][1]= s;
            }
        }
    }else{
        boolean goodSpin = false;
        for(int j = 0;j<nfix;j++){
            goodSpin = false;
            while(!goodSpin){
                u = (int) (Ran.nextDouble()*N) ;
                // if spin hasnt been fixed then continue
                if(currFixed.get(u)== 0){goodSpin=true;}    
            }

            nfixed++;
            fixed[j][0]= u;
            currFixed.put((u), 1);
            if(type.contains("balanced") && ((j%2)==1)){
                fixed[j][1]= (-1)*s;
            }else{                    
                fixed[j][1]= s;
            }
        }
    }
    
    DirAndFileStructure dir = new DirAndFileStructure();
    
    String fName = dir.getFixedDirectory()+"fixed"+type+".txt";
    
    // Regularly just name after type but for postfix case use the postfix
    if(!postFix.equalsIgnoreCase("")){
        fName = dir.getFixedDirectory()+"fixed"+postFix+".txt";
    }

    // Make file just in case
    try {
        PrintStream out = new PrintStream(new FileOutputStream(
            fName));
        out.println();
        out.close();
    } catch (FileNotFoundException e) {
        e.printStackTrace();
    }

    try {
        PrintStream out = new PrintStream(new FileOutputStream(
            fName,false));
        int sumCheck = 0;
        for (int i = 0; i < fixed.length; i++){
                if(fixed[i][0]==0 && fixed[i][1]==0 && fixed[i][2]==0){}else{
                            if(type.equalsIgnoreCase("balancedrandom")){sumCheck = sumCheck+fixed[i][2];}
                out.println(fixed[i][0]+ "     "+fixed[i][1]);}
        }
        if(type.equalsIgnoreCase("balancedrandom")){System.out.println("Sum of Fixed: "+sumCheck);}

        out.println();

        out.close();
    } catch (FileNotFoundException e) {
        e.printStackTrace();
    }

    System.out.println("Done with fixing "+nfixed+" of config  "+type+"  in file : "+fName+"  !");
}
/**
*       makeMovieSet makes a set of configurations that are compatible with 
*   the example movie simulation. It also makes one configuration for every type
*   of fixed configuration supported
*/
public static void makeMovieSet(){
    MakeFixedConfig.makeFixedConfig2D("square", -1, 140, 140, 256,8,8);
    MakeFixedConfig.makeFixedConfig2D("circle", -1, 140, 140, 256,9,64);
    MakeFixedConfig.makeFixedConfig2D("crack", -1, 140, 140, 256,9,64);
    MakeFixedConfig.makeFixedConfig2D("cross", -1, 140, 140, 256,32,32);	
    MakeFixedConfig.makeFixedConfig2D("random", -1, 140, 140, 256,32,64);
}

/**
*       movieSetExists checks if a compatible set exists to run the movie examples.
*   It also checks for one configuration existing for each of the nonrandom 
*   fixed spin configurations/
* 
* @return  true if complete compatible set already exists
*/
public static boolean movieSetExists(){
    boolean setHere=true;
    DirAndFileStructure dir = new DirAndFileStructure();
    String fName = dir.getFixedDirectory()+"fixed";
    File f;
    f = new File((fName+"cross.txt"));
    if(!f.exists()){setHere=false;}
    f = new File((fName+"circle.txt"));
    if(!f.exists()){setHere=false;}
    f = new File((fName+"square.txt"));
    if(!f.exists()){setHere=false;}
    f = new File((fName+"crack.txt"));
    if(!f.exists()){setHere=false;}

    return setHere;
}

/**
*       makeConfig2D queries the user through the command line for all 
*   the parameters necessary to make a fixed configuration of the input type.
* 
* @param type - type of configuration to make
*/
public static void makeConfig2D(String type){
    int s;
    int x;
    int y;
    int l1;
    int l2;
    int L;
    
    Scanner scan = new Scanner(System.in);
    
    System.out.println("Spin of the fixed spins?");
    s = scan.nextInt();
    System.out.println("Center x coordinate of configuration?");
    x = scan.nextInt();
    System.out.println("Center y coordinate of configuration?");
    y = scan.nextInt();
    System.out.println("Length of lattice?");
    L = scan.nextInt();
    System.out.println("First parameter of configuration? ( D for circle or L of square)");
    l1 = scan.nextInt();
    System.out.println("Second parameter of configuration? (L for crack)");
    l2 = scan.nextInt();
    
    MakeFixedConfig.makeFixedConfig2D(type, s,x,y,L,l1,l2);	
}

public String typeIntercationQuilt(){
    return "InteractionQuilt";
}
public String typeIntercationQuiltSections(){
    return "InteractionQuiltSections";
}
public String typeRandom(){
    return "InteractionQuilt";
}
public String typeRandomBalanced(){
    return "randombalanced";
}
public String typeRandomBalancedConserved(){
    return "randombalancedconserved";
}

public String typeSquare(){
    return "square";
}
public String typeBox(){
    return "box";
}
public String typeCrack(){
    return "crack";
}
public String typeRectangle(){
    return "rectangle";
}
public String typeCircle(){
    return "circle";
}

// test the class
public static void main(String[] args) {
    String response="";
    Scanner scan = new Scanner(System.in);
    //MakeFixedConfig.makeMovieSet();
    boolean doneRes = true;

    
    while(doneRes){
	System.out.println("Type of configuration? (circle,crack,square,random,box,cross)");
	response =scan.next();
	if(response.equalsIgnoreCase("circle") || response.equalsIgnoreCase("cross")
			|| response.equalsIgnoreCase("crack")||response.equalsIgnoreCase("box")
			||response.equalsIgnoreCase("rectangle")||response.equalsIgnoreCase("square")
                        || response.equalsIgnoreCase("random")|| response.contains("balancedrandom")){
		doneRes=false;
	}else{}
    }
	
    // type of fixed, spin initial, x, y, L of sytem, L1 and L2 of geometry
    MakeFixedConfig.makeConfig2D(response);
    //cricle 69 box 68 square 64 crack 68 cross 67
    // L1 is radius of circle , L2 is length of crack
    // L2 foro random is amount fixed
    }

}
