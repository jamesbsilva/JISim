package Backbone.Util;

/**
*    @(#) ImgOpen
*/  

import AnalysisAndVideoBackend.VideoMakerCV;
import Backbone.Visualization.SquareLattice2D;
import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import javax.imageio.ImageIO;

/**
*      ImgOpen is a class to deal with opening image files and turning them into
*   data and does some data analysis to feature these functions.
* 
* <br>
* 
* @author      James Silva <jbsilva @ bu.edu>                 
* @since       2012-03
*/
public class ImgOpen {
    private int L; private int range; private int minScale = 5; private int minImageSize = 600;
    private int instId = 0; private boolean clusterMode=false;
    private DirAndFileStructure dir; private ParameterBank param;
    private int imageTextOffset = 90; private boolean imageTextOffsetSet = false;
    private boolean normalizeAvgImage = true;
    private SquareLattice2D sq = new SquareLattice2D(1,1,0.0,true,false);
	
    /**
    * 
    * @param id - Instance Id
    */
    public ImgOpen(int id){
            this(id,"");
    }

    /**
    * 
    * @param id - Instance Id
    * @param paramPost - postscript for parameter bank file 
    */
    public ImgOpen(int id, String paramPost){
        dir = new DirAndFileStructure();instId=id;param = new ParameterBank(paramPost);
        L = param.L; range = param.R;
    }

    /**
    *      setClusterMode changes image interpretation into cluster mapping.
    * 
    * @param st - cluster mode on or off (false).
    */
    public void setClusterMode(boolean st){ clusterMode = st; }
    
    /**
     *      setSystemL sets system size L.
     * @param l - system size 
     */
    public void setSystemL(int l){ L = l; }
    
    /**
     *      setSystemR sets system range R.
     * @param r - system size 
     */
    public void setSystemR(int r){ range = r; }
    
    /**
     *      setImgTextOffset sets image text offset variable.
     * @param off - offset value 
     */
    public void setImgTextOffset(int off){ imageTextOffset=off; }
    
    /**
     *      getImgTextOffset gets image text offset variable. 
     */
    public int getImgTextOffset(){ return imageTextOffset; }
    
    /**
    *    compileDirIntoAvgArray takes one of the Saved Config directories
    *   and averages out all the images in this directory into an float array.
    * 
    * @param nfixed - number fixed for the Saved Configs directory. 
    *       directory in form Saved Configs/Fixed nfixed/
    * @return average spin array
    */
    public float[][] compileDirIntoAvgArray(int nfixed){
        File[] images = getSavedConfigFileList(nfixed,("L-"+L));
        return compileDirIntoAvgArray(images);
    }
    
    
    public float[][] compileDirIntoAvgArray(String directory){
        File f = new File(directory);
        File[] images = f.listFiles();
        return compileDirIntoAvgArray(images);
    }
    
    public float[][] compileDirIntoAvgArray(File[] images){
        float[][] lat = new float[L][L]; int[][] sumLat = new int[L][L];int NumImg = images.length;
        imageTextOffset = calculateImageTextOffset(getBufferedImage(images[0]));

        for( int k = 0; k < NumImg; k++ ){
            if(clusterMode){
                sumLat = updateClusterArrayUsingImg(getBufferedImage(images[k]), sumLat);
            }else{
                sumLat = updateArrayUsingImg(getBufferedImage(images[k]), sumLat);
            }
        }

        double frac = 0; double sumMax = -1; double sumMin = 1;
        for( int i = 0; i < L; i++ ){ for( int j = 0; j < L; j++ ){
            frac = ((double)sumLat[i][j])/((double)NumImg);
            if(frac > sumMax){sumMax=frac;}
            if(frac < sumMin){sumMin=frac;}
        }}
        System.out.println("Max Sum : "+sumMin);
        //Normalize to fit in range 0-1
        for( int i = 0; i < L; i++ ){ for( int j = 0; j < L; j++ ){
            // Fit into -1/2 to 1/2
            frac = ((double)sumLat[i][j])/((double)NumImg);
            if(normalizeAvgImage){
                lat[i][j] = (float)((frac-sumMin)/(2.0f*(sumMax-sumMin))+0.5f);
            }else{
                lat[i][j] = (float)(frac)/(2.0f);
                lat[i][j] = lat[i][j]+0.5f;
            }        
        }}

        return lat;
    }
	 
    /**
    *       makeDirIntoAvgImg takes one of the Saved Config directories
    *   and averages out all the images in this directory into an image.
    * 
    * @param nfixed - number fixed for the Saved Configs directory. 
    *       directory in form Saved Configs/Fixed nfixed/
    * @param fname - output image filename
    */
    public void makeDirIntoAvgImg(int nfixed,String fname){
        float[][] lat;
        
        // compile directory
        lat = compileDirIntoAvgArray(nfixed);
        
        // make and save image
        DataSaver d = new DataSaver(instId,param.postFix);
        int n = getSavedConfigNumImages(nfixed,("L-"+L));
        String outString = "nfixed: "+ nfixed+"    nruns: "+n;
        d.saveFloatArrayAsImgInImages(lat, fname,outString);
    }

    /**
    *       makeDirIntoAvgImg takes a directory and averages out all the 
    *   images in this directory into an image.
    * 
    * @param dir - directory to compile
    * @param fname - output image filename
    */
    public void makeDirIntoAvgImg(String dir,String fname){
        float[][] lat;
        
        // compile directory
        lat = compileDirIntoAvgArray(dir);
        
        // make and save image
        DataSaver d = new DataSaver(instId,param.postFix);
        d.saveFloatArrayAsImgInImages(lat, fname);
    }
    
    /**
    *       makeMovieOfDirAvg takes one of the Saved Config directories
    *   and averages out all the images in this directory and makes a movie of 
    *   process.
    * 
    * @param nfixed - number fixed for the Saved Configs directory. 
    *       directory in form Saved Configs/Fixed nfixed/
    * @param numOfCaps - number of images to process
    */
    public void makeMovieOfDirAvg(int nfixed, int numOfCaps){
        clusterMode = false;
        makeMovieOfDirAvg(nfixed,numOfCaps);
    }
    
    /**
    *       makeMovieOfDirAvg takes a directory (defaults to interpret using cluster mode)
    *   and averages out all the images in this directory and makes a movie of 
    *   process.
    * 
    * @param directory - directory to make a movie of averaging process
    * @param numOfCaps - number of images to process
    * @param cluster - movie of cluster mode
    */
    public void makeMovieOfDirAvg(String directory, int numOfCaps,boolean cluster){
        clusterMode = cluster;
        makeMovieOfDirAvg(0,numOfCaps,directory);
    }
    
    
    /**
    *       makeMovieOfDirAvg takes a directory (defaults to interpret using cluster mode)
    *   and averages out all the images in this directory and makes a movie of 
    *   process.
    * 
    * @param directory - directory to make a movie of averaging process
    * @param numOfCaps - number of images to process
    */
    public void makeMovieOfDirAvg(String directory, int numOfCaps){
        clusterMode = true;
        makeMovieOfDirAvg(0,numOfCaps,directory);
    }
    
    /**
    *       makeMovieOfDirAvg takes one of the Saved Config directories
    *   and averages out all the images in this directory and makes a movie of 
    *   process.
    * 
    * @param nfixed - number fixed for the Saved Configs directory. 
    *       directory in form Saved Configs/Fixed nfixed/
    * @param numOfCaps - number of images to process
    * @param directory - directory to compile if in cluster mode
    */
    private void makeMovieOfDirAvg(int nfixed, int numOfCaps, String directory){
        // Make Video
        VideoMakerCV vid=null; BufferedImage img; DataSaver dsave = new DataSaver(); File[] images;
        if(clusterMode){
            File f = new File(directory); images = f.listFiles();
        }else{
            images = getSavedConfigFileList(nfixed,("L-"+L));
        }
        int NumImg = images.length;

        float[][] lat = new float[L][L];
        int framesPerCapture = 1;
        if( numOfCaps > NumImg ){ numOfCaps = NumImg; System.out.println("ImgOpen Error: Not enough images for captures will only do max possible");}

        int VideoLength = numOfCaps * framesPerCapture;
        
        String type = "averaging"+nfixed;
        vid = new VideoMakerCV(0,VideoLength,type);
        vid.setFramesPerCapture(framesPerCapture);

        int[][] sumLat = new int[L][L]; int k=0;
        imageTextOffset = calculateImageTextOffset(getBufferedImage(images[0]));        
        
        while( k < numOfCaps ){
            if(clusterMode){
                sumLat = updateClusterArrayUsingImg(getBufferedImage(images[k]), sumLat);
            }else{
                sumLat = updateArrayUsingImg(getBufferedImage(images[k]), sumLat);
            }

            if( k % 50 == 0 ){ System.out.println("ImgOpen Images Processed: "+k); }

            int sum;
            //Normalize to fit in range 0-1
            for( int i = 0; i < L; i++ ){ for( int j = 0; j < L; j++ ){
                // Fit into -1/2 to 1/2
                sum = sumLat[i][j];
                lat[i][j] = (sum)/(2.0f*(k+1));
                lat[i][j] = lat[i][j]+0.5f;
            }}	 
            vid.addLatticeFrame(dsave.getImageFromArray(lat));
            k++;
        }
        
        // write the video
        vid.writeVideo();
    }
    
    public void processImagesForFixedDensity(String directory){
        File f = new File(directory);
        File[] images = f.listFiles(); int NumImg = images.length;
        imageTextOffset = calculateImageTextOffset(getBufferedImage(images[0]));
        imageTextOffsetSet=true;
        System.out.println("Image Offset For Text: "+imageTextOffset);
        for( int k = 0; k < NumImg; k++ ){
            createFixedDensityImage(getBufferedImage(images[k]), images[k].getAbsolutePath(),true,false);
        }
    }
    
    public void createFixedDensityImage(String fname){
        File f = new File(fname);
        createFixedDensityImage(getBufferedImage(f),fname,true,true);
    }
        
    public void createFixedDensityImage(BufferedImage img, String name,boolean normalizeMax,boolean saveText){
        int[][] lat = new int[L][L]; int[][] fixed = new int[L][L]; float[][] fixDensity = new float[L][L];
        if(!imageTextOffsetSet){
            imageTextOffset = calculateImageTextOffset(img);
            System.out.println("ImageOffset: "+imageTextOffset);
        }
        lat = convertImgToSpinIntArray(img); fixed = convertImgToFixIntArray(img);
        int Q = ( range == 0 ) ? 4 : (2*range+1)*(2*range+1)-1;
        
        int x; int y; int sum; int sumMax = 0; int sumMin = Q; int sumIJ=0;
        int sumF = 0; int sumF2 = 0;
        for(int i = 0 ; i < L;i++){for(int j = 0 ;j < L;j++){
            sum = 0;
            if(fixed[i][j]==1){sumIJ += lat[i][j];}
            if(range==0){
                x= ((i-1+L)%L);y=j;if(fixed[x][y]==1){sum +=lat[x][y];}
                x= ((i+1+L)%L);y=j;if(fixed[x][y]==1){sum +=lat[x][y];}
                x= i;y=((j-1+L)%L);if(fixed[x][y]==1){sum +=lat[x][y];}
                x= i;y=((j+1+L)%L);if(fixed[x][y]==1){sum +=lat[x][y];}
            }else{
                for(int u = 0 ; u < (2*range+1);u++){for(int v = 0 ;v < (2*range+1);v++){
                    x = ((i-range+u+L)%L); y = ((j-range+v+L)%L);
                    if( fixed[x][y] == 1 ){
                        if( x != i || y != j ){ sum += lat[x][y]; }
                    }
                }}
            }
            if(sumMax < sum){sumMax = sum;} if(sumMin > sum){sumMin = sum;}
            fixDensity[i][j] = sum;
            sumF += Math.abs(sum); sumF2 += sum*sum;
        }}
        double mean = ((double)sumF/(double)(L*L));
        double fano = ((double)sumF2/(double)(L*L))- mean*mean;
        fano = fano/mean;
        System.out.println("Sum of Fixed Spins: "+sumIJ+"      max : "+sumMax+"      min : "+sumMin+"     fano : "+fano);
        
        String[] set = name.split("/");
        String root = set[set.length-1];
        String rootSubset="";
        
        if(!(name.substring(0,1)).contains("/")){
            for(int u =0;u<(set.length-1);u++){ rootSubset += set[u]+"/"; }
        }
        
        String dir = name.substring(0,name.length()-root.length());         
        name = name.substring(name.length()-root.length(),name.length());
        set = name.split(".png");
        root = set[set.length-1];
        String edLat = root; DataSaver d = new DataSaver(instId,param.postFix); d.setL(L);
        // normalize
        for( int i = 0 ; i < L; i++ ){ for( int j = 0 ; j < L; j++ ){
            if(saveText){
                double[] density = new double[3];
                density[0] = i; density[1] = j; density[2] = fixDensity[i][j];
                d.saveDoubleArrayData(density, rootSubset+edLat+"-FieldDensity.txt");
            }
            if(normalizeMax){
                fixDensity[i][j] = ((fixDensity[i][j]-sumMin)/(2.0f*(sumMax-sumMin))+0.5f);
            }else{
                fixDensity[i][j] = (fixDensity[i][j]/(4.0f*Q))+0.5f;
            }
        }}
        // make and save image
        d.saveFloatArrayAsImgInImages(fixDensity,(rootSubset+edLat+"-FixDensity.png"),"");
    }
    public void createSpinDensityImage(BufferedImage img, String name,boolean normalizeMax,boolean saveText){
        int[][] lat = new int[L][L]; int[][] fixed = new int[L][L]; float[][] fixDensity = new float[L][L];
        if( !imageTextOffsetSet ){
            imageTextOffset = calculateImageTextOffset(img);
            System.out.println("ImageOffset: "+imageTextOffset);
        }
        lat = convertImgToSpinIntArray(img); fixed = convertImgToFixIntArray(img);
        int Q = ( range == 0 ) ? 4 : (2*range+1)*(2*range+1)-1;
        
        int x; int y; double sum; double sumAll; double sumMax = 0; double sumMin = Q; int sumIJ=0;
        for( int i = 0 ; i < L; i++ ){ for( int j = 0 ;j < L; j++ ){
            sum = 0; sumAll = 0;
            if( fixed[i][j] == 1 ){ sumIJ += lat[i][j]; }
            if( range == 0 ){
                x= ((i-1+L)%L);y=j;if(fixed[x][y]==1){sum +=lat[x][y];}
                x= ((i+1+L)%L);y=j;if(fixed[x][y]==1){sum +=lat[x][y];}
                x= i;y=((j-1+L)%L);if(fixed[x][y]==1){sum +=lat[x][y];}
                x= i;y=((j+1+L)%L);if(fixed[x][y]==1){sum +=lat[x][y];}
            }else{
                for( int u = 0 ; u < (2*range+1); u++ ){ for(int v = 0 ; v < (2*range+1); v++ ){
                    x = ((i-range+u+L)%L); y = ((j-range+v+L)%L);
                    if( fixed[x][y] != 1 ){
                        if( x != i || y != j ){
                            sum += lat[x][y];
                            sumAll++;
                        }
                    }
                }}
            }
            if(sumMax < sum){sumMax = sum;} if(sumMin > sum){sumMin = sum;}
            fixDensity[i][j] = (float)(sum/sumAll);
        }}
        System.out.println("Sum of Fixed Spins: "+sumIJ);
        
        String[] set = name.split("/");
        String root = set[set.length-1];
        String rootSubset="";
        
        if(!(name.substring(0,1)).contains("/")){
            for( int u = 0; u < (set.length-1); u++ ){
                rootSubset += set[u]+"/";
            }
        }
        
        String dir = name.substring(0,name.length()-root.length());         
        name = name.substring(name.length()-root.length(),name.length());
        set = name.split(".png");
        root = set[set.length-1];
        String edLat = root;
        DataSaver d = new DataSaver(instId,param.postFix); d.setL(L);
        // normalize
        for(int i = 0 ; i < L;i++){for(int j = 0 ;j < L;j++){
            if(saveText){
                double[] density = new double[3];
                density[0] = i; density[1] = j;
                density[2] = fixDensity[i][j];
                d.saveDoubleArrayData(density, rootSubset+edLat+"-SpinDensity.txt");
            }
            if(normalizeMax){
                fixDensity[i][j] = (float)((fixDensity[i][j]-sumMin)/(2.0f*(sumMax-sumMin))+0.5f);
            }else{
                fixDensity[i][j] = (fixDensity[i][j]/(4.0f*Q))+0.5f;
            }
        }}
        
        // make and save image
        d.saveFloatArrayAsImgInImages(fixDensity,(rootSubset+edLat+"-SpinDensity.png"),"");
    }
    
    /**
    *       getBufferedImage opens the image file and returns the buffered image
    *   object.
    * 
    * @param pic - image file to process
    * @return bufferedimage object
    */
    public BufferedImage getBufferedImage(File pic){
        BufferedImage image = null;
        try {
            image = ImageIO.read(pic);
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
            System.err.print(pic.getName()+"  : file not found / retrieved.");
        }
        return image;
    }

    /**
    *       getSavedConfigFileList gets the list of filenames of directory
    *   that are located in directory of saved configs with given fixed number.
    * 
    * @param nfixed - number fixed for the Saved Configs directory. 
    *       directory in form Saved Configs/Fixed nfixed/ 
    * @return Filename of images in directory
    */
    public File[] getSavedConfigFileList(int nfixed, String configSavePost){
        // Find amount of files
        String cdirect = dir.getSavedConfigDirectory(instId)+"Fixed_"+nfixed+"_"+configSavePost+"/";
        System.out.println("ImgOpen Taking files from directory : "+cdirect);
        File path = new File(cdirect);
        File[] files = path.listFiles();
        return files;
    }

    /**
    *       getSavedConfigNumImages gets the number of filenames of directory
    *   that are located in directory of saved configs with given fixed number.
    * 
    * @param nfixed - number fixed for the Saved Configs directory. 
    *       directory in form Saved Configs/Fixed nfixed/
    * @return 
    */
    private int getSavedConfigNumImages(int nfixed, String configSavePost){
        // Find amount of files
        String cdirect = dir.getSavedConfigDirectory(instId)+"Fixed_"+nfixed+"_"+configSavePost+"/";
        File path = new File(cdirect);
        File[] files = path.listFiles();
        return files.length;
    }
	
    /**
    *       getSpinArray takes a full pathfile name and turns the image file
    *   into an array of the spin values.
    * 
    * @param fname - full filename of image to turn into spin array
    * @return spin value array
    */
    public int[][] getSpinArray(String fname){
        File f = new File(fname); BufferedImage img = getBufferedImage(f);
        imageTextOffset = calculateImageTextOffset(img);
        return convertImgToSpinIntArray(img);
    }
	
    /**
    *       getFixedArray takes a full pathfile name and turns the image file
    *   into an array of the fixed spin boolean 1 or 0 values.
    * 
    * @param fname - full filename of image to turn into spin array
    * @return fixed 1 or 0 array
    */  
    public int[][] getFixedArray(String fname){
        File f = new File(fname); BufferedImage img = getBufferedImage(f);
        imageTextOffset = calculateImageTextOffset(img);
        return convertImgToFixIntArray(img);
    }
		
    /**
    *       convertImgToSpinIntArray converts a given image into an array
    *   of the spin values.
    * 
    * @param img - image to process
    * @return spin value array
    */
    private int[][] convertImgToSpinIntArray(BufferedImage img){
        // Scale default by 5
        int scale = minScale; int sizeImage = L*scale;
        while(sizeImage < minImageSize){scale = scale*2;sizeImage = scale*L;}
        int clr; int[][] lattice = new int[L][L];
        //System.out.println("L: "+L+"  scale : "+scale+"    H : "+img.getHeight()+"    W : "+img.getWidth()+"    txtoff: "+imageTextOffset);
        for( int i = 0; i < L; i++ ){ for( int j = 0; j < L; j++ ){
            // Getting pixel color by position x=100 and y=40 
            clr =  img.getRGB((int)Math.floor((i+0.5)*scale),(int)Math.floor((j+0.5)*scale)+ imageTextOffset);
            lattice[i][j] = getSpinFromRGB(clr);
            if( lattice[i][j] == -2222 ){System.out.println("ImgOpen Error : Spin Pixel color out of range. | "+clr);}
        }}        
        return lattice;		
    }
	
    /**
    *       convertImgToFixIntArray converts a given image into an array
    *   of the fixed spin 1 or 0 values.
    * 
    * @param img - image to process
    * @return fixed 1 or 0 array
    */
    private int[][] convertImgToFixIntArray(BufferedImage img){
        // Scale default by 5
        int scale = minScale; int sizeImage = L*scale; int clr; int[][] lattice = new int[L][L];
        while(sizeImage < minImageSize){scale = scale*2;sizeImage = scale*L;}
        
        for( int i = 0; i < L; i++ ){ for( int j = 0; j < L; j++ ){
            // Getting pixel color by position x=100 and y=40 
            clr =  img.getRGB((int)((i+0.5)*scale),(int)((j+0.5)*scale)+ imageTextOffset);
            lattice[i][j] = getFixedFromRGB(clr);
            if(lattice[i][j]==-2222){System.out.println("ImgOpen Error : Pixel fixed color out of range. | "+clr);}
        }}
        return lattice;
    }

    /**
    *       updateArrayUsingImg updates an integer array with the values
    *   from a given image.
    * 
    * @param img - image to process
    * @return updated additive spin array
    */
    private int[][] updateArrayUsingImg(BufferedImage img, int[][] lat){
        // Scale default by 5
        int clr; int scale = minScale; int sizeImage = L*scale;
        while(sizeImage < minImageSize){scale = scale*2;sizeImage = scale*L;}
        boolean useFix = false;
        
        for(int i=0;i<L;i++){ for(int j=0;j<L;j++){
            // Getting pixel color by position x=100 and y=40 
            clr =  img.getRGB((int)Math.floor((i+0.5)*scale),(int)Math.floor((j+0.5)*scale)+ imageTextOffset);
            if(useFix){
                lat[i][j] = lat[i][j]+getSpinFromRGB(clr);
            }else{
                lat[i][j] = (getFixedFromRGB(clr) == 1) ? lat[i][j] : lat[i][j]+getSpinFromRGB(clr);
            }
                      
        }}		
        return lat;
    }

    
    /**
    *       updateArrayUsingImg updates an integer array with the values
    *   from a given image.
    * 
    * @param img - image to process
    * @return updated additive spin array
    */
    private int[][] updateClusterArrayUsingImg(BufferedImage img, int[][] lat){
        // Scale default by 5
        int scale = minScale; int clr; int sizeImage = L*scale;
        while(sizeImage < minImageSize){scale = scale*2;sizeImage = scale*L;}
        
        for( int i = 0; i < L; i++ ){ for( int j = 0; j < L; j++ ){
            // Getting pixel color by position x=100 and y=40 
            clr =  img.getRGB(i*scale,j*scale+ imageTextOffset);
            lat[i][j] = lat[i][j]+getClusterFromRGB(clr);       
        }}		
        return lat;
    }

    /**
    *       getClusterFromRGB converts and RGB color code to a spin value.
    * 
    * @param col - RGB code of color
    * @return cluster value corresponding to RGB code
    */
    private int getClusterFromRGB(int col){
        Color jCol; int c1; int cluster = 0;
        jCol = Color.ORANGE;c1 = jCol.getRGB();
        if(col == c1){ cluster= 1; }
        return cluster;
    }
    
    /**
    *       getSpinFromRGB converts and RGB color code to a spin value.
    * 
    * @param col - RGB code of color
    * @return spin value corresponding to RGB code
    */
    private int getSpinFromRGBold(int col){
        Color jCol; int c1; int c2; int spin = -2222;

        jCol = Color.CYAN;c1 =jCol.getRGB();jCol = Color.RED;c2 = jCol.getRGB();
        if(col == c1 || col == c2){spin= -1;}
        jCol = Color.YELLOW;c1 =jCol.getRGB();jCol = Color.GRAY;c2 = jCol.getRGB();
        if(col == c1 || col == c2){spin= 1;}
        jCol = Color.WHITE;c1 =jCol.getRGB();jCol = Color.MAGENTA;c2 = jCol.getRGB();
        if(col == c1 || col == c2){spin= 0;}

        return spin;
    }

    /**
    *       getSpinFromRGB converts and RGB color code to a spin value.
    * 
    * @param col - RGB code of color
    * @return spin value corresponding to RGB code
    */
    private int getSpinFromRGB(int col){
        Color jCol; int c1; int c2; int spin = -2222;

        //jCol = Color.RED;c1 =jCol.getRGB();jCol = Color.BLACK;c2 = jCol.getRGB();
        jCol = sq.getColor(-1);c1 =jCol.getRGB();jCol = sq.getColorFixed(-1);c2 = jCol.getRGB();
        if(col == c1 || col == c2){spin= -1;}
        //jCol = Color.BLUE;c1 =jCol.getRGB();jCol = Color.WHITE;c2 = jCol.getRGB();
        jCol = sq.getColor(1);c1 =jCol.getRGB();jCol = sq.getColorFixed(1);c2 = jCol.getRGB();
        if(col == c1 || col == c2){spin= 1;}
        //jCol = Color.YELLOW;c1 =jCol.getRGB();jCol = Color.MAGENTA;c2 = jCol.getRGB();
        jCol = sq.getColor(0);c1 =jCol.getRGB();jCol = sq.getColorFixed(0);c2 = jCol.getRGB();
        if(col == c1 || col == c2){spin= 0;}

        return spin;
    }
    
    /**
    *       getFixedFromRGB converts and RGB color code to a fixed spin 1 or 0 value.
    * 
    * @param col - RGB code of color
    * @return fixed spin value corresponding to RGB code
    */
    private int getFixedFromRGB(int col){
        Color jCol; int c1 = 0; int c2 = 0; int c3 = 0; int c4 = 0;
        int fixed = -2222;

        //jCol = Color.WHITE;c1 =jCol.getRGB();jCol = Color.YELLOW;c2 = jCol.getRGB();
        //jCol = Color.CYAN;c3 =jCol.getRGB();jCol = Color.BLACK;c4 = jCol.getRGB();
        jCol = sq.getColor(0);c1 =jCol.getRGB();jCol = sq.getColor(-1);c2 = jCol.getRGB();
        jCol = sq.getColor(1);c3 =jCol.getRGB();jCol = Color.BLACK;c4 = jCol.getRGB();
        if(col == c1 || col == c2 || col == c3){fixed= 0;}
        //jCol = Color.RED;c1 =jCol.getRGB();jCol = Color.GRAY;c2 = jCol.getRGB();
        //jCol = Color.MAGENTA;c3 = jCol.getRGB();jCol = Color.BLUE;c4 = jCol.getRGB();
        jCol = sq.getColorFixed(1);c1 =jCol.getRGB();jCol = sq.getColorFixed(-1);c2 = jCol.getRGB();
        jCol = sq.getColorFixed(0);c3 =jCol.getRGB();jCol = Color.MAGENTA;c4 = jCol.getRGB();
        if(col == c1 || col == c2 || col == c3){fixed= 1;}

        return fixed;
    }
    

    /**
    *       calculateImageTextOffset determines the length of any offset in given
    *   image which may have been added for text.
    * 
    * @param img - image to process
    * @return integer value of offset
    */
    public int calculateImageTextOffset(BufferedImage img){
        int off=0; int clr; int val; int k=0; 
        if( img.getHeight() == img.getWidth() ){
            return 0;
        }else{
            return img.getHeight()-img.getWidth();
        }
    }
	
    // test the class
    public static void main(String[] args) {
        // TODO Auto-generated method stub
        ImgOpen dat = new ImgOpen(0);
        /*int[][] lat = dat.getSpinArray("/home/j2/test.png");
        int[][] fix = dat.getFixedArray("/home/j2/test.png");
        dataSaver d = new dataSaver(0);
        d.saveIntArrayAsImg(lat, fix, "out.png");*/
        //int nfixed = 10;
        dat.setSystemR(8); dat.setSystemL(128); int nfixed = 0;
        //dat.processImagesForFixedDensity("/home/j2/JISim/Images/NucleationT/iqsm1b2");
        //dat.makeDirIntoAvgImg("/home/j2/Research/Nov29Nucleation/Untitled Folder 2/NucleationT/random1n/x-50-hfield-0.75/snapshot-seed-471450023/","ER"+nfixed+".png");
        //dat.makeDirIntoAvgImg("/home/j2/Research/Nov29Nucleation/Untitled Folder 2/NucleationT/random1n/x-0","ER"+0+".png");
        //dat.makeDirIntoAvgImg("/home/j2/Research/Nov29Nucleation/Untitled Folder 2/NucleationT/random1n/x-48-hfield-0.75/snapshot-seed-1185570770","ER"+48+".png");
        //dat.makeDirIntoAvgImg("/home/j2/Research/Nov29Nucleation/Untitled Folder 2/NucleationT/random1n/x-50-hfield-0.75/snapshot-seed-471450023","ER"+50+".png");
        //dat.makeDirIntoAvgImg("/home/j2/Research/Nov29Nucleation/Untitled Folder 2/NucleationT/random1n/x-52-hfield-0.75/snapshot-seed-972092929","ER"+52+".png");
        //dat.processImagesForFixedDensity("/home/j2/Research/Nov29Nucleation/Untitled Folder 2/NucleationT/random1n/x-52-hfield-0.75/snapshot-seed-972092929");
        // 0.42657439446
        //dat.createFixedDensityImage("/home/j2/Research/Nov29Nucleation/Untitled Folder 2/NucleationT/random1n/x-50-hfield-0.75/snapshot-seed-471450023/Config-t-11071-Run-211-field-0.75-temp-1.8-R-8-Seed-569073069.png");
        // 0.42311418685
        //dat.createFixedDensityImage("/home/j2/Research/Nov29Nucleation/Untitled Folder 2/NucleationT/random1n/x-48-hfield-0.75/snapshot-seed-1185570770/Config-t-339-Run-233-field-0.75-temp-1.8-R-8-Seed-2097349315.png");
        // 0.4092733564
        //dat.createFixedDensityImage("/home/j2/Research/Nov29Nucleation/Untitled Folder 2/NucleationT/random1n/x-52-hfield-0.75/snapshot-seed-972092929/Config-t-320-Run-97-field-0.75-temp-1.8-R-8-Seed-1754198076.png");
        //dat.createFixedDensityImage("/home/j2/Research/Nov29Nucleation/Untitled Folder 2/NucleationT/random1n/x-72-hfield-0.55/snapshot-seed-26713362/Config-t-10-Run-60-field-0.55-temp-1.8-R-8-Seed-1745294767.png");
        //dat.createFixedDensityImage("/home/j2/Research/Nov29Nucleation/Untitled Folder 2/NucleationT/random1n/x-74-hfield-0.55/snapshot-seed-106691818/Config-t-4135-Run-141-field-0.55-temp-1.8-R-8-Seed-1489776271.png");
        //dat.makeMovieOfDirAvg("/home/j2/Research/AnalysisScratchSpace/Oct30Img/StableRandomR5/Config",4,false);
        //dat.processImagesForFixedDensity("/home/j2/Research/Nov29/NucleationT/iqsm1p");
        //dat.createFixedDensityImage("/home/j2/JISim/Images/test9.png");
        //0  0.9261363636363636
        //50  0.1926605504587156  -27
        //48   0.3130434782608696 -28
        //52  0.4166666666666667  -32
        System.out.println("ImgOpen ______________________________");
        System.out.println("ImgOpen Done!!!");
        System.out.println("ImgOpen ______________________________");
    }
}
