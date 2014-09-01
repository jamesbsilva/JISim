package Backbone.System;


/**
 * 
 *    @(#) Neighbors Class 
 */  
 /** 
 *   Determines the neighbors for i,j.
 *  <br>
 * 
 *  @param i,j - the index for the lattice location
 *  @param N - System Size  (256 default)
 *  @param D - Dimension (2 default)
 *  @param Geo - Geometry,  2 is rectangular
 *  @param Neighbors
 * 
 * @author      James Silva <jbsilva @ bu.edu>                 
 * @since       2011-09    
 */

public class Neighbors {
    private int L = 256;
    private int D = 2;          // dimension
    private int Geo = 2;
    private LatticeMagInt sLat;
    private boolean useLongR = false;
    private int R = 1;
    private int neighSum;

    public Neighbors(LatticeMagInt inLat){
        sLat = inLat;
        L = sLat.getLength();
        D = sLat.getDimension();
        Geo = sLat.getGeo();

        R = sLat.getRange();

        if(R>1){useLongR=true;}else{useLongR=false;}
    }


    public int getLastNeighborSum(){return neighSum;}

    /**
    * getNeighbors                          (1)
    *
    * get the neighbors of a point on a lattice
    *
    */
    public int[] getNeighbors(int i, int j, int k){
        int[] neighbors = null;
        neighSum = 0;

        if (D==2 && Geo ==2){neighbors = getNextSq(i,j,0);}
        if (D==3 && Geo ==2){neighbors = getNextSq(i,j,k);}

        return neighbors;
    }

    private int[] getNextSq(int i, int j, int k){
        int [] neighbors;

        neighSum = 0;

        if(!useLongR){
        int sN = 2*D;
            neighbors = new int[sN];

        // Get neighbors deal with periodic conditions, Calculate Sum while getting values
        neighbors[0] = sLat.getValue(((i+1)%L),j,k);neighSum+=neighbors[0];
        neighbors[1] = sLat.getValue(i,((j+1)%L),k);neighSum+=neighbors[1];
        neighbors[2] = sLat.getValue(((i+L-1)%L),j,k);neighSum+=neighbors[2];
        neighbors[3] = sLat.getValue(i,((j+L-1)%L),k);neighSum+=neighbors[3];
        // 3d neighbors
        if(D==3){
                neighbors[4] = sLat.getValue(i,j,((k+1)%L));neighSum+=neighbors[4];
                neighbors[5] = sLat.getValue(i,j,((k+L-1)%L));neighSum+=neighbors[5];
                }

        // 3D Case of square range
        }else if (D==3){

        int sN = (2*R+1)*(2*R+1)*(2*R+1)-1;
        neighbors = new int[sN];int u;int v;int z;
        // use Sn as index
        sN = 0;

        for (int m=0;m<(2*R+1);m++){for (int n=0;n<(2*R+1);n++){for (int p=0;p<(2*R+1);p++){
            u= ((i-R+m+L)%L);
            v= ((j-R+n+L)%L);
            z= ((k-R+p+L)%L);

            if(u==i && v==j && z==k){}else{
            neighbors[sN] = sLat.getValue(u,v,z);
            //if(neighbors[sN]>1){System.out.println("spin >1");}
            neighSum+=neighbors[sN];
            sN++;}
        }}}

        // 2D Case of square long range	
        }else{

        int sN = (2*R+1)*(2*R+1)-1;
        neighbors = new int[sN];int u;int v;

        // use Sn as index
        sN = 0;

        for (int m=0;m<(2*R+1);m++){for (int n=0;n<(2*R+1);n++){
            u= ((i-R+m+L)%L);
            v= ((j-R+n+L)%L);

            if(u==i && v==j){}else{
            neighbors[sN] = sLat.getValue(u,v,k);neighSum+=neighbors[sN];
            sN++;}
            }}

        }

        return  neighbors;
    }




    public void printNeighbors(int i,int j,int k){
        int[] neighbors= getNeighbors(i,j,k);
        for(int m = 0;m < neighbors.length;m++){
            System.out.println(neighbors[m]);
        }
    }

    // test the class
    public static void main(String[] args) {
        // make L=4 square 2d lattice initialized to -1
        SimpleLattice lat = new SimpleLattice(4,3,2,2,-1);
        int i=2;
        int j=1;
        int k=0;
        int L = lat.getLength();

        System.out.println(lat.getValue(((i+1)%L),j,k));
        System.out.println(lat.getValue(i,((j+1)%L),k));
        System.out.println(lat.getValue(((i+L-1)%L),j,k));
        System.out.println(lat.getValue(i,((j+L-1)%L),k));
        Neighbors neigh = new Neighbors(lat);
        System.out.println("______________________________________________");
        neigh.printNeighbors(i, j, k);
        System.out.println("______________________________________________");
        System.out.println("Sum :"+neigh.getLastNeighborSum());
        System.out.println("______________________________________________");
        lat.setValue(2,2,0,3,5);
        neigh.printNeighbors(i, j, k);
        System.out.println("______________________________________________");
        System.out.println("Sum :"+neigh.getLastNeighborSum());
        System.out.println("Done!");
    }
}
