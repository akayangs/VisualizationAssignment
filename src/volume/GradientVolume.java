/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volume;
import volvis.RaycastRenderer;
/**
 *
 * @author michel
 */
public class GradientVolume {
    RaycastRenderer RR;

    public GradientVolume(Volume vol) {
        volume = vol;
        dimX = vol.getDimX();
        dimY = vol.getDimY();
        dimZ = vol.getDimZ();
        data = new VoxelGradient[dimX * dimY * dimZ];
        compute();
        maxmag = -1.0;
    }
    double getMagnitude(double[] vector){
        return Math.sqrt(vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2]);
    }
   
    public VoxelGradient getGradient(int x, int y, int z) {
        return data[x + dimX * (y + dimY * z)];
    }

    
    public void setGradient(int x, int y, int z, VoxelGradient value) {
        data[x + dimX * (y + dimY * z)] = value;
    }

    public void setVoxel(int i, VoxelGradient value) {
        data[i] = value;
    }

    public VoxelGradient getVoxel(int i) {
        return data[i];
    }
    
    public int getDimX() {
        return dimX;
    }

    public int getDimY() {
        return dimY;
    }

    public int getDimZ() {
        return dimZ;
    }

    private void compute() {

        // this just initializes all gradients to the vector (0,0,0)
        for (int i=0; i<data.length; i++) {
            data[i] = zero;
        }                
    }
    
    public double getgradientmagnitude(double[] coord) {
        double[] gradient = new double[3];
        double[] coord1 = new double[3];
        double[] coord2 = new double[3];
        
        //X
        coord1[0] = coord[0] + 1;  //x+1
        coord1[1] = coord[1];
        coord1[2] = coord[2];
        
        coord2[0] = coord[0] - 1;  //x-1
        coord2[1] = coord[1];
        coord2[2] = coord[2];
                
        gradient[0] = 0.5 * (RR.getVoxel2(coord1) - RR.getVoxel2(coord2));
        
//y
        coord1[0] = coord[0]-1;
        coord1[1] = coord[1]+1;//y+1
        coord1[2] = coord[2];
        
        coord2[0] = coord[0]+1;
        coord2[1] = coord[1] -1;//y-1
        coord2[2] = coord[2];
        gradient[1] = 0.5 * (RR.getVoxel2(coord1) - RR.getVoxel2(coord2));
        
        //z
        coord1[0] = coord[0];
        coord1[1] = coord[1]-1;
        coord1[2] = coord[2]+1;//z+1
        
        coord2[0] = coord[0];
        coord2[1] = coord[1]+1;
        coord2[2] = coord[2]-1;//z-1
        
        gradient[2] = 0.5 * (RR.getVoxel2(coord1) - RR.getVoxel2(coord2));
        
        return getMagnitude(gradient);
 }
 
   //Calculate gradient Nearest Neighbour
    public double getgradientmagnitude2(double[] coord) {
        double[] gradient = new double[3];
        double[] coord1 = new double[3];
        double[] coord2 = new double[3];
        
        //X
        coord1[0] = coord[0] + 1;  //x+1
        coord1[1] = coord[1];
        coord1[2] = coord[2];
        
        coord2[0] = coord[0] - 1;  //x-1
        coord2[1] = coord[1];
        coord2[2] = coord[2];
                
        gradient[0] = 0.5 * (RR.getVoxel(coord1) - RR.getVoxel(coord2));
        
//y
        coord1[0] = coord[0] - 1;
        coord1[1] = coord[1]+1;//y+1
        coord1[2] = coord[2];
        
        coord2[0] = coord[0] + 1;
        coord2[1] = coord[1] -1;//y-1
        coord2[2] = coord[2];
        gradient[1] = 0.5 * (RR.getVoxel(coord1) - RR.getVoxel(coord2));
        
        //z
        coord1[0] = coord[0];
        coord1[1] = coord[1] - 1;
        coord1[2] = coord[2]+1;//z+1
        
        coord2[0] = coord[0];
        coord2[1] = coord[1] + 1;
        coord2[2] = coord[2]-1;//z-1
        
        gradient[2] = 0.5 * (RR.getVoxel(coord1) - RR.getVoxel(coord2));
        
        return getMagnitude(gradient);
 }
    
    public double getMaxGradientMagnitude() {
        if (maxmag >= 0) {
            return maxmag;
        } else {
            double magnitude = data[0].mag;
            for (int i=0; i<data.length; i++) {
                magnitude = data[i].mag > magnitude ? data[i].mag : magnitude;
            }   
            maxmag = magnitude;
            return magnitude;
        }
    }
    
    private int dimX, dimY, dimZ;
    private VoxelGradient zero = new VoxelGradient();
    VoxelGradient[] data;
    Volume volume;
    double maxmag;
}
