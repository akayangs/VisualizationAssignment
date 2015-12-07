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
   
    public VoxelGradient getGradient(double x, double y, double z) {
        double[] gradient = new double[3];
        double[] coord1 = new double[3];
        double[] coord2 = new double[3];
        
        gradient[0] = x;
        gradient[1] = y;
        gradient[2] = z;
        
        //X
        coord1[0] = dimX+ 1;  //x+1
        coord1[1] = dimY;
        coord1[2] = dimZ;
        
        coord2[0] = dimX - 1;  //x-1
        coord2[1] = dimY;
        coord2[2] = dimZ;
                
        x = 0.5 * (getVoxel2(coord1) - getVoxel2(coord2));
        
//y
        coord1[0] = dimX-1;
        coord1[1] = dimY+1;//y+1
        coord1[2] = dimZ;
        
        coord2[0] = dimX+1;
        coord2[1] = dimY -1;//y-1
        coord2[2] = dimZ;
        y = 0.5 * (getVoxel2(coord1) - getVoxel2(coord2));
        
        //z
        coord1[0] = dimX;
        coord1[1] = dimY-1;
        coord1[2] = dimZ+1;//z+1
        
        coord2[0] = dimX;
        coord2[1] = dimY+1;
        coord2[2] = dimZ-1;//z-1
        
        z = 0.5 * (getVoxel2(coord1) - getVoxel2(coord2));
      

        return getMagnitude(gradient);
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
double getVoxel2(double[] coord) {
        
        //XYZ coordinate
        double x = coord[0];
        double y = coord[1];
        double z = coord[2];

        // Get the box in xyz
        int x0 = (int) Math.floor(coord[0]);
        int y0 = (int) Math.floor(coord[1]);
        int z0 = (int) Math.floor(coord[2]);
        int x1 = x0+1;
        int y1 = y0+1;
        int z1 = z0+1;
        //System.out.println("x = " + x + " y = " + y + " z = " +z);
        
        if ((x0 >= 0) && (x1 < volume.getDimX()) && (y0 >= 0) && (y1 < volume.getDimY())
                && (z0 >= 0) && (z1 < volume.getDimZ())) {
            double T_voxel =  volume.getVoxel(x0,y0,z0) * (x1-x) * (y1-y) * (z1-z) +
                            volume.getVoxel(x1,y0,z0) * (x-x0) * (y1-y) * (z1-z) + 
                            volume.getVoxel(x0,y1,z0) * (x1 - x) * (y-y0) * (z1 - z) +
                            volume.getVoxel(x0,y0,z1) * (x1 - x) * (y1 - y) * (z-z0) +
                            volume.getVoxel(x1,y0,z1) * (x-x0) * (y1 - y) * (z-z0) +
                            volume.getVoxel(x0,y1,z1) * (x1 - x) * (y-y0) * (z-z0) +
                            volume.getVoxel(x1,y1,z0) * (x-x0) * (y-y0) * (z1 - z) +
                            volume.getVoxel(x1,y1,z1) * (x-x0) * (y-y0) * (z-z0);
            		//System.out.println("T_VOXEL = " + T_voxel);					
            return T_voxel;
        } else {
            return 0;
        }
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
