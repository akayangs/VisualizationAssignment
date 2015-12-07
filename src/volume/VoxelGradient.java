/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volume;

/**
 *
 * @author michel
 */
public class VoxelGradient {

    public double x, y, z;
    public double mag;
    
    public VoxelGradient() {
        x = y = z = mag = 0.0f;
    }
    
    public VoxelGradient(double gx, double gy, double gz) {
        x = gx;
        y = gy;
        z = gz;
        mag = (double) Math.sqrt(x*x + y*y + z*z);
    }
    
}
