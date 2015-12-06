/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volvis;

import com.jogamp.opengl.util.texture.Texture;
import com.jogamp.opengl.util.texture.awt.AWTTextureIO;
import gui.RaycastRendererPanel;
import gui.TransferFunction2DEditor;
import gui.TransferFunctionEditor;
import java.awt.image.BufferedImage;
import java.awt.Color;
import javax.media.opengl.GL;
import javax.media.opengl.GL2;
import util.TFChangeListener;
import util.VectorMath;
import volume.GradientVolume;
import volume.Volume;
import volume.VoxelGradient;

/**
 *
 * @author michel
 */
public class RaycastRenderer extends Renderer implements TFChangeListener {

    private Volume volume = null;
    private GradientVolume gradients = null;
    RaycastRendererPanel panel;
    TransferFunction tFunc;
    TransferFunctionEditor tfEditor;
    TransferFunction2DEditor tfEditor2D;
    public boolean MIP = false;
    public boolean slicer = false;
    public boolean compositing = false;
    public boolean transfer2D = false;
    
    public RaycastRenderer() {
        panel = new RaycastRendererPanel(this);
        panel.setSpeedLabel("0");
    }

    public void setVolume(Volume vol) {
        System.out.println("Assigning volume");
        volume = vol;

        System.out.println("Computing gradients");
        gradients = new GradientVolume(vol);

        // set up image for storing the resulting rendering
        // the image width and height are equal to the length of the volume diagonal
        int imageSize = (int) Math.floor(Math.sqrt(vol.getDimX() * vol.getDimX() + vol.getDimY() * vol.getDimY()
                + vol.getDimZ() * vol.getDimZ()));
        if (imageSize % 2 != 0) {
            imageSize = imageSize + 1;
        }
        image = new BufferedImage(imageSize, imageSize, BufferedImage.TYPE_INT_ARGB);
        // create a standard TF where lowest intensity maps to black, the highest to white, and opacity increases
        // linearly from 0.0 to 1.0 over the intensity range
        tFunc = new TransferFunction(volume.getMinimum(), volume.getMaximum());
        
        // uncomment this to initialize the TF with good starting values for the orange dataset 
        //tFunc.setTestFunc();
        
        
        tFunc.addTFChangeListener(this);
        tfEditor = new TransferFunctionEditor(tFunc, volume.getHistogram());
        
        tfEditor2D = new TransferFunction2DEditor(volume, gradients);
        tfEditor2D.addTFChangeListener(this);

        System.out.println("Finished initialization of RaycastRenderer");
    }

    public RaycastRendererPanel getPanel() {
        return panel;
    }

    public TransferFunction2DEditor getTF2DPanel() {
        return tfEditor2D;
    }
    
    public TransferFunctionEditor getTFPanel() {
        return tfEditor;
    }
     

    short getVoxel(double[] coord) {

        if (coord[0] < 0 || coord[0] > volume.getDimX() || coord[1] < 0 || coord[1] > volume.getDimY()
                || coord[2] < 0 || coord[2] > volume.getDimZ()) {
            return 0;
        }

        int x = (int) Math.floor(coord[0]);
        int y = (int) Math.floor(coord[1]);
        int z = (int) Math.floor(coord[2]);

        return volume.getVoxel(x, y, z);
    }


    void slicer(double[] viewMatrix) {

        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();

        
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                        + volumeCenter[0];
                pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                        + volumeCenter[1];
                pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                        + volumeCenter[2];

                int val = getVoxel(pixelCoord);
                
                // Map the intensity to a grey value by linear scaling
                voxelColor.r = val/max;
                voxelColor.g = voxelColor.r;
                voxelColor.b = voxelColor.r;
                voxelColor.a = val > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                // Alternatively, apply the transfer function to obtain a color
                // voxelColor = tFunc.getColor(val);
                
                
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
            }
        }

    }


    private void drawBoundingBox(GL2 gl) {
        gl.glPushAttrib(GL2.GL_CURRENT_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glColor4d(1.0, 1.0, 1.0, 1.0);
        gl.glLineWidth(1.5f);
        gl.glEnable(GL.GL_LINE_SMOOTH);
        gl.glHint(GL.GL_LINE_SMOOTH_HINT, GL.GL_NICEST);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glDisable(GL.GL_LINE_SMOOTH);
        gl.glDisable(GL.GL_BLEND);
        gl.glEnable(GL2.GL_LIGHTING);
        gl.glPopAttrib();

    }

    @Override
    public void visualize(GL2 gl) {


        if (volume == null) {
            return;
        }

        drawBoundingBox(gl);

        gl.glGetDoublev(GL2.GL_MODELVIEW_MATRIX, viewMatrix, 0);

        long startTime = System.currentTimeMillis();
        
        if (MIP == true) {
            MIP(viewMatrix);
        }
        if (slicer == true) {
            slicer(viewMatrix);
        }
        if (compositing == true) {
            Compositing(viewMatrix);
        }
        
        long endTime = System.currentTimeMillis();
        double runningTime = (endTime - startTime);
        panel.setSpeedLabel(Double.toString(runningTime));

        Texture texture = AWTTextureIO.newTexture(gl.getGLProfile(), image, false);

        gl.glPushAttrib(GL2.GL_LIGHTING_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        // draw rendered image as a billboard texture
        texture.enable(gl);
        texture.bind(gl);
        double halfWidth = image.getWidth() / 2.0;
        gl.glPushMatrix();
        gl.glLoadIdentity();
        gl.glBegin(GL2.GL_QUADS);
        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        gl.glTexCoord2d(0.0, 0.0);
        gl.glVertex3d(-halfWidth, -halfWidth, 0.0);
        gl.glTexCoord2d(0.0, 1.0);
        gl.glVertex3d(-halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 1.0);
        gl.glVertex3d(halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 0.0);
        gl.glVertex3d(halfWidth, -halfWidth, 0.0);
        gl.glEnd();
        texture.disable(gl);
        texture.destroy(gl);
        gl.glPopMatrix();

        gl.glPopAttrib();


        if (gl.glGetError() > 0) {
            System.out.println("some OpenGL error: " + gl.glGetError());
        }

    }
    private BufferedImage image;
    private double[] viewMatrix = new double[4 * 4];

    @Override
    public void changed() {
        for (int i=0; i < listeners.size(); i++) {
            listeners.get(i).changed();
        }
    }
    
    // get a voxel from the volume data by trilinear interpolation    
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
double getmagnitude(double[] vector){
   return Math.sqrt(vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2]);
}   
 
  //Calculate gradient trilinear interpolation
 double getgradientmagnitude(double[] coord) {
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
                
        gradient[0] = 0.5 * (getVoxel2(coord1) - getVoxel2(coord2));
        
//y
        coord1[0] = coord[0]-1;
        coord1[1] = coord[1]+1;//y+1
        coord1[2] = coord[2];
        
        coord2[0] = coord[0]+1;
        coord2[1] = coord[1] -1;//y-1
        coord2[2] = coord[2];
        gradient[1] = 0.5 * (getVoxel2(coord1) - getVoxel2(coord2));
        
        //z
        coord1[0] = coord[0];
        coord1[1] = coord[1]-1;
        coord1[2] = coord[2]+1;//z+1
        
        coord2[0] = coord[0];
        coord2[1] = coord[1]+1;
        coord2[2] = coord[2]-1;//z-1
        
        gradient[2] = 0.5 * (getVoxel2(coord1) - getVoxel2(coord2));
        
        return getmagnitude(gradient);
 }
 
   //Calculate gradient Nearest Neighbour
 double getgradientmagnitude2(double[] coord) {
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
                
        gradient[0] = 0.5 * (getVoxel(coord1) - getVoxel(coord2));
        
//y
        coord1[0] = coord[0] - 1;
        coord1[1] = coord[1]+1;//y+1
        coord1[2] = coord[2];
        
        coord2[0] = coord[0] + 1;
        coord2[1] = coord[1] -1;//y-1
        coord2[2] = coord[2];
        gradient[1] = 0.5 * (getVoxel(coord1) - getVoxel(coord2));
        
        //z
        coord1[0] = coord[0];
        coord1[1] = coord[1] - 1;
        coord1[2] = coord[2]+1;//z+1
        
        coord2[0] = coord[0];
        coord2[1] = coord[1] + 1;
        coord2[2] = coord[2]-1;//z-1
        
        gradient[2] = 0.5 * (getVoxel(coord1) - getVoxel(coord2));
        
        return getmagnitude(gradient);
 }
 
    double CalculateOpacity(double fx, double alphan, double fvn, double alphanplus1, double fvnplus1,double gradmagnitude){
    double opacity;
    if ((fx >= fvn) && (fx<=fvnplus1)){
        
        opacity = gradmagnitude * ( (alphanplus1*(fx - fvn)/(fvnplus1 - fvn))
                    + (alphan*(fvnplus1 - fx)/(fvnplus1-fvn)));
        return opacity;
    }else{
        return 0;
    }
 }
 
     public void MIP(double[] viewMatrix) {
        //use lower resolution for responsiveness
        int scale = 2;
        int imageheight = image.getHeight()/scale;
        int imagewidth = image.getWidth()/scale;
        int imageCenter = image.getWidth()/2;   
        //System.out.println("j0 = " + (imageheight/2 - imageCenter) + " jn = " + (imageheight/2 + imageCenter));
        // clear image and put image in the center
        for (int j = imageCenter- imageheight/2; j < imageheight/2 + imageCenter; j++) {
            for (int i = imageCenter- imagewidth/2; i < imagewidth/2 + imageCenter; i++) {
              image.setRGB(i, j, 0);
            }
        }
        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        
        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);
        int DimZ = volume.getDimZ();
        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        int sampling_distance = 2;
        for (int j = imageCenter - imageheight/2; j < imageheight/2 + imageCenter; j++) {
            for (int i = imageCenter - imagewidth/2; i < imagewidth/2 + imageCenter; i++) {
                //Cast a ray!! 
                //System.out.println("j = " + (j-imageheight/2 - imageCenter)*scale + " i = " + (i-imagewidth/2 - imageCenter)*scale);
                int maxvalue = 0;
                for (int k = 0; k < DimZ;k++){
                   //(i - (imageCenter - imagewidth/2)) shift i = 0
                   pixelCoord[0] = uVec[0] * ((i - (imageCenter - imagewidth/2))*scale - imageCenter) + vVec[0] * ((j - (imageCenter - imageheight/2))*scale - imageCenter)
                        + volumeCenter[0] + k * sampling_distance * viewVec[0];
                   pixelCoord[1] = uVec[1] * ((i - (imageCenter - imagewidth/2))*scale - imageCenter) + vVec[1] * ((j - (imageCenter - imageheight/2))*scale - imageCenter)
                        + volumeCenter[1] + k * sampling_distance * viewVec[1];
                   pixelCoord[2] = uVec[2] * ((i - (imageCenter - imagewidth/2))*scale - imageCenter) + vVec[2] * ((j - (imageCenter - imageheight/2))*scale- imageCenter)
                        + volumeCenter[2] + k * sampling_distance * viewVec[2];
                   int val = (int) getVoxel2(pixelCoord);
                   if(val > maxvalue ) maxvalue = val; //MIP
                   if(pixelCoord[2]>= DimZ) break;
                }
                // Apply the transfer function to obtain a color
                //System.out.println("max value = " + maxvalue);
                TFColor voxelColor = tFunc.getColor(maxvalue);
				
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
            }
        }
    }
    
     void Compositing (double[] viewMatrix) {
        //Back to front
        //use lower resolution for responsiveness
        int scale = 3;
        int imageheight = image.getHeight()/scale;
        int imagewidth = image.getWidth()/scale;
        int imageCenter = image.getWidth()/2;   
        //System.out.println("j0 = " + (imageheight/2 - imageCenter) + " jn = " + (imageheight/2 + imageCenter));
        // clear image and put image in the center
        for (int j = imageCenter- imageheight/2; j < imageheight/2 + imageCenter; j++) {
            for (int i = imageCenter- imagewidth/2; i < imagewidth/2 + imageCenter; i++) {
              image.setRGB(i, j, 0);
            }
        }

        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);
        int DimZ = volume.getDimZ();
        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        int sampling_distance = 4;
        for (int j = imageCenter - imageheight/2; j < imageheight/2 + imageCenter; j++) {
            for (int i = imageCenter - imagewidth/2; i < imagewidth/2 + imageCenter; i++) {
                //Cast a ray!! 
                //System.out.println("j = " + (j-imageheight/2 - imageCenter)*scale + " i = " + (i-imagewidth/2 - imageCenter)*scale);
                double accumulatedColorR = 0;
                double accumulatedColorG = 0;
                double accumulatedColorB = 0;
                double accumulatedopacity = 0;
                for (int k = 0; k < DimZ;k++){

                   //(i - (imageCenter - imagewidth/2)) shift i = 0
                   pixelCoord[0] = uVec[0] * ((i - (imageCenter - imagewidth/2))*scale - imageCenter) + vVec[0] * ((j - (imageCenter - imageheight/2))*scale - imageCenter)
                        + volumeCenter[0] + k * sampling_distance * viewVec[0];
                   pixelCoord[1] = uVec[1] * ((i - (imageCenter - imagewidth/2))*scale - imageCenter) + vVec[1] * ((j - (imageCenter - imageheight/2))*scale - imageCenter)
                        + volumeCenter[1] + k * sampling_distance * viewVec[1];
                   pixelCoord[2] = uVec[2] * ((i - (imageCenter - imagewidth/2))*scale - imageCenter) + vVec[2] * ((j - (imageCenter - imageheight/2))*scale- imageCenter)
                        + volumeCenter[2] + k * sampling_distance * viewVec[2];
                   int val = (int) getVoxel2(pixelCoord);
                   
                   // Apply the transfer function to obtain a color
                   TFColor voxelColor = tFunc.getColor(val);
                   
                   // Levoy's Compositing front to end
                   if (voxelColor.a > 0) {
                     accumulatedColorR = voxelColor.r * voxelColor.a + accumulatedColorR * (1-voxelColor.a) ;
                     accumulatedColorG = voxelColor.g * voxelColor.a+ accumulatedColorG * (1-voxelColor.a) ;
                     accumulatedColorB = voxelColor.b * voxelColor.a+ accumulatedColorB * (1-voxelColor.a) ;
                     accumulatedopacity = accumulatedopacity * (1 - voxelColor.a) + voxelColor.a;
                   }
                   
                   if(pixelCoord[2]>= DimZ) break;
                   if(accumulatedopacity>=1) break;
                }
                	   
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = accumulatedopacity <= 1.0 ? (int) Math.floor(accumulatedopacity * 255) : 255;
                int c_red = accumulatedColorR <= 1.0 ? (int) Math.floor(accumulatedColorR * 255) : 255;
                int c_green = accumulatedColorG <= 1.0 ? (int) Math.floor(accumulatedColorG * 255) : 255;
                int c_blue = accumulatedColorB <= 1.0 ? (int) Math.floor(accumulatedColorB * 255) : 255;

                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
      
                image.setRGB(i, j, pixelColor);
            }
        }
    }
           void Region (double[] viewMatrix,double alphan, double fvn, double alphanplus1, double fvnplus1) {
		// Region based surface using local normalization one per ray
		
        //use lower resolution for responsiveness
        int scale = 1;
        int imageheight = image.getHeight()/scale;
        int imagewidth = image.getWidth()/scale;
        int imageCenter = image.getWidth()/2;   
        //System.out.println("j0 = " + (imageheight/2 - imageCenter) + " jn = " + (imageheight/2 + imageCenter));
        // clear image and put image in the center
        for (int j = imageCenter- imageheight/2; j < imageheight/2 + imageCenter; j++) {
            for (int i = imageCenter- imagewidth/2; i < imagewidth/2 + imageCenter; i++) {
              image.setRGB(i, j, 0);
            }
        }

        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);
        int DimZ = volume.getDimZ();
        // sample on a plane through the origin of the volume data
        //double max = volume.getMaximum();
        int sampling_distance = 4;
        for (int j = imageCenter - imageheight/2; j < imageheight/2 + imageCenter; j++) {
            for (int i = imageCenter - imagewidth/2; i < imagewidth/2 + imageCenter; i++) {
                //Cast a ray!! 
                //System.out.println("j = " + (j-imageheight/2 - imageCenter)*scale + " i = " + (i-imagewidth/2 - imageCenter)*scale);
                double accumulatedColorR = 0;
                double accumulatedColorG = 0;
                double accumulatedColorB = 0;
                double accumulatedopacity = 0;
                double totalGradient = 0;
                int[] buffer = new int[DimZ];
                double[] gradientbuffer = new double[DimZ];
                int index=0;
                double gradmax = 0;
                double gradmin = 999;
                for (int k = 0; k < DimZ;k++){
                   //(i - (imageCenter - imagewidth/2)) shift i = 0
                   pixelCoord[0] = uVec[0] * ((i - (imageCenter - imagewidth/2))*scale - imageCenter) + vVec[0] * ((j - (imageCenter - imageheight/2))*scale - imageCenter)
                        + volumeCenter[0] + k * sampling_distance * viewVec[0];
                   pixelCoord[1] = uVec[1] * ((i - (imageCenter - imagewidth/2))*scale - imageCenter) + vVec[1] * ((j - (imageCenter - imageheight/2))*scale - imageCenter)
                        + volumeCenter[1] + k * sampling_distance * viewVec[1];
                   pixelCoord[2] = uVec[2] * ((i - (imageCenter - imagewidth/2))*scale - imageCenter) + vVec[2] * ((j - (imageCenter - imageheight/2))*scale- imageCenter)
                        + volumeCenter[2] + k * sampling_distance * viewVec[2];
                   int val = (int) getVoxel2(pixelCoord);
                   index=k;
                   buffer[index] = val;
                   double gradientmagnitude = getgradientmagnitude(pixelCoord);
                   gradientbuffer[index] = gradientmagnitude;
                   totalGradient = totalGradient + gradientmagnitude;
                   //Calculate opacity in the buffer
                   if(gradmax<gradientmagnitude) gradmax = gradientmagnitude;
                   if(gradmin>gradientmagnitude) gradmin = gradientmagnitude;
                   if(pixelCoord[2]>= DimZ) break;
                }
                
                for (int k = 0; k < index;k++){
                   // Apply the transfer function to obtain a color
                   TFColor voxelColor = tFunc.getColor(buffer[k]);
                   
                   //Get normalized gradient magnitude (Feature scaling)
                   double gradmagnitude = (gradientbuffer[k] - gradmin) / (gradmax - gradmin);                   
                 
                    voxelColor.a = 3*CalculateOpacity(buffer[k],alphan,fvn,alphanplus1,fvnplus1,gradmagnitude);

                   // Levoy's Compositing back to front
				   
                   if (voxelColor.a > 0) {
                     accumulatedColorR = voxelColor.r * voxelColor.a + accumulatedColorR * (1-voxelColor.a) ;
                     accumulatedColorG = voxelColor.g * voxelColor.a + accumulatedColorG * (1-voxelColor.a) ;
                     accumulatedColorB = voxelColor.b * voxelColor.a + accumulatedColorB * (1-voxelColor.a) ;
                     accumulatedopacity = accumulatedopacity * (1 - voxelColor.a) + voxelColor.a;
                   }
                   if(accumulatedopacity>=1) break;				   
		
                }
                // BufferedImage expects a pixel color packed as ARGB in an int
                //System.out.println("accumulatedrgba: " + accumulatedopacity + " " + accumulatedColorR + " " + accumulatedColorG + " " + accumulatedColorB + " ");
                int c_alpha = accumulatedopacity <= 1.0 ? (int) Math.floor(accumulatedopacity * 255) : 255;
                int c_red = accumulatedColorR <= 1.0 ? (int) Math.floor(accumulatedColorR * 255) : 255;
                int c_green = accumulatedColorG <= 1.0 ? (int) Math.floor(accumulatedColorG * 255) : 255;
                int c_blue = accumulatedColorB <= 1.0 ? (int) Math.floor(accumulatedColorB * 255) : 255;
                //System.out.println("argb: " + c_alpha + " " + c_red + " " + c_green + " " + c_blue + " ");
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
      
                image.setRGB(i, j, pixelColor);
            }
        }
    }

	 void Region2 (double[] viewMatrix,double alphan, double fvn, double alphanplus1, double fvnplus1) {
		// Region based surface using global normalization 
		
        //use lower resolution for responsiveness
        int scale = 2;
        int imageheight = image.getHeight()/scale;
        int imagewidth = image.getWidth()/scale;
        int imageCenter = image.getWidth()/2;   
		int DimZ = volume.getDimZ();
		double gradmax = 0;
        double gradmin = 999;
		
		// vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);
        int sampling_distance = 4;
		
        // image is square
        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);
        // clear image and put image in the center
        int countTot = 0;
        double totalGradient = 0;
        double[] grad = new double[DimZ];
        for (int j = imageCenter- imageheight/2; j < imageheight/2 + imageCenter; j++) {
            for (int i = imageCenter- imagewidth/2; i < imagewidth/2 + imageCenter; i++) {
              image.setRGB(i, j, 0);
				// Find global gradient minimum and maximum
			     for (int k = 0; k < DimZ;k++){
					pixelCoord[0] = uVec[0] * ((i - (imageCenter - imagewidth/2))*scale - imageCenter) + vVec[0] * ((j - (imageCenter - imageheight/2))*scale - imageCenter) + volumeCenter[0] + k * sampling_distance * viewVec[0];
                   pixelCoord[1] = uVec[1] * ((i - (imageCenter - imagewidth/2))*scale - imageCenter) + vVec[1] * ((j - (imageCenter - imageheight/2))*scale - imageCenter)
                        + volumeCenter[1] + k * sampling_distance * viewVec[1];
                   pixelCoord[2] = uVec[2] * ((i - (imageCenter - imagewidth/2))*scale - imageCenter) + vVec[2] * ((j - (imageCenter - imageheight/2))*scale- imageCenter)
                        + volumeCenter[2] + k * sampling_distance * viewVec[2];
                   if(gradmax<getgradientmagnitude(pixelCoord)) gradmax = getgradientmagnitude(pixelCoord);
                   if(gradmin>getgradientmagnitude(pixelCoord)) gradmin = getgradientmagnitude(pixelCoord);                 
                   if(pixelCoord[2]>= DimZ) break;
                }
            }
        }

        // sample on a plane through the origin of the volume data
        //double max = volume.getMaximum();
        
        for (int j = imageCenter - imageheight/2; j < imageheight/2 + imageCenter; j++) {
            for (int i = imageCenter - imagewidth/2; i < imagewidth/2 + imageCenter; i++) {
                //Cast a ray!! 
                //System.out.println("j = " + (j-imageheight/2 - imageCenter)*scale + " i = " + (i-imagewidth/2 - imageCenter)*scale);
                double accumulatedColorR = 0;
                double accumulatedColorG = 0;
                double accumulatedColorB = 0;
                double accumulatedopacity = 0;
                int[] buffer = new int[DimZ];
                double[] gradientbuffer = new double[DimZ];
                int index=0;
                for (int k = 0; k < DimZ;k++){
                   //(i - (imageCenter - imagewidth/2)) shift i = 0
                   pixelCoord[0] = uVec[0] * ((i - (imageCenter - imagewidth/2))*scale - imageCenter) + vVec[0] * ((j - (imageCenter - imageheight/2))*scale - imageCenter)
                        + volumeCenter[0] + k * sampling_distance * viewVec[0];
                   pixelCoord[1] = uVec[1] * ((i - (imageCenter - imagewidth/2))*scale - imageCenter) + vVec[1] * ((j - (imageCenter - imageheight/2))*scale - imageCenter)
                        + volumeCenter[1] + k * sampling_distance * viewVec[1];
                   pixelCoord[2] = uVec[2] * ((i - (imageCenter - imagewidth/2))*scale - imageCenter) + vVec[2] * ((j - (imageCenter - imageheight/2))*scale- imageCenter)
                        + volumeCenter[2] + k * sampling_distance * viewVec[2];
                   int val = (int) getVoxel2(pixelCoord);
                   index=k;
                   buffer[index] = val;
                   double gradientmagnitude = getgradientmagnitude(pixelCoord);
                   gradientbuffer[index] = gradientmagnitude;
                   //Calculate opacity in the buffer
                   if(pixelCoord[2]>= DimZ) break;
                }
                
                for (int k = 0; k < index;k++){
                   // Apply the transfer function to obtain a color
                   TFColor voxelColor = tFunc.getColor(buffer[k]);
                   
                   //Get normalized gradient magnitude (Feature scaling)
                  double gradmagnitude = (gradientbuffer[k] - gradmin) / (gradmax - gradmin);
                   

                    voxelColor.a = 3 * CalculateOpacity(buffer[k],alphan,fvn,alphanplus1,fvnplus1,gradmagnitude);

                   // Levoy's Compositing back to front
                   if (voxelColor.a > 0) {
                     accumulatedColorR = voxelColor.r * voxelColor.a + accumulatedColorR * (1-voxelColor.a) ;
                     accumulatedColorG = voxelColor.g * voxelColor.a+ accumulatedColorG * (1-voxelColor.a) ;
                     accumulatedColorB = voxelColor.b * voxelColor.a+ accumulatedColorB * (1-voxelColor.a) ;
                     accumulatedopacity = accumulatedopacity * (1 - voxelColor.a) + voxelColor.a;
                   }
                   if(accumulatedopacity>=1) break;
                }
                // BufferedImage expects a pixel color packed as ARGB in an int
                //System.out.println("accumulatedrgba: " + accumulatedopacity + " " + accumulatedColorR + " " + accumulatedColorG + " " + accumulatedColorB + " ");
                int c_alpha = accumulatedopacity <= 1.0 ? (int) Math.floor(accumulatedopacity * 255) : 255;
                int c_red = accumulatedColorR <= 1.0 ? (int) Math.floor(accumulatedColorR * 255) : 255;
                int c_green = accumulatedColorG <= 1.0 ? (int) Math.floor(accumulatedColorG * 255) : 255;
                int c_blue = accumulatedColorB <= 1.0 ? (int) Math.floor(accumulatedColorB * 255) : 255;
                //System.out.println("argb: " + c_alpha + " " + c_red + " " + c_green + " " + c_blue + " ");
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
      
                image.setRGB(i, j, pixelColor);
            }
        }
    }
}
