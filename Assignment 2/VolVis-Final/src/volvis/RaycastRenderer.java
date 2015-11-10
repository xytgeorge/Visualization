/*/
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volvis;

import com.jogamp.opengl.util.texture.Texture;
import com.jogamp.opengl.util.texture.awt.AWTTextureIO;
import gui.RaycastRendererPanel;
import gui.TransferFunctionEditor;
import java.awt.image.BufferedImage;
import javax.media.opengl.GL2;
import util.TFChangeListener;
import util.VectorMath;
import volume.Volume;
import java.awt.AWTException;
import javax.media.opengl.GLContext;



/**
 *
 * @author michel
 */
public class RaycastRenderer extends Renderer implements TFChangeListener {

    private Volume volume = null;
    RaycastRendererPanel panel;
    TransferFunction tFunc;
    TransferFunctionEditor tfEditor;
    //private int step = 5;
    public static int activeRenderer = 0;//using as flags
    
    private GL2 glCanvas;
    private double scalar = 1;
    private Visualization visualization;
    
    public void setVisualization(Visualization visualization) {
        this.visualization = visualization;
    }
    
    public void update(){
        this.visualization.update();
    }

    public RaycastRenderer() {
        panel = new RaycastRendererPanel(this);
        panel.setSpeedLabel("0");
    }

    public void setVolume(Volume vol) {
        volume = vol;

        // set up image for storing the resulting rendering
        // the image width and height are equal to the length of the volume diagonal
        int imageSize = (int) Math.floor(Math.sqrt(vol.getDimX() * vol.getDimX() + vol.getDimY() * vol.getDimY()
                + vol.getDimZ() * vol.getDimZ()));
        if (imageSize % 2 != 0) {
            imageSize = imageSize + 1;
        }
        image = new BufferedImage(imageSize, imageSize, BufferedImage.TYPE_INT_ARGB);
        tFunc = new TransferFunction(volume.getMinimum(), volume.getMaximum(), vol);
        tFunc.addTFChangeListener(this);
        tfEditor = new TransferFunctionEditor(tFunc, volume.getHistogram(), this);
        panel.setTransferFunctionEditor(tfEditor);
    }
    
    public void resetTransferFunction(Volume vol){
        tFunc = new TransferFunction(volume.getMinimum(), volume.getMaximum(), vol);
        tFunc.addTFChangeListener(this);
        tfEditor = new TransferFunctionEditor(tFunc, volume.getHistogram(),this);
        panel.setTransferFunctionEditor(tfEditor);
    }

    @Override
    public void changed() {
        for (int i = 0; i < listeners.size(); i++) {
            listeners.get(i).changed();
        }
    }

    public RaycastRendererPanel getPanel() {
        return panel;
    }
    
    short getVoxel(double xi,double yi,double zi) {

        int x = (int) xi;
        int y = (int) yi;
        int z = (int) zi;

        if ((x >= 0) && (x < volume.getDimX()) && (y >= 0) && (y < volume.getDimY())
                && (z >= 0) && (z < volume.getDimZ())) {
            return volume.getVoxel(x, y, z);
        } else {
            return 0;
        }
    }
    
    short getVoxel(double[] coord) {

        double x = coord[0];
        double y = coord[1];
        double z = coord[2];

        if ((x >= 0) && (x < volume.getDimX()) && (y >= 0) && (y < volume.getDimY())
                && (z >= 0) && (z < volume.getDimZ())) {

            int x0 = (int) Math.floor(x);
            int y0 = (int) Math.floor(y);
            int z0 = (int) Math.floor(z);

            int x1 = (int) Math.ceil(x);
            int y1 = (int) Math.ceil(y);
            int z1 = (int) Math.ceil(z);
            // interplotate along x-axis
            int xd = (int) ((x - x0) / (x1 - x0));
            int yd = (int) ((y - y0) / (y1 - y0));
            int zd = (int) ((z - z0) / (z1 - z0));
            // interploate along y-axis
            try {
                double c00 = volume.getVoxel(x0, y0, z0) * (1 - xd) + volume.getVoxel(x1, y0, z0) * xd;
                double c10 = volume.getVoxel(x0, y1, z0) * (1 - xd) + volume.getVoxel(x1, y1, z0) * xd;
                double c01 = volume.getVoxel(x0, y0, z1) * (1 - xd) + volume.getVoxel(x1, y0, z1) * xd;
                double c11 = volume.getVoxel(x0, y1, z1) * (1 - xd) + volume.getVoxel(x1, y1, z1) * xd;

                double c0 = c00 * (1 - yd) + c10 * yd;
                double c1 = c01 * (1 - yd) + c11 * yd;
                // interploate along z-axis
                return (short) (c0 * (1 - zd) + c1 * zd);
            } catch (Exception e) {
                return 0;
            }
        } else {
            return 0;
        }
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
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                        + volumeCenter[0];
                pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                        + volumeCenter[1];
                pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                        + volumeCenter[2];

                int val = getVoxel(pixelCoord);
                // Apply the transfer function to obtain a color
                TFColor voxelColor = tFunc.getColor(val);

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

    double RB_alpha(double val) {
        int fn = (int)Math.floor(val);
        int fn_1 = (int)Math.ceil(val);
        double alpha = 0.0;
        TFColor cn = tFunc.getColor(fn);
        TFColor cn_1 = tFunc.getColor(fn_1);
        if((val >= fn) && (val <= fn_1)) {
            alpha =  (cn_1.a * ((val - fn) / (fn_1 - fn)) + 
                cn.a * ((fn_1 - val) / (fn_1 - fn)));
        }
        else {alpha = 0;}
        return alpha;
    }
    
    void MIP(double[] viewMatrix, GL2 gl) {

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
        int imageSize = (int) Math.floor(Math.sqrt(volume.getDimX() * volume.getDimX() 
                        + volume.getDimY() * volume.getDimY() + volume.getDimZ() * volume.getDimZ()))/2;
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
            int maxValue = Integer.MIN_VALUE;

                for (int k = - imageSize; k < imageSize; k = k + step)
                {
                    pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                        + volumeCenter[0] + viewVec[0] * k;
                    pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                        + volumeCenter[1] + viewVec[1] * k;
                    pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                        + volumeCenter[2] + viewVec[2] * k;
                    
                    int val = getVoxel(pixelCoord);
                    
                    if(val > maxValue)
                    {
                        maxValue = val;
                    }
                }

                // Apply the transfer function to obtain a color
                TFColor voxelColor = tFunc.getColor(maxValue);
               // voxelColor.a = RB_alpha(maxValue);

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
    
    void backFront(double[] viewMatrix,GL2 gl) {

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
        int c_alpha = 0;
        int c_red = 0; 
        int c_green = 0;
        int c_blue = 0;
        double s_alpha = 0;
        double s_red = 0;
        double s_green = 0;
        double s_blue = 0;
        int k=0;
        int pixelColor;
        int val=-1;
        double alpha=0;
        TFColor voxelColor;
        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        int imageSize = (int) Math.floor(Math.sqrt(volume.getDimX() * volume.getDimX() 
                        + volume.getDimY() * volume.getDimY() + volume.getDimZ() * volume.getDimZ()))/2;
 
            for (int j = 0; j < image.getHeight(); j++) {
                for (int i = 0; i < image.getWidth(); i++) {

                     s_alpha = 0;
                     s_red = 0; 
                     s_green = 0;
                     s_blue = 0;
               
                    for (k = - imageSize; k < imageSize; k = k + step){
                        
                   
                        pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                        + volumeCenter[0] + viewVec[0] * k;
                        
                        pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                        + volumeCenter[1] + viewVec[1] * k;
                        
                        pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                        + volumeCenter[2] + viewVec[2] * k;
                        
                        val = getVoxel(pixelCoord);
                        
                        voxelColor = tFunc.getColor(val);
                        
                        alpha = voxelColor.a;
                        
                        s_alpha = s_alpha + ( 1 - s_alpha ) * alpha;
                        s_red = s_red * ( 1 - alpha ) + alpha * voxelColor.r;
                        s_green = s_green * ( 1 - alpha ) + alpha * voxelColor.g;
                        s_blue = s_blue * ( 1 - alpha ) + alpha * voxelColor.b;
                        
                        //s_alpha = RB_alpha(val);
                        
                    }
                         // BufferedImage expects a pixel color packed as ARGB in an int
       
                        c_alpha = (s_alpha) <= 1.0 ? (int) Math.floor(s_alpha * 255) : 255;
                        c_red = (s_red) <= 1.0 ? (int) Math.floor(s_red * 255) : 255;
                        c_green = (s_green) <= 1.0 ? (int) Math.floor(s_green * 255) : 255;
                        c_blue = (s_blue <= 1.0) ? (int) Math.floor(s_blue * 255) : 255;
                        pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                       
                       image.setRGB(i, j, pixelColor);
 
                }
            }
    }
    
    void frontBack(double[] viewMatrix,GL2 gl) {

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
        int c_alpha = 0;
        int c_red = 0; 
        int c_green = 0;
        int c_blue = 0;
        double s_alpha = 0;
        double s_red = 0;
        double s_green = 0;
        double s_blue = 0;
        double alpha = 0;
        int k=0;
        int pixelColor;
        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        int imageSize = (int) Math.floor(Math.sqrt(volume.getDimX() * volume.getDimX() 
                        + volume.getDimY() * volume.getDimY() + volume.getDimZ() * volume.getDimZ()))/2;
 
            for (int j = 0; j < image.getHeight(); j++) {
                for (int i = 0; i < image.getWidth(); i++) {

                     s_alpha = 0;
                     s_red = 0; 
                     s_green = 0;
                     s_blue = 0;
                    
                     for (k = - imageSize; k < imageSize; k = k + step){

                        pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                        + volumeCenter[0] + viewVec[0] * k;
                        pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                        + volumeCenter[1] + viewVec[1] * k;
                        pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                        + volumeCenter[2] + viewVec[2] * k;

                        int val = getVoxel(pixelCoord);
                        TFColor voxelColor = tFunc.getColor(val);

                        alpha = voxelColor.a;
                        
                        s_alpha = s_alpha + ( 1 - s_alpha ) * voxelColor.a;
                        s_red = s_red + voxelColor.r * ( 1 - s_alpha ) * alpha;
                        s_green = s_green + voxelColor.g * ( 1 - s_alpha ) * alpha;
                        s_blue = s_blue + voxelColor.b * ( 1 - s_alpha ) * alpha;
                     
                    }
                         // BufferedImage expects a pixel color packed as ARGB in an int
       
                        c_alpha = (s_alpha) <= 1.0 ? (int) Math.floor(s_alpha * 255) : 255;
                        c_red = (s_red) <= 1.0 ? (int) Math.floor(s_red * 255) : 255;
                        c_green = (s_green) <= 1.0 ? (int) Math.floor(s_green * 255) : 255;
                        c_blue = (s_blue <= 1.0) ? (int) Math.floor(s_blue * 255) : 255;
                        pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                        image.setRGB(i, j, pixelColor);
                    
                }
            }
    }
    
    void OpacityWeighting(double[] viewMatrix, GL2 gl) {

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
        int c_alpha = 0;
        int c_red = 0; 
        int c_green = 0;
        int c_blue = 0;
        double s_alpha = 0;
        double s_red = 0;
        double s_green = 0;
        double s_blue = 0;
        double s = 0;
        int k = 0;
        int pixelColor;
        int val = -1;
        double is = 0;
        double js = 0;
        double ks = 0;
        double x = 0;
        double y = 0;
        double z = 0;
        double alpha = 0;
        double deltaf = 0;
        double deltafMax = 0;
        double deltafNormalized = 0;
        TFColor voxelColor;
        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        int imageSize = (int) Math.floor(Math.sqrt(volume.getDimX() * volume.getDimX() 
                        + volume.getDimY() * volume.getDimY() + volume.getDimZ() * volume.getDimZ()))/2;
 
            for (int j = 0; j < image.getHeight(); j++) {
                for (int i = 0; i < image.getWidth(); i++) {

                     s_alpha = 0;
                     s_red = 0; 
                     s_green = 0;
                     s_blue = 0;
                     is=0;
                     js=0;
                     ks=0;
                     
                     for (k = - imageSize; k < imageSize; k = k + step){
                        
                        is=uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                            + volumeCenter[0] + viewVec[0] * k;
                        pixelCoord[0] = is;
                        
                        js=uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                            + volumeCenter[1] + viewVec[1] * k;
                        pixelCoord[1] = js;
                        
                        ks=uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                            + volumeCenter[2] + viewVec[2] * k;
                        pixelCoord[2] = ks;
                        
                        val = getVoxel(pixelCoord);
                        voxelColor = tFunc.getColor(val);
                        
                        x=(tFunc.getColor(getVoxel(is+1,js,ks)).a-tFunc.getColor(getVoxel(is-1,js,ks)).a)/2;
                        y=(tFunc.getColor(getVoxel(is,js+1,ks)).a-tFunc.getColor(getVoxel(is,js-1,ks)).a)/2;
                        z=(tFunc.getColor(getVoxel(is,js,ks+1)).a-tFunc.getColor(getVoxel(is,js,ks-1)).a)/2;
                      
                        deltaf = Math.sqrt( x * x + y * y + z * z);
                        
                        if(deltaf > deltafMax)
                            deltafMax = deltaf;
                     }
                     
                     for (k = - imageSize; k < imageSize; k = k + step){
                        
                        is=uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                            + volumeCenter[0] + viewVec[0] * k;
                        pixelCoord[0] = is;
                        
                        js=uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                            + volumeCenter[1] + viewVec[1] * k;
                        pixelCoord[1] = js;
                        
                        ks=uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                            + volumeCenter[2] + viewVec[2] * k;
                        pixelCoord[2] = ks;
                        
                        val = getVoxel(pixelCoord);
                        voxelColor = tFunc.getColor(val);
                        
                        x=(tFunc.getColor(getVoxel(is+1,js,ks)).a-tFunc.getColor(getVoxel(is-1,js,ks)).a)/2;
                        y=(tFunc.getColor(getVoxel(is,js+1,ks)).a-tFunc.getColor(getVoxel(is,js-1,ks)).a)/2;
                        z=(tFunc.getColor(getVoxel(is,js,ks+1)).a-tFunc.getColor(getVoxel(is,js,ks-1)).a)/2;
                       
                        s = tFunc.gradientWeighting(val);
            
                        deltaf = Math.sqrt( x * x + y * y + z * z);
                        deltafNormalized = deltaf / deltafMax;
                        
                        alpha = s * deltafNormalized;
                        
                        s_alpha = s_alpha + ( 1 - s_alpha ) * alpha;
                        s_red = s_red + voxelColor.r * ( 1 - s_alpha ) * alpha;
                        s_green = s_green + voxelColor.g * ( 1 - s_alpha ) * alpha;
                        s_blue = s_blue + voxelColor.b * ( 1 - s_alpha ) * alpha;
                      
                    }
                         // BufferedImage expects a pixel color packed as ARGB in an int
       
                        c_alpha = (s_alpha) <= 1.0 ? (int) Math.floor(s_alpha * 255) : 255;
                        c_red = (s_red) <= 1.0 ? (int) Math.floor(s_red * 255) : 255;
                        c_green = (s_green) <= 1.0 ? (int) Math.floor(s_green * 255) : 255;
                        c_blue = (s_blue <= 1.0) ? (int) Math.floor(s_blue * 255) : 255;
                        pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                       
                       image.setRGB(i, j, pixelColor);  
                }
            }
    }

    private void drawBoundingBox(GL2 gl) {
        gl.glPushAttrib(GL2.GL_CURRENT_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glColor4d(1.0, 1.0, 1.0, 1.0);
        gl.glLineWidth(1.5f);
        gl.glEnable(GL2.GL_LINE_SMOOTH);
        gl.glHint(GL2.GL_LINE_SMOOTH_HINT, GL2.GL_NICEST);
        gl.glEnable(GL2.GL_BLEND);
        gl.glBlendFunc(GL2.GL_SRC_ALPHA, GL2.GL_ONE_MINUS_SRC_ALPHA);

        gl.glBegin(GL2.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL2.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL2.GL_LINE_LOOP);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL2.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL2.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL2.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glDisable(GL2.GL_LINE_SMOOTH);
        gl.glDisable(GL2.GL_BLEND);
        gl.glEnable(GL2.GL_LIGHTING);
        gl.glPopAttrib();

    }

    @Override
    public void visualize(GL2 gl) {
        glCanvas=gl;

        if (volume == null) {
            return;
        }

        drawBoundingBox(gl);

        gl.glGetDoublev(GL2.GL_MODELVIEW_MATRIX, viewMatrix, 0);

        long startTime = System.currentTimeMillis();
        
          if (activeRenderer == 0) 
          {
            slicer(viewMatrix);
          } 
          else if (activeRenderer == 1) 
          {
            MIP(viewMatrix,gl);
          } 
          else if (activeRenderer == 2) 
          {
            backFront(viewMatrix,gl);
          }
          else if (activeRenderer == 3) 
          {
            frontBack(viewMatrix,gl);
          }
          else if (activeRenderer == 4) 
          {
            OpacityWeighting(viewMatrix,gl);
          }
        
        long endTime = System.currentTimeMillis();
        double runningTime = (endTime - startTime);
        panel.setSpeedLabel(Double.toString(runningTime));

        Texture texture = AWTTextureIO.newTexture(gl.getGLProfile(), image, false);

        gl.glPushAttrib(GL2.GL_LIGHTING_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glEnable(GL2.GL_BLEND);
        gl.glBlendFunc(GL2.GL_SRC_ALPHA, GL2.GL_ONE_MINUS_SRC_ALPHA);

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

    public void setScalar(double scalar) {
        this.scalar = scalar;  
    }
    
    private BufferedImage image;
    private double[] viewMatrix = new double[4 * 4];
}
