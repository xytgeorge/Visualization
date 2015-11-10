/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volvis;

import java.awt.Color;
import java.util.ArrayList;
import util.TFChangeListener;
import volume.Volume;

/**
 *
 * @author michel
 */
public class TransferFunction {

    private ArrayList<TFChangeListener> listeners = new ArrayList<TFChangeListener>();
    private boolean opacityWeighting = false;
    Volume vol;
    double rVal  = 1.0;
    // Construct a default grey-scale transfer function over the scalar range min - max.
    // The opacity increases linearly from 0.0 to 1.0
    public TransferFunction(short min, short max, Volume volume) {
        vol = volume;
        sMin = min;
        sMax = max;
        sRange = sMax - sMin;
        controlPoints = new ArrayList<ControlPoint>();

        controlPoints.add(new ControlPoint(min, new TFColor(0.0, 0.0, 0.0, 0.0)));
        controlPoints.add(new ControlPoint(max, new TFColor(1.0, 1.0, 1.0, 1.0)));

        LUTsize = sRange;
        LUT = new TFColor[LUTsize];

        buildLUT();

    }

    public int getMinimum() {
        return sMin;
    }

    public int getMaximum() {
        return sMax;
    }

    public ArrayList<ControlPoint> getControlPoints() {
        return controlPoints;
    }

    public TFColor getColor(int value) {
        return LUT[computeLUTindex(value)];
    }
    
    public void setColor(int value, TFColor col) {
        LUT[computeLUTindex(value)] = col;
    }

    public int addControlPoint(int value, double r, double g, double b, double a) {
        if (value < sMin || value > sMax) {
            return -1;
        }
        a = Math.floor(a * 100) / 100.0;

        ControlPoint cp = new ControlPoint(value, new TFColor(r, g, b, a));
        int idx = 0;
        while (idx < controlPoints.size() && controlPoints.get(idx).compareTo(cp) < 0) {
            idx++;
        }


        if (controlPoints.get(idx).compareTo(cp) == 0) {
            controlPoints.set(idx, cp);
        } else {
            controlPoints.add(idx, cp);
        }

        buildLUT();
        return idx;
    }

    public void removeControlPoint(int idx) {
        controlPoints.remove(idx);
        buildLUT();
    }

    public void updateControlPointScalar(int index, int s) {
        controlPoints.get(index).value = s;
        buildLUT();
    }

    public void updateControlPointAlpha(int index, double alpha) {
        alpha = Math.floor(alpha * 100) / 100.0;
        controlPoints.get(index).color.a = alpha;
        buildLUT();
    }

    public void updateControlPointColor(int idx, Color c) {
        ControlPoint cp = controlPoints.get(idx);
        cp.color.r = c.getRed() / 255.0;
        cp.color.g = c.getGreen() / 255.0;
        cp.color.b = c.getBlue() / 255.0;
        buildLUT();
    }

    // Add a changelistener, which will be notified is the transfer function changes
    public void addTFChangeListener(TFChangeListener l) {
        if (!listeners.contains(l)) {
            listeners.add(l);
        }
    }

    // notify the change listeners
    public void changed() {
        for (int i = 0; i < listeners.size(); i++) {
            listeners.get(i).changed();
        }
    }

    private int computeLUTindex(int value) {
        int idx = ((LUTsize - 1) * (value - sMin)) / sRange;
        return idx;
    }

    private void buildLUT() {

        for (int i = 1; i < controlPoints.size(); i++) {
            ControlPoint prev = controlPoints.get(i - 1);
            ControlPoint next = controlPoints.get(i);
            //System.out.println(prev.value + " " + prev.color + " -- " + next.value + " " + next.color);
            double range = next.value - prev.value;
            for (int k = prev.value; k <= next.value; k++) {
                double frac = (k - prev.value) / range;
                TFColor newcolor = new TFColor();
                newcolor.r = prev.color.r + frac * (next.color.r - prev.color.r);
                newcolor.g = prev.color.g + frac * (next.color.g - prev.color.g);
                newcolor.b = prev.color.b + frac * (next.color.b - prev.color.b);
                newcolor.a = prev.color.a + frac * (next.color.a - prev.color.a);
                LUT[computeLUTindex(k)] = newcolor;
            }

        }
          if(opacityWeighting){
            gradientWeighting();
        }

    }
    
     public void enableOpacityWeighting(){
            opacityWeighting = true;
            buildLUT();
            changed();
        }
    public void disableOpacityWeighting(){
            opacityWeighting = false;
            buildLUT();
            changed();
        }
    
    
    public void gradientWeighting(){
        for(int i = 0; i < LUTsize; i++){
            double alphaFinal = -1;
            double alpha = 0;
            for (int j = 1; j < controlPoints.size(); j++) {
                double Fv = controlPoints.get(j).value; //the value where we want it to be exactly alphaV
                double alphaV = controlPoints.get(j).color.a;
                ///
                double rVoxelThickness = Fv/rVal;
                
                int[] histogram = vol.getHistogram();
                int[] histogramDf = new int[histogram.length];
                int[] histogramDDf = new int[histogramDf.length];
                for(int q = 1; q < histogram.length -1; q++){
                    histogramDf[q-1] = (histogram[q+1]-histogram[q])/1;
                }
                for(int q = 1; q < histogramDf.length -1; q++){
                    histogramDDf[q-1] = (histogramDf[q+1]-histogramDf[q])/1;
                }
                for(int q = 1; q < histogramDf.length-1; q++){
                    if ((histogramDf[q+1] < histogramDf[q]) && (histogramDf[q-1] < histogramDf[q])){
                        //i is max
                        if((Math.signum(histogramDDf[q+1]) != Math.signum(histogramDDf[q-1])) && histogramDDf[q] == 0){
                            //determine opacity as per paper eqn.3

                            if((Math.abs(histogramDf[q]) == 0) && (histogram[q] == Fv)){ alpha = 1;}

                            else if((Math.abs(histogramDf[q]) > 0) &&
                                    ((histogram[q] - rVoxelThickness * Math.abs(histogramDf[q])) <= Fv) &&
                                    (Fv <= histogram[q] + rVoxelThickness * Math.abs(histogramDf[q]))){

                                alpha = 1 - (1/rVoxelThickness) * Math.abs((Fv - histogram[q])/(Math.abs(histogramDf[q])));
                            }
                            else{alpha = 0;}
                        }
                    }
                }
                alpha = alpha * alphaV;
            if (alphaFinal < 0.0){ alphaFinal = (1 - alpha);}
            else{alphaFinal = alphaFinal * (1 - alpha);}
                
                ///
                
            }
            alphaFinal = 1 - alphaFinal;
            //set alpha for this point based on all the CP's
            TFColor tCol = getColor(i);
            tCol.a = alphaFinal;
            setColor(i, tCol);
        }
    }
    
    public void setRValue(Integer r){
        this.rVal = r;
        //gradientWeighting();
        changed();
    }
    

    public class ControlPoint implements Comparable<ControlPoint> {

        public int value;
        public TFColor color;

        public ControlPoint(int v, TFColor c) {
            value = v;
            color = c;
        }

        @Override
        public int compareTo(ControlPoint t) {
            return (value < t.value ? -1 : (value == t.value ? 0 : 1));
        }

        @Override
        public String toString() {
            return new String("(" + value + ") -> " + color.toString());
        }
    }
    private short sMin, sMax;
    private int sRange;
    private TFColor[] LUT;
    private int LUTsize = 4095;
    private ArrayList<ControlPoint> controlPoints;
}
