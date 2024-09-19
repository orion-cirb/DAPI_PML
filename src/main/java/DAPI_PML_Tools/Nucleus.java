package DAPI_PML_Tools;

import ij.ImagePlus;
import java.util.HashMap;
import mcib3d.geom2.Object3DComputation;
import mcib3d.geom2.Object3DInt;
import mcib3d.geom2.Objects3DIntPopulation;
import mcib3d.geom2.measurements.MeasureIntensity;
import mcib3d.geom2.measurements.MeasureVolume;
import mcib3d.image3d.ImageHandler;

/**
 * @author Héloïse Monnet @ ORION-CIRB
 */
public class Nucleus {
    
    public Object3DInt nucleus;
    public Objects3DIntPopulation pmlFoci;
    
    public Nucleus(Object3DInt nucleus, Objects3DIntPopulation pmlFoci) {
        this.nucleus = nucleus;
        this.pmlFoci = pmlFoci;
    }
    
    public HashMap<String, Double> computeParams(ImagePlus imgPml, double bgPml, double pixArea) {
        HashMap<String, Double> params = new HashMap<>();
        
        // Nucleus
        params.put("nucId", (double) nucleus.getLabel());
        double nucArea = new MeasureVolume(nucleus).getVolumeUnit();
        params.put("nucArea", nucArea);
        
        // PML foci
        double[] pmlFociParams = computeFociParams(imgPml, bgPml);
        params.put("pmlFociNb", (double) pmlFoci.getNbObjects());
        params.put("pmlFociArea", pmlFociParams[0]);
        params.put("pmlFociMeanInt", pmlFociParams[1]);
        params.put("pmlFociTotInt", pmlFociParams[2]);
        
        // PML diffuse
        double[] pmlDiffuseParams = computeFociDiffuseParams(imgPml, bgPml, nucArea, pixArea);
        params.put("pmlDiffuseArea", pmlDiffuseParams[0]);
        params.put("pmlDiffuseMeanInt", pmlDiffuseParams[1]);
        params.put("pmlDiffuseTotInt", pmlDiffuseParams[2]);
                   
        return(params);
    }
    
    
    /**
     * Compute parameters of PML foci in nucleus:
     * - area
     * - background-corrected mean intensity 
     * - background-corrected raw integrated density
     */
    public double[] computeFociParams(ImagePlus imgPml, double bgPml) {
        ImageHandler imhPml = ImageHandler.wrap(imgPml);
        
        ImageHandler imhFoci = imhPml.createSameDimensions();
        for (Object3DInt foci: pmlFoci.getObjects3DInt())
            foci.drawObject(imhFoci, 255);
        Object3DInt objFoci = new Object3DInt(imhFoci);
        
        double area = new MeasureVolume(objFoci).getVolumeUnit();
        double meanInt = new MeasureIntensity(objFoci, imhPml).getValueMeasurement(MeasureIntensity.INTENSITY_AVG) - bgPml;
        double totInt = new MeasureIntensity(objFoci, imhPml).getValueMeasurement(MeasureIntensity.INTENSITY_SUM) - bgPml * new MeasureVolume(objFoci).getVolumePix();
        double[] params = {area, meanInt, totInt};
        
        imhFoci.closeImagePlus();
        return(params);
    }
    
    
    /**
     * Compute parameters of PML diffuse signal in nucleus:
     * - area
     * - background-corrected mean intensity 
     * - background-corrected raw integrated density
     */
    public double[] computeFociDiffuseParams(ImagePlus imgPml, double bgPml, double nucArea, double pixArea) {
        ImageHandler imhFoci = ImageHandler.wrap(imgPml.duplicate());
        
        double dilFociArea = 0;
        for (Object3DInt foci: pmlFoci.getObjects3DInt()) {
            Object3DInt dilatedObj = new Object3DComputation(foci).getObjectDilated(2f, 2f, 0);
            dilatedObj.drawObject(imhFoci, 0);
            dilFociArea += new MeasureVolume(dilatedObj).getVolumeUnit();
        }
        
        double area = nucArea - dilFociArea;
        double areaInPix = area / pixArea;
        double totInt = new MeasureIntensity(nucleus, imhFoci).getValueMeasurement(MeasureIntensity.INTENSITY_SUM) - bgPml * areaInPix;
        double meanInt = totInt / areaInPix;
        double[] params = {area, meanInt, totInt};
        imhFoci.closeImagePlus();
        return(params);
    }
    
}
