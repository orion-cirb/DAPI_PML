package DNA_Utils;


import Cellpose.CellposeSegmentImgPlusAdvanced;
import Cellpose.CellposeTaskSettings;
import DNA_PML_StardistOrion.StarDist2D;
import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.WaitForUserDialog;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.plugin.Duplicator;
import ij.plugin.ImageCalculator;
import ij.plugin.RGBStackMerge;
import ij.process.ImageProcessor;
import java.awt.Color;
import java.awt.Font;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import javax.swing.ImageIcon;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom2.BoundingBox;
import mcib3d.geom2.Object3DInt;
import mcib3d.geom2.Objects3DIntPopulation;
import mcib3d.geom2.measurements.MeasureCentroid;
import mcib3d.geom2.measurements.MeasureIntensity;
import mcib3d.geom2.measurements.MeasureVolume;
import mcib3d.image3d.ImageHandler;
import net.haesleinhuepf.clij.clearcl.ClearCLBuffer;
import net.haesleinhuepf.clij2.CLIJ2;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.scijava.util.ArrayUtils;


 /*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose dots_Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author phm
 */
public class DNA_Processing {
    
    // DNA filter size
    private double minDNA = 0.004;
    private double maxDNA = 1;
    // Nucleus filter
    private double minNuc = 50;
    private double maxNuc = Double.MAX_VALUE;
    // PML filter size
    private double minPML = 0.05;
    private double maxPML = 70;
    public String dna = "";

    
    private Calibration cal = new Calibration();
    private double pixVol = 0;
    public CLIJ2 clij2 = CLIJ2.getInstance();
    
    
    // Stardist
    public Object syncObject = new Object();
    public final double stardistPercentileBottom = 2;
    public final double stardistPercentileTop = 99.8;
    public final double stardistProbFociThresh = 0.2;
    public final double stardistOverlayFociThresh = 0.25;
    public File stardistModelsPath = new File(IJ.getDirectory("imagej")+File.separator+"models");
    public String stardistModel = "";
    public String stardistOutput = "Label Image"; 
    
    // Cellpose
    public int cellPoseDiameter = 50;
    public String cellPoseModel = "cyto2";
    private final String cellposeEnvDirPath = (IJ.isLinux()) ? "/opt/miniconda3/envs/cellpose" : System.getProperty("user.home")+"\\miniconda3\\envs\\CellPose";
    private double stitchThreshold = 0.25;
    public double distMax = 1.5;
    public final ImageIcon icon = new ImageIcon(this.getClass().getResource("/Orion_icon.png"));
    
    
     /**
     * check  installed modules
     * @return 
     */
    public boolean checkInstalledModules() {
        // check install
        ClassLoader loader = IJ.getClassLoader();
        try {
            loader.loadClass("net.haesleinhuepf.clij2.CLIJ2");
        } catch (ClassNotFoundException e) {
            IJ.log("CLIJ not installed, please install from update site");
            return false;
        }
        try {
            loader.loadClass("mcib3d.geom.Object3D");
        } catch (ClassNotFoundException e) {
            IJ.log("3D ImageJ Suite not installed, please install from update site");
            return false;
        }
        return true;
    }
    
    /**
     * Find image type
     */
    public String findImageType(File imagesFolder) {
        String ext = "";
        String[] files = imagesFolder.list();
        for (String name : files) {
            String fileExt = FilenameUtils.getExtension(name);
            switch (fileExt) {
                case "nd" :
                    return fileExt;
                case "czi" :
                  return fileExt;
                case "lif"  :
                    return fileExt;
                case "isc2" :
                   return fileExt;
                default :
                   ext = fileExt;
                   break; 
            }
        }
        return(ext);
    }
    
    /*
    Find starDist models in Fiji models folder
    */
    private String[] findStardistModels() throws IOException {
        FilenameFilter filter = (dir, name) -> name.startsWith("generic");
        File[] modelList = stardistModelsPath.listFiles(filter);
        for (File f : modelList) {
            if (f.exists())
               FileUtils.deleteDirectory(f); 
        }
        filter = (dir, name) -> name.endsWith(".zip");
        modelList = stardistModelsPath.listFiles(filter);
        String[] models = new String[modelList.length];
        for (int i = 0; i < modelList.length; i++) {
            models[i] = modelList[i].getName();
        }
        return(models);
    }
    
    
    /**
     * Find images in folder
     * @param imagesFolder
     * @param imageExt
     * @return 
     */
    public ArrayList<String> findImages(String imagesFolder, String imageExt) {
        File inDir = new File(imagesFolder);
        String[] files = inDir.list();
        if (files == null) {
            System.out.println("No Image found in "+imagesFolder);
            return null;
        }
        ArrayList<String> images = new ArrayList();
        for (String f : files) {
            // Find images with extension
            String fileExt = FilenameUtils.getExtension(f);
            if (fileExt.equals(imageExt))
                images.add(imagesFolder + File.separator + f);
        }
        return(images);
    }
       
     /**
     * Find channels name
     * @param imageName
     * @return 
     * @throws loci.common.services.DependencyException
     * @throws loci.common.services.ServiceException
     * @throws loci.formats.FormatException
     * @throws java.io.IOException
     */
    public String[] findChannels (String imageName, IMetadata meta, ImageProcessorReader reader, boolean bioformat) throws DependencyException, ServiceException, FormatException, IOException {
        int chs = reader.getSizeC();
        String[] channels = new String[chs];
        String imageExt =  FilenameUtils.getExtension(imageName);
        switch (imageExt) {
            case "nd" :
                for (int n = 0; n < chs; n++) 
                {
                    if (meta.getChannelID(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelName(0, n);
                    if (!bioformat) {
                        channels[n] = channels[n].replace("_", "-");
                        channels[n] = "w"+(n+1)+channels[n];
                    }
                }
                break;
            case "lif" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelName(0, n);
                break;
            case "czi" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelFluor(0, n);
                break;
            case "ics" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelExcitationWavelength(0, n).value().toString();
                break;    
            default :
                for (int n = 0; n < chs; n++)
                    channels[0] = Integer.toString(n);
        }
        return(channels);         
    }
    
    /**
     * Find image calibration
     * @param meta
     * @return 
     */
    public Calibration findImageCalib(IMetadata meta, ImageProcessorReader reader) {
        cal = new Calibration();  
        // read image calibration
        cal.pixelWidth = meta.getPixelsPhysicalSizeX(0).value().doubleValue();
        cal.pixelHeight = cal.pixelWidth;
        if (meta.getPixelsPhysicalSizeZ(0) != null)
            cal.pixelDepth = meta.getPixelsPhysicalSizeZ(0).value().doubleValue();
        else
            cal.pixelDepth = 1;
        cal.setUnit("microns");
        return(cal);
    }
    
    public Calibration getCalib()
    {
        System.out.println("x cal = " +cal.pixelWidth+", z cal=" + cal.pixelDepth);
        return cal;
    }
    
    
    /**
    * 
    * @param FileResults
    * @param resultsFileName
    * @param header
    * @return 
    */
    public BufferedWriter writeHeaders(String outFileResults, String header) throws IOException {
       FileWriter FileResults = new FileWriter(outFileResults, false);
       BufferedWriter outPutResults = new BufferedWriter(FileResults); 
       outPutResults.write(header);
       outPutResults.flush();
       return(outPutResults);
    }    
    
    /**
     * Dialog 
     */
    public int[] dialog(String[] channels) throws IOException {
        String[] models = findStardistModels();
        String[] chNames = {"Nucleus", "DNA/PML"};
        String[] dotsTypes = {"DNA", "PML"};
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.setInsets​(0, 80, 0);
        gd.addImage(icon);
        gd.addMessage("Channels selection", Font.getFont("Monospace"), Color.blue);
        for (String chName : chNames) {
            gd.addChoice(chName+" : ", channels, channels[0]);
        }
        gd.addMessage("DNA/PML parameters", Font.getFont("Monospace"), Color.blue);
        gd.addChoice("Dots type : ", dotsTypes, dna);
        gd.addMessage("StarDist model for foci detection", Font.getFont("Monospace"), Color.blue);
        if (models.length > 0) {
            gd.addChoice("StarDist model :",models, models[0]);
        }
        else {
            gd.addMessage("No StarDist model found in Fiji !!", Font.getFont("Monospace"), Color.red);
            gd.addFileField("StarDist model :", stardistModel);
        }
        gd.addNumericField("Min DNA/PML size (µm3) : ", minDNA, 3);
        gd.addNumericField("Max DNA/PML size (µm3) : ", maxDNA, 3);
        gd.addMessage("Image calibration", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("Calibration xy (µm)  :", cal.pixelWidth, 3);
        gd.addNumericField("Calibration z (µm)  :", cal.pixelDepth, 3);
        gd.showDialog();
        int[] chChoices = new int[chNames.length];
        for (int n = 0; n < chChoices.length; n++) {
            chChoices[n] = ArrayUtils.indexOf(channels, gd.getNextChoice());
        }
        dna = gd.getNextChoice();
       
        
        if (models.length > 0) {
            stardistModel = gd.getNextChoice();
        }
        else {
            stardistModel = gd.getNextString();
        }
        if (dna.equals("DNA")) {
            minDNA = gd.getNextNumber();
            maxDNA = gd.getNextNumber();
        }
        else {
            minPML = gd.getNextNumber();
            maxPML = gd.getNextNumber();
        }
        
        cal.pixelWidth = gd.getNextNumber();
        cal.pixelHeight = cal.pixelWidth;
        cal.pixelDepth = gd.getNextNumber();
        pixVol = cal.pixelWidth*cal.pixelHeight*cal.pixelDepth;
        if (gd.wasCanceled())
                chChoices = null;
        return(chChoices);
    }
    
    
   /* Median filter 
     * Using CLIJ2
     * @param ClearCLBuffer
     * @param sizeXY
     * @param sizeZ
     */ 
    public ClearCLBuffer median_filter(ClearCLBuffer  imgCL, double sizeXY, double sizeZ) {
        ClearCLBuffer imgCLMed = clij2.create(imgCL);
        clij2.mean3DBox(imgCL, imgCLMed, sizeXY, sizeXY, sizeZ);
        clij2.release(imgCL);
        return(imgCLMed);
    }
    
     /**
     * Remove object with size < min and size > max
     * @param pop
     * @param min
     * @param max
     */
    public void popFilterSize(Objects3DIntPopulation pop, double min, double max) {
        pop.getObjects3DInt().removeIf(p -> (new MeasureVolume(p).getVolumeUnit() < min) || (new MeasureVolume(p).getVolumeUnit() > max));
        pop.resetLabels();
    }
    
    
    /*
     * Remove objects present in only one z slice from population 
     */
    public Objects3DIntPopulation zFilterPop (Objects3DIntPopulation pop) {
        Objects3DIntPopulation popZ = new Objects3DIntPopulation();
        for (Object3DInt obj : pop.getObjects3DInt()) {
            int zmin = obj.getBoundingBox().zmin;
            int zmax = obj.getBoundingBox().zmax;
            if (zmax != zmin)
                popZ.addObject(obj);
        }
        return popZ;
    }
   
    /**
     * Look for all 3D cells in a Z-stack: 
     * - apply CellPose in 2D slice by slice 
     * - let CellPose reconstruct cells in 3D using the stitch threshold parameters
     * @param img
     * @return 
     * @throws java.io.IOException
     */
    public Objects3DIntPopulation cellposeDetection(ImagePlus img) throws IOException{
        boolean reSize = false;
        double resizeFactor = 0.25f;
        ImagePlus imgResized;
        if (img.getWidth() > 1024) {
            reSize = true;
            imgResized = img.resize((int)(img.getWidth()*resizeFactor), (int)(img.getHeight()*resizeFactor), 1, "none");
        }
        else
            imgResized = new Duplicator().run(img);

        // Define CellPose settings
        CellposeTaskSettings settings = new CellposeTaskSettings(cellPoseModel, 1, cellPoseDiameter, cellposeEnvDirPath);
        settings.setStitchThreshold(stitchThreshold);
        settings.useGpu(true);
       
        // Run CellPose
        CellposeSegmentImgPlusAdvanced cellpose = new CellposeSegmentImgPlusAdvanced(settings, imgResized);
        ImagePlus imgOut = cellpose.run();
        ImagePlus imgNuc = (reSize) ? imgOut.resize(img.getWidth(), img.getHeight(), "none") : imgOut;
        imgNuc.setCalibration(cal);
        // Get cells as a population of objects
        Objects3DIntPopulation pop = new Objects3DIntPopulation(ImageHandler.wrap(imgNuc));
        System.out.println(pop.getNbObjects() + " CellPose detections");
        flush_close(imgNuc);
        flush_close(imgOut);
        // Remove cell with only one Z
        zFilterPop(pop);
        // Filter cells by size
        popFilterSize(pop, minNuc, maxNuc);
        return(pop);
    } 
    
     /**
     * Apply StarDist 2D slice by slice
     * Label detections in 3D
     * @return objects population
     */
    public Objects3DIntPopulation stardistPop(ImagePlus img, String type) throws IOException{
        // Resize image to be in a StarDist-friendly scale
       ImagePlus imgIn = null;
       float factor = 0.5f;
       boolean resize = false;
       double minVol = 0, maxVol = 0;
       String model = ""; 
       if (img.getWidth() > 1024) {
           imgIn = img.resize((int)(img.getWidth()*factor), (int)(img.getHeight()*factor), 1, "none");
           resize = true;
       } 
       else
           imgIn = new Duplicator().run(img);
      
        switch (dna) {
            case "DNA" :
                minVol = minDNA;
                maxVol = maxDNA;
             break;
            case "PML" :
                minVol = minPML;
                maxVol = maxPML;
             break;
        }
        model = stardistModel;
        // StarDist
        File starDistModelFile = new File(stardistModelsPath+File.separator+model);
        StarDist2D star = new StarDist2D(syncObject, starDistModelFile);
        star.loadInput(imgIn);
        star.setParams(stardistPercentileBottom, stardistPercentileTop, stardistProbFociThresh, stardistOverlayFociThresh, stardistOutput);
        star.run();
        flush_close(imgIn);

        // Label detections in 3D
        ImagePlus imgOut = (resize) ? star.getLabelImagePlus().resize(img.getWidth(), img.getHeight(), 1, "none") : star.getLabelImagePlus();       
        ImagePlus imgLabels = star.associateLabels(imgOut);
        imgLabels.setCalibration(cal); 
        flush_close(imgOut);
        Objects3DIntPopulation pop = new Objects3DIntPopulation(ImageHandler.wrap(imgLabels));
        System.out.println("Filter size min="+minVol+" max="+maxVol+" of "+pop.getNbObjects()+" foci");
        popFilterSize(pop, minVol, maxVol);
       flush_close(imgLabels);
       return(pop);
   }
    
     /**
     * Find foci inside nucleus
     * @param fociPop
     * @param nucObj
     */
    public Objects3DIntPopulation fociInNucleus(Objects3DIntPopulation fociPop, Object3DInt nucObj) {
        Objects3DIntPopulation fociNucPop = new Objects3DIntPopulation();
        for (Object3DInt foci : fociPop.getObjects3DInt()) {
            MeasureCentroid fociCenter = new MeasureCentroid(foci);
            if (nucObj.contains(fociCenter.getCentroidRoundedAsVoxelInt())){
               fociNucPop.addObject(foci); 
            }
        }
       return(fociNucPop);    }
    
    
    
     /**
     * Label object
     * @param pop
     * @param img 
     * @param fontSize 
     */
    public void labelObject(Objects3DIntPopulation pop, ImagePlus img, int fontSize) {
        if (IJ.isMacOSX())
            fontSize *= 3;
        for (Object3DInt obj : pop.getObjects3DInt()) {
            BoundingBox bb = obj.getBoundingBox();
            int z = bb.zmax - bb.zmin;
            int x = bb.xmin;
            int y = bb.ymin;
            img.setSlice(z);
            ImageProcessor ip = img.getProcessor();
            ip.setFont(new Font("SansSerif", Font.PLAIN, fontSize));
            ip.setColor(255);
            ip.drawString(Integer.toString((int)obj.getLabel()), x, y);
            img.updateAndDraw();
        }
    }
    
    
// Flush and close images
    public void flush_close(ImagePlus img) {
        img.flush();
        img.close();
    } 
    
    
    /**
     * find diffus intensity in nucleus
     * @param img
     * @param dots
     * @param nucleus
     * @return 
     */
    public double getDiffus(ImagePlus img, Objects3DIntPopulation fociPop, Object3DInt nucleus) {
        ImageHandler imgDiffus = ImageHandler.wrap(img.duplicate());
        for (Object3DInt obj : fociPop.getObjects3DInt())
            obj.drawObject(imgDiffus, 0);
        double diffusInt = new MeasureIntensity(nucleus, imgDiffus).getValueMeasurement(MeasureIntensity.INTENSITY_SUM);
        return(diffusInt);
    }
    
     /**
     * Find sum volume of objects  
     * @param dotsPop
     * @return vol
     */
    
    public double findPopVolume (Objects3DIntPopulation dotsPop) {
        IJ.showStatus("Findind object's volume");
        List<Double[]> results = dotsPop.getMeasurementsList(new MeasureVolume().getNamesMeasurement());
        double sum = results.stream().map(arr -> arr[1]).reduce(0.0, Double::sum);
        return(sum);
    }
    
     /**
     * Find sum intensity of objects  
     * @param dotsPop
     * @param img
     * @return intensity
     */
    
    public double findPopIntensity (Objects3DIntPopulation dotsPop, ImagePlus img) {
        IJ.showStatus("Findind object's intensity");
        ImageHandler imh = ImageHandler.wrap(img);
        double sumInt = 0;
        for(Object3DInt obj : dotsPop.getObjects3DInt()) {
            MeasureIntensity intMes = new MeasureIntensity(obj, imh);
            sumInt +=  intMes.getValueMeasurement(MeasureIntensity.INTENSITY_SUM);
        }
        return(sumInt);
    }
    
    
    /**
     * Find foci associated or not to all nucleus
     * 
     * @param fociPop
     * @param nucPop
     * @param img
     * @return 
     */
   public Objects3DIntPopulation findFociInNuclei(Objects3DIntPopulation fociPop, Objects3DIntPopulation nucPop, ImagePlus img) {
        ImageHandler imh = ImageHandler.wrap(img).createSameDimensions();
        fociPop.drawInImage(imh);
        ImageHandler imhCopy = imh.duplicate();
        for (Object3DInt obj : nucPop.getObjects3DInt()) {
            obj.drawObject(imh, 0);
        }

        ImageCalculator imgCal = new ImageCalculator();
        ImagePlus imgSub = imgCal.run("subtract stack create", imhCopy.getImagePlus(), imh.getImagePlus());
        ImageHandler imhSub = ImageHandler.wrap(imgSub);
        Objects3DIntPopulation pop = new Objects3DIntPopulation(imhSub);
        imh.closeImagePlus();
        imhCopy.closeImagePlus();
        return(pop);  
   }
   

    // Save objects image
    public void saveObjects (ImagePlus img, Objects3DIntPopulation nucPop, Objects3DIntPopulation fociPop, String imageName) {
        ImageHandler imgNucObjects = ImageHandler.wrap(img).createSameDimensions();
        ImageHandler imgFociObjects = imgNucObjects.duplicate();
        for (Object3DInt obj : nucPop.getObjects3DInt())
            obj.drawObject(imgNucObjects, 255);
        labelObject(nucPop, imgNucObjects.getImagePlus(), 30);
        for (Object3DInt obj : fociPop.getObjects3DInt())
            obj.drawObject(imgFociObjects, 255);
        ImagePlus[] imgColors = {imgFociObjects.getImagePlus(), null, imgNucObjects.getImagePlus()};
        ImagePlus imgObjects = new RGBStackMerge().mergeHyperstacks(imgColors, false);
        FileSaver ImgObjectsFile = new FileSaver(imgObjects);
        ImgObjectsFile.saveAsTiff(imageName);
        imgNucObjects.closeImagePlus();
        imgFociObjects.closeImagePlus();
        flush_close(imgObjects);
    }
}
