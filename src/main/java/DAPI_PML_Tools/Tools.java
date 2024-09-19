package DAPI_PML_Tools;

import DAPI_PML_Tools.Cellpose.CellposeSegmentImgPlusAdvanced;
import DAPI_PML_Tools.Cellpose.CellposeTaskSettings;
import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.Roi;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.RGBStackMerge;
import ij.plugin.filter.Analyzer;
import ij.process.AutoThresholder;
import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import javax.swing.ImageIcon;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom2.Object3DInt;
import mcib3d.geom2.Objects3DIntPopulation;
import mcib3d.geom2.Objects3DIntPopulationComputation;
import mcib3d.geom2.measurements.MeasureCentroid;
import mcib3d.geom2.measurements.MeasureIntensity;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageLabeller;
import net.haesleinhuepf.clij.clearcl.ClearCLBuffer;
import net.haesleinhuepf.clij2.CLIJ2;
import org.apache.commons.io.FilenameUtils;


/**
 * @author Héloïse Monnet @ ORION-CIRB
 */
public class Tools {
    
    public final ImageIcon icon = new ImageIcon(this.getClass().getResource("/Orion_icon.png"));
    private final String helpUrl = "https://github.com/orion-cirb/DAPI_PML";
    
    private CLIJ2 clij2 = CLIJ2.getInstance();
    
    private String[] chDialog = new String[]{"DAPI nuclei", "PML foci"};
    private Calibration cal;
    public double pixArea;
   
    // Nuclei detection with Cellpose
    private String cellposeEnvDir = IJ.isWindows()? System.getProperty("user.home")+File.separator+"miniconda3"+File.separator+"envs"+File.separator+"CellPose" : "/opt/miniconda3/envs/cellpose";
    private final String cellposeModelPath = IJ.isWindows()? System.getProperty("user.home")+"\\.cellpose\\models\\" : "";
    private String cellposeModel = "cyto2";
    private int cellposeDiam = 100;
    private double minAreaNuc = 50; // µm2
    private double maxAreaNuc = 550; // µm2
    
    // PML foci detection with DoG + Thresholding
    private final double dogSigma1 = 1;
    private final double dogSigma2 = 3;
    private String thMethod = "Triangle";
    private double minAreaFoci = 0.05; // µm2
    private double maxAreaFoci = 3; // µm2
    

    /**
     * Display a message in the ImageJ console and status bar
     */
    public void print(String log) {
        System.out.println(log);
        IJ.showStatus(log);
    }
    
    
    /**
     * Flush and close an image
     */
    public void closeImage(ImagePlus img) {
        img.flush();
        img.close();
    }
    
    
    /**
     * Check that needed modules are installed
     */
    public boolean checkInstalledModules() {
        ClassLoader loader = IJ.getClassLoader();
        try {
            loader.loadClass("net.haesleinhuepf.clij2.CLIJ2");
        } catch (ClassNotFoundException e) {
            IJ.log("CLIJ not installed, please install from update site");
            return false;
        }
        try {
            loader.loadClass("mcib3d.geom2.Object3DInt");
        } catch (ClassNotFoundException e) {
            IJ.log("3D ImageJ Suite not installed, please install from update site");
            return false;
        }
        return true;
    }
    
    
    /**
     * Get extension of the first image found in the folder
     */
    public String findImageType(File imagesFolder) {
        String ext = "";
        String[] files = imagesFolder.list();
        for (String name : files) {
            String fileExt = FilenameUtils.getExtension(name);
            switch (fileExt) {
                case "nd" :
                   ext = fileExt;
                   break;
                case "nd2" :
                   ext = fileExt;
                   break;
                case "lif"  :
                    ext = fileExt;
                    break;
                case "czi" :
                   ext = fileExt;
                   break;
                case "ics" :
                    ext = fileExt;
                    break;
                case "ics2" :
                    ext = fileExt;
                    break;
                case "lsm" :
                    ext = fileExt;
                    break;
                case "tif" :
                    ext = fileExt;
                    break;
                case "tiff" :
                    ext = fileExt;
                    break;
            }
        }
        return(ext);
    }
    
    
    /**
     * Get images with given extension in folder
     */
    public ArrayList<String> findImages(String imagesFolder, String imageExt) {
        File inDir = new File(imagesFolder);
        String[] files = inDir.list();
        ArrayList<String> images = new ArrayList();
        for (String f : files) {
            String fileExt = FilenameUtils.getExtension(f);
            if (fileExt.equals(imageExt) && !f.startsWith("."))
                images.add(imagesFolder + f);
        }
        Collections.sort(images);
        return(images);
    }
       
    
    /**
     * Get image calibration
     */
    public void findImageCalib(IMetadata meta) {
        cal = new Calibration();
        cal.pixelWidth = meta.getPixelsPhysicalSizeX(0).value().doubleValue();
        cal.pixelHeight = cal.pixelWidth;
        if (meta.getPixelsPhysicalSizeZ(0) != null)
            cal.pixelDepth = meta.getPixelsPhysicalSizeZ(0).value().doubleValue();
        else
            cal.pixelDepth = 1;
        cal.setUnit("microns");
        System.out.println("XY calibration = " + cal.pixelWidth + ", Z calibration = " + cal.pixelDepth);
    }
    
    
    /**
     * Get channels name and add None at the end of channels list
     * @throws loci.common.services.DependencyException
     * @throws loci.common.services.ServiceException
     * @throws loci.formats.FormatException
     * @throws java.io.IOException
     */
    public String[] findChannels(String imageName, IMetadata meta, ImageProcessorReader reader) throws DependencyException, ServiceException, FormatException, IOException {
        int chs = reader.getSizeC();
        String[] channels = new String[chs];
        String imageExt =  FilenameUtils.getExtension(imageName);
        switch (imageExt) {
            case "nd" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = (meta.getChannelName(0, n).toString().equals("")) ? Integer.toString(n) : meta.getChannelName(0, n).toString();
                break;
            case "nd2" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = (meta.getChannelName(0, n).toString().equals("")) ? Integer.toString(n) : meta.getChannelName(0, n).toString();
                break;
            case "lif" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null || meta.getChannelName(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelName(0, n).toString();
                break;
            case "czi" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = (meta.getChannelFluor(0, n).toString().equals("")) ? Integer.toString(n) : meta.getChannelFluor(0, n).toString();
                break;
            case "ics" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = meta.getChannelEmissionWavelength(0, n).value().toString();
                break;    
            case "ics2" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = meta.getChannelEmissionWavelength(0, n).value().toString();
                break; 
            default :
                for (int n = 0; n < chs; n++)
                    channels[n] = Integer.toString(n);
        }
        return(channels);     
    }
    
    
    /**
     * Generate dialog box
     */
    public String[] dialog(String[] chMeta) {
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.setInsets​(0, 70, 0);
        gd.addImage(icon);
        
        gd.addMessage("Channels", Font.getFont("Monospace"), Color.blue);
        for (int n = 0; n < chDialog.length; n++)
            gd.addChoice(chDialog[n]+": ", chMeta, chMeta[n]);
        
        gd.addMessage("Nuclei detection", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("Min nucleus area (µm2): ", minAreaNuc, 2);
        gd.addNumericField("Max nucleus area (µm2): ", maxAreaNuc, 2);
        
        gd.addMessage("Foci detection", Font.getFont("Monospace"), Color.blue);
        String[] thMethods = AutoThresholder.getMethods();
        gd.addChoice("Thresholding method: ",thMethods, thMethod);
        gd.addNumericField("Min foci area (µm2): ", minAreaFoci, 3);
        gd.addNumericField("Max foci area (µm2): ", maxAreaFoci, 3);
        
        gd.addMessage("Image calibration", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("XY calibration (µm): ", cal.pixelHeight, 4);
        gd.addHelp(helpUrl);
        gd.showDialog();
        
        String[] chOrder = new String[chDialog.length];
        for (int n = 0; n < chOrder.length; n++)
            chOrder[n] = gd.getNextChoice();
        
        minAreaNuc = gd.getNextNumber();
        maxAreaNuc = gd.getNextNumber();
        
        thMethod = gd.getNextChoice();
        minAreaFoci = gd.getNextNumber();
        maxAreaFoci = gd.getNextNumber();

        cal.pixelHeight = cal.pixelWidth = gd.getNextNumber();
        cal.pixelDepth = 1;
        pixArea = cal.pixelHeight*cal.pixelWidth;
        
        if (gd.wasCanceled())
            chOrder = null;
        return(chOrder);
    }
    
    
    /**
     * Detect objects in 2D using Cellpose
     */
    public ImagePlus cellposeDetection(ImagePlus imgIn) {
        // Define CellPose settings
        CellposeTaskSettings settings = new CellposeTaskSettings(cellposeModelPath+cellposeModel, 1, cellposeDiam, cellposeEnvDir);
        settings.useGpu(true);

        // Run Cellpose
        ImagePlus img = imgIn.duplicate();
        CellposeSegmentImgPlusAdvanced cellpose = new CellposeSegmentImgPlusAdvanced(settings, img);
        ImagePlus imgOut = cellpose.run();
        
        closeImage(img);
        return(imgOut);
    }
    
    
    public double computeBackground(ImagePlus maskPop, ImagePlus img) {
        ImagePlus maskBg = maskPop.duplicate();
        IJ.setRawThreshold(maskBg, 1, 65535);
        IJ.run(maskBg, "Convert to Mask", "");
        IJ.run(maskBg, "Invert", "");
        IJ.run(maskBg, "Create Selection", "");
        Roi roiBg = maskBg.getRoi();
        
        ImagePlus imgBg = img.duplicate();
        imgBg.setRoi(roiBg);
        
        ResultsTable rt = new ResultsTable();
        Analyzer analyzer = new Analyzer(imgBg, Analyzer.MEDIAN, rt);
        analyzer.measure();
        double bg = rt.getValue("Median", 0);
        
        closeImage(maskBg);
        closeImage(imgBg);
        return(bg);
    }
    
    
    /**
     * Get population of objects from labelled mask and filter them out by area
     */
    public Objects3DIntPopulation filterPop(ImagePlus mask) {
        Objects3DIntPopulation pop = new Objects3DIntPopulation(ImageHandler.wrap(mask));
        pop = new Objects3DIntPopulationComputation(pop).getExcludeBorders(ImageHandler.wrap(mask), false);
        System.out.println("Nb nuclei detected: "+pop.getNbObjects());
        pop = new Objects3DIntPopulationComputation(pop).getFilterSize(minAreaNuc/pixArea, maxAreaNuc/pixArea);
        System.out.println("Nb nuclei remaining after size filtering: "+ pop.getNbObjects());
        pop.resetLabels();
        return(pop);
    }
    
    
    /**
     * Detect dots with DoG filtering + automatic thresholding + fill holes
     */
    public Objects3DIntPopulation fociDetection(ImagePlus imgIn) {
        ImagePlus imgDOG = DOG(imgIn, dogSigma1, dogSigma2);
        ImagePlus imgBin = threshold(imgDOG, thMethod);
        ImagePlus imgFill = fillHoles(imgBin);
        imgFill.setCalibration(cal);
        
        Objects3DIntPopulation pop = getPopFromImage(imgFill);
        System.out.println("Nb PML foci detected: "+pop.getNbObjects());
        pop = new Objects3DIntPopulationComputation(pop).getFilterSize(minAreaFoci/pixArea, maxAreaFoci/pixArea);
        System.out.println("Nb PML foci remaining after size filtering: "+ pop.getNbObjects());
        pop.resetLabels();
        
        closeImage(imgDOG);
        closeImage(imgBin);
        closeImage(imgFill);
        return(pop);
    }

    
    /**
     * Difference of Gaussians filtering using CLIJ2
     */ 
    public ImagePlus DOG(ImagePlus img, double size1, double size2) {
        ClearCLBuffer imgCL = clij2.push(img);
        ClearCLBuffer imgCLDOG = clij2.create(imgCL);
        clij2.differenceOfGaussian2D(imgCL, imgCLDOG, size1, size1, size2, size2);
        ImagePlus imgDOG = clij2.pull(imgCLDOG);
        clij2.release(imgCL);
        clij2.release(imgCLDOG);
        return(imgDOG);
    }
    
    
    /**
     * Automatic thresholding using CLIJ2
     */
    public ImagePlus threshold(ImagePlus img, String thMed) {
        ClearCLBuffer imgCL = clij2.push(img);
        ClearCLBuffer imgCLBin = clij2.create(imgCL);
        clij2.automaticThreshold(imgCL, imgCLBin, thMed);
        ImagePlus imgBin = clij2.pull(imgCLBin);
        clij2.release(imgCL);
        clij2.release(imgCLBin);
        return(imgBin);
    }
    
    
    /**
     * Fill holes using CLIJ2
     */
    public ImagePlus fillHoles(ImagePlus img) {
        ClearCLBuffer imgCL = clij2.push(img);
        ClearCLBuffer imgCLBin = clij2.create(imgCL);
        clij2.binaryFillHoles(imgCL, imgCLBin);
        ImagePlus imgBin = clij2.pull(imgCLBin);
        clij2.release(imgCL);
        clij2.release(imgCLBin);
        return(imgBin);
    }
    
        
    /**
     * Get population of objects from binary mask
     */
    public Objects3DIntPopulation getPopFromImage(ImagePlus img) {
        ImageInt labels = new ImageLabeller().getLabels(ImageHandler.wrap(img));
        Objects3DIntPopulation pop = new Objects3DIntPopulation(labels);
        labels.closeImagePlus();
        return(pop);
    }
    
    
    /**
     * Return list of Nucleus with their respective population of PML foci
     */
    public ArrayList<Nucleus> colocalizeNucFoci(Objects3DIntPopulation nucPop, Objects3DIntPopulation pmlPop) {
        ArrayList<Nucleus> nuclei = new ArrayList();
        
        nucPop.getObjects3DInt().stream().forEach(nuc -> {
            Objects3DIntPopulation pmlInNucPop = new Objects3DIntPopulation();
            pmlPop.getObjects3DInt().stream()
                .filter(pml -> nuc.contains(new MeasureCentroid(pml).getCentroidRoundedAsVoxelInt()))
                .forEach(pml -> {pmlInNucPop.addObject(pml);});
                
            nuclei.add(new Nucleus(nuc, pmlInNucPop));
        });
        
        return(nuclei);
    }
    
    
    /**
     * Draw results
     */
    public void drawResults(ArrayList<Nucleus> nuclei, ImagePlus imgPml, ImagePlus imgDapi, String outDir, String imgName) {
        ImageHandler imhNuc = ImageHandler.wrap(imgPml).createSameDimensions();
        ImageHandler imhPml = imhNuc.createSameDimensions();
        
        for (Nucleus nucleus: nuclei) {
            nucleus.nucleus.drawObject(imhNuc);
            for (Object3DInt pml: nucleus.pmlFoci.getObjects3DInt())
                pml.drawObject(imhPml, 255);
        }
        
        IJ.run(imhNuc.getImagePlus(), "glasbey on dark", "");
        IJ.run(imhNuc.getImagePlus(), "Enhance Contrast", "saturated=0.35");
        IJ.resetMinAndMax(imgDapi);
        IJ.run(imhPml.getImagePlus(), "Green", "");
        IJ.run(imhPml.getImagePlus(), "Enhance Contrast", "saturated=0.35");
        IJ.run(imgPml, "Enhance Contrast", "saturated=0.35");

        ImagePlus[] imgColors = {imhNuc.getImagePlus(), null, imhPml.getImagePlus(), imgDapi, imgPml};
        ImagePlus imgObjects = new RGBStackMerge().mergeHyperstacks(imgColors, false);
        imgObjects.setCalibration(cal);
        new FileSaver(imgObjects).saveAsTiff(outDir + imgName + ".tif");
        
        imhNuc.closeImagePlus();
        imhPml.closeImagePlus();
        closeImage(imgObjects);
    }
}
