/*
 * Analyze DNA or PML in nucleus
 * 
 * Author Philippe Mailly
 */


/* 
* Images on local machine
*/


import DNA_Utils.DNA_Processing;
import ij.*;
import ij.measure.Calibration;
import ij.plugin.PlugIn;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.logging.Level;
import java.util.logging.Logger;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.common.services.ServiceFactory;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.formats.services.OMEXMLService;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom.Objects3DPopulation;
import java.util.ArrayList;
import java.util.Collections;
import loci.plugins.BF;
import loci.plugins.in.ImporterOptions;
import java.io.FilenameFilter;
import mcib3d.geom.Object3D;
import mcib3d.geom2.Object3DInt;
import mcib3d.geom2.Objects3DIntPopulation;
import mcib3d.geom2.measurements.MeasureIntensity;
import mcib3d.geom2.measurements.MeasureVolume;
import org.apache.commons.io.FilenameUtils;


public class DNA_PML_Morph implements PlugIn {

    private String imageDir = "";
    public static String outDirResults = "";
    private static Calibration cal = new Calibration();
    private File inDir;
    private String[] chsName;
    public BufferedWriter outPutDNAResultsGlobal; 
    
    DNA_Processing proc = new DNA_Processing();

    /**
     * 
     * @param arg
     */
    @Override
    public void run(String arg) {
        
        final boolean canceled = false;
        try {
            if (canceled) {
                IJ.showMessage(" Pluging canceled");
                return;
            }
            if (!proc.checkInstalledModules()) {
                IJ.showMessage(" Pluging canceled");
                return;
            }
            imageDir = IJ.getDirectory("Images folder");
            if (imageDir == null) {
                return;
            }
            inDir = new File(imageDir);
            String imageExt = proc.findImageType(inDir);
            ArrayList<String> imageFiles = proc.findImages(imageDir, imageExt);
            if (imageFiles == null) {
                return;
            }
            
            
            // Reset foreground and background
            IJ.run("Colors...", "foreground=white background=black");
            
            // create OME-XML metadata store of the latest schema version
            ServiceFactory factory;
            factory = new ServiceFactory();
            OMEXMLService service = factory.getInstance(OMEXMLService.class);
            IMetadata meta = service.createOMEXMLMetadata();
            ImageProcessorReader reader = new ImageProcessorReader();
            reader.setMetadataStore(meta);
            Collections.sort(imageFiles);
            // Find channel names , calibration
            reader.setId(imageFiles.get(0));
            cal = proc.findImageCalib(meta, reader);
            chsName = proc.findChannels(imageFiles.get(0), meta, reader, true);
            int[] channelIndex = proc.dialog(chsName);
            cal = proc.getCalib();
            if (channelIndex == null)
                return;
            
            // create output folder
            outDirResults = imageDir +"Results"+ File.separator;
            File outDir = new File(outDirResults);
            if (!Files.exists(Paths.get(outDirResults))) {
                outDir.mkdir();
            }
            /*
            * Write headers results for results file
            */
            // Global file for DNA/PML results
            String resultsName = "Results_"+proc.dna+".xls";
            String header = "ImageName\t#Nucleus\tNucleus Volume\t"+proc.dna+" dot number\t"+proc.dna+" Total IntDensity"
                                +"\t"+proc.dna+" Diffuse IntDensity\t"+proc.dna+" Total dot Volume\n";
            outPutDNAResultsGlobal = proc.writeHeaders(outDirResults+resultsName, header); 
            for (String f : imageFiles) {
                String rootName = FilenameUtils.getBaseName(f);
                reader.setId(f);
                int series = reader.getSeries();
                reader.setSeries(series);
                series = reader.getSeriesCount();  
                for (int s = 0; s < series; s++) {
                    reader.setSeries(s);
                    String seriesName = meta.getImageName(s);
                    ImporterOptions options = new ImporterOptions();
                    options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                    options.setId(f);
                    options.setSplitChannels(true);
                    options.setQuiet(true);
                    options.setSeriesOn(s, true);
                    
                        
                    /*
                    * Open channels
                    */
                    // DAPI channel
                    int dapiCh = channelIndex[0];                        
                    options.setCBegin(s, dapiCh);
                    options.setCEnd(s, dapiCh);
                    System.out.println("-- Series : "+ seriesName);
                    System.out.println("Opening Nucleus channel");
                    ImagePlus imgNucOrg= BF.openImagePlus(options)[0]; 
                    // Find nucleus
                    Objects3DIntPopulation nucPop = proc.cellposeDetection(imgNucOrg);
                    int totalNucPop = nucPop.getNbObjects();
                    System.out.println("Detected nucleus after filtering = "+totalNucPop);
                        
                    // DNA/PML channel
                    int dnaPmlCh = channelIndex[1];
                    options.setCBegin(s, dnaPmlCh);
                    options.setCEnd(s, dnaPmlCh);System.out.println("Opening dna/pml channel");
                    ImagePlus imgDnaPmlOrg = BF.openImagePlus(options)[0];
                     // Find Dna/Pml
                    Objects3DIntPopulation dnaPmlPop = proc.stardistPop(imgDnaPmlOrg, proc.dna);
                    System.out.println("Total "+proc.dna+" pop = "+ dnaPmlPop.getNbObjects()); 
                    
                    // Find dna/pml inside all nucleus
                    Objects3DIntPopulation dnaPmlInAllNucleus = proc.findFociInNuclei(dnaPmlPop, nucPop, imgNucOrg);
                    System.out.println("DNA/PML pop inside nucleus = "+ dnaPmlInAllNucleus.getNbObjects()); 
                    // Find dna/pml inside each nucleus
                    for (Object3DInt nucObj : nucPop.getObjects3DInt()) {
                        double nucVol = new MeasureVolume(nucObj).getVolumeUnit();
                        Objects3DIntPopulation dnaPmlInNucleus = proc.fociInNucleus(dnaPmlInAllNucleus, nucObj);
                        double fociVol = proc.findPopVolume(dnaPmlInNucleus);
                        double diffusInt = proc.getDiffus(imgDnaPmlOrg, dnaPmlInNucleus, nucObj);
                        double fociInt = proc.findPopIntensity(dnaPmlInNucleus, imgDnaPmlOrg);
                        outPutDNAResultsGlobal.write(rootName+"\t"+nucObj.getLabel()+"\t"+nucVol+"\t"+dnaPmlInNucleus.getNbObjects()+"\t"+fociInt+"\t"
                                +diffusInt+"\t"+fociVol+"\n");
                        outPutDNAResultsGlobal.flush();
                    }   

                    proc.saveObjects(imgNucOrg, nucPop, dnaPmlInAllNucleus, outDirResults+rootName+"_"+proc.dna+"_Objects.tif");

                    proc.flush_close(imgNucOrg);
                    proc.flush_close(imgDnaPmlOrg);
                    options.setSeriesOn(s, false);
                }
            }
            outPutDNAResultsGlobal.close();
            IJ.showStatus("Process done");
        }   catch (IOException | DependencyException | ServiceException | FormatException ex) {
            Logger.getLogger(DNA_PML_Morph.class.getName()).log(Level.SEVERE, null, ex);
        } 
    }
}