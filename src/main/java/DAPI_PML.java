import DAPI_PML_Tools.Nucleus;
import DAPI_PML_Tools.Tools;
import ij.*;
import ij.plugin.PlugIn;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.SimpleDateFormat;
import java.util.logging.Level;
import java.util.logging.Logger;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.plugins.util.ImageProcessorReader;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import loci.formats.MetadataTools;
import loci.plugins.BF;
import loci.plugins.in.ImporterOptions;
import mcib3d.geom2.Objects3DIntPopulation;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.ArrayUtils;

/**
 * @author Héloïse Monnet @ ORION-CIRB
 */
public class DAPI_PML implements PlugIn {
    
    private Tools tools = new Tools();

    public void run(String arg) {
        try {
            if (!tools.checkInstalledModules()) {
                return;
            }
            
            String imgDir = IJ.getDirectory("Choose images folder");
            if (imgDir == null) {
                return;
            }
            
            // Find extension of first image in input folder
            String fileExt = tools.findImageType(new File(imgDir));
            // Find all images with corresponding extension in folder
            ArrayList<String> imgFiles = tools.findImages(imgDir, fileExt);
            if (imgFiles.isEmpty()) {
                IJ.showMessage("ERROR", "No image found with " + fileExt + " extension in " + imgDir + " folder");
                return;
            }
            
            // Instantiate metadata and reader
            IMetadata meta = MetadataTools.createOMEXMLMetadata();
            ImageProcessorReader reader = new ImageProcessorReader();
            reader.setMetadataStore(meta);
            reader.setId(imgFiles.get(0));
            
            // Find image calibration
            tools.findImageCalib(meta);

            // Find channels name
            String[] chMeta = tools.findChannels(imgFiles.get(0), meta, reader);
            
            // Generate dialog box
            String[] chOrder = tools.dialog(chMeta);
            if (chOrder == null) {
                return;
            }
            
            // Create output folder
            String outDir = imgDir + File.separator + "Results_" + new SimpleDateFormat("yyyy-MM-dd_HH-mm-ss").format(new Date()) + File.separator;
            if (!Files.exists(Paths.get(outDir))) {
                new File(outDir).mkdir();
            }
           
            // Write header in results file
            FileWriter fwResults = new FileWriter(outDir + "results.csv", false);
            BufferedWriter results = new BufferedWriter(fwResults);
            results.write("Image name\tPML background noise\tNucleus ID\tNucleus area (µm2)\tPML foci number"
                          + "\tPML foci total area (µm2)\tPML foci bg-corr mean intensity\tPML foci bg-corr raw integrated density"
                          + "\tPML diffuse area (µm2)\tPML diffuse bg-corr mean intensity\tPML diffuse bg-corr raw integrated density\n");
            results.flush();
            
            for (String file: imgFiles) {
                String imgName = FilenameUtils.getBaseName(file);
                System.out.println("--- ANALYZING IMAGE " + imgName + " ---");
                reader.setId(file);
                
                ImporterOptions options = new ImporterOptions();
                options.setId(file);
                options.setSplitChannels(true);
                options.setQuiet(true);
                options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                    
                // Open DAPI channel
                tools.print("- Opening DAPI channel -");
                int chIndex = ArrayUtils.indexOf(chMeta, chOrder[0]);
                ImagePlus imgDapi = BF.openImagePlus(options)[chIndex];
                
                // Open PML channel
                tools.print("- Opening PML channel -");
                chIndex = ArrayUtils.indexOf(chMeta, chOrder[1]);
                ImagePlus imgPml = BF.openImagePlus(options)[chIndex];
                
                // Detect DAPI nuclei
                tools.print("- Detecting DAPI nuclei -");
                ImagePlus maskDapi = tools.cellposeDetection(imgDapi);
                double pmlBg = tools.computeBackground(maskDapi, imgPml);
                Objects3DIntPopulation nucPop = tools.filterPop(maskDapi);
                
                // Detect PML foci
                tools.print("- Detecting PML foci -");
                Objects3DIntPopulation pmlPop = tools.fociDetection(imgPml);

                // Get nuclei with their respective population of PML foci
                tools.print("- Getting PML foci for each nucleus -");
                ArrayList<Nucleus> nuclei = tools.colocalizeNucFoci(nucPop, pmlPop);
               
                // Write results
                tools.print("- Writing results -");
                for (Nucleus nucleus: nuclei) {
                    HashMap<String, Double> p = nucleus.computeParams(imgPml, pmlBg, tools.pixArea);
                    results.write(imgName+"\t"+pmlBg+"\t"+p.get("nucId").intValue()+"\t"+p.get("nucArea")+"\t"+p.get("pmlFociNb").intValue()+"\t"+
                                  p.get("pmlFociArea")+"\t"+p.get("pmlFociMeanInt")+"\t"+p.get("pmlFociTotInt")+"\t"+
                                  p.get("pmlDiffuseArea")+"\t"+p.get("pmlDiffuseMeanInt")+"\t"+p.get("pmlDiffuseTotInt")+"\n");
                    results.flush();
                }
                
                // Draw results
                tools.print("- Drawing results -");
                tools.drawResults(nuclei, imgPml, imgDapi, outDir, imgName);

                tools.closeImage(imgDapi);
                tools.closeImage(imgPml);
                tools.closeImage(maskDapi);
            }
            results.close();
            tools.print("--- All done! ---");
        } catch (DependencyException | ServiceException | IOException | FormatException ex) {
            Logger.getLogger(DAPI_PML.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}