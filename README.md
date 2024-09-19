# DAPI_PML

* **Developed by:** Héloïse
* **Developed for:** Adèle
* **Team:** De Thé
* **Date:** September 2024
* **Software:** Fiji

### Images description

3D images taken with a x63 objective.

2 channels:
  1. *DAPI:* DAPI nuclei
  2. *Alexa Fluor 488:* PML foci

### Plugin description

* Detect DAPI nuclei with Cellpose
* Segment PML foci with DoG filtering + automatic thresholding
* Compute background noise as the median intensity of pixels outside nuclei
* For each nucleus, provide:
  * PML foci number, area and background-corrected intensity
  * PML diffuse signal (aka outside foci) background-corrected intensityDetect nuclei with Cellpose

### Dependencies

* **3DImageSuite** Fiji plugin
* **CLIJ** Fiji plugin
* **Cellpose** conda environment + *cyto2* model

### Version history

Version 2 released on September 19, 2024.
