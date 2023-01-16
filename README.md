# DNA_PML_MORPH

* **Developed for:** Domitille
* **Team:** De Thé
* **Date:** January 2023
* **Software:** Fiji

### Images description

3D images taken with a x63 objective

2 channels:
  1. *Alexa Fluor 405:* DAPI 
  2. *Alexa Fluor 488:* DNA/PML

### Plugin description

* Detect nuclei with Cellpose
* Find DNA/PML foci with Stardist
* For each nucleus, count foci number, volume and intensity
* For each nucleus, compute also diffuse intensity (ie. nucleus intensity in DNA/PML channel after setting intensity of detected DNA/PML foci to zero)

### Dependencies

* **3DImageSuite** Fiji plugin
* **CLIJ** Fiji plugin
* **Cellpose** conda environment + *cyto2* model
* **Stardist** conda environment + model selected by user

### Version history

Version 1 released on January 13, 2023.
