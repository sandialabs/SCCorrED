# SCCorrED

Sandia Cross Correlation for Elastic Deformations (SCCorrED) is a MATLAB code for High Resolution Electron Back Scatter Diffraction (HR-EBSD) analysis developed at Sandia National Laboratories. 

## Description

SCCorrED is meant as an HR-EBSD alternative compatible with the popular MTEX EBSD software with more modularity to encourage further development of the technique. This is a work in progress and there are many features planned that are still being implemented, however the main functionality of HREBSD analysis for elastic strains and bulk GND density have been implemented.  

## Getting Started

### Dependencies

* MATLAB R2021b or later
* MTEX
* Signal processing toolbox
* Parallel processing toolbox (optional, required only for parallel processing)
* EMSoft (optional, required for pattern matching for absolute strain calculations)

### Installing

* Clone/download this repository
* Create a mtexHREBSD directory in `<`user directory`>`/.config directory (may need to create .config) and create "mtexHREBSDConfig.json" file with the following:
'''
{
	"tensorPath": "<location of MTEX's tensor directory>",
	"EMdataPath": "<location of EMSoft's data directory>", 
	"EMsoftPath": "<location of EMSoft's main directory>",
	"materialPath": "<location of MTEX's material directory>"
}
'''

### Executing program

Two live scripts have been provided in the "Example" directory, along with a small data set, one for single processing and one for parallel processing:
* exampleDriver.mlx
* exampleDriverParallel.mlx
The relevant live script will show how to execute the code. 


## Authors

[Will Gilliland](wggilli@sandia.gov)
[Tim Ruggles](truggle@sandia.gov)
[Thomas Bennett](tbenne@sandia.gov) 

## Version History

* 0.1
    * Initial Release

## License

This project is licensed under the MIT License - see the LICENSE.txt file for details.  The files in Code/+patterns concerning pattern handling are from OpenXY, another open-source HR-EBSD package.  The MIT license for OpenXY is included in that folder as well.
