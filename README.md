# Microscope_QC

Protocols for collecting data, and Python/MATLAB code for analyzing data for quality control at MPI-BI microscopes. 

For questions, comments, and suggestions, please reach out to: imaging@bi.mpg.de

## Protocols

The protocols in this repository are written specifically for the microscopes at MPI-BI, but may be adapted for use elsewhere.
Currently, protocols are for:

Leica: STELLARIS 5, SP8 confocals

Zeiss: LSM900 confocal

Nikon: CSU-W1 SoRa spinning disk

## Analysis

### Installation
#### For Python: Clone the repository and install dependencies:

git clone https://github.com/MPIBI-ImagingFacility/Microscope_QC.git

cd microscope_qc

pip install -r requirements.txt

#### For MATLAB (PSF analysis):
You will need a MATLAB license. Code runs in MATLAB R2023a; newer versions may not work.

You will also need to download the bioformats_package.jar file and move it to the MATLAB_PSF local file, using this link: https://www.openmicroscopy.org/bio-formats/downloads/




### Usage

Open a Python interpreter:

- Run the script containing functions for analysis:
python QC_functions.py

- Open the "QC_plots_xxx" file corresponding to the protocol used (currently: Leica, Zeiss, or Nikon)

- Edit the line below to match the path to a specific date's data:

#Input the path to QC folder for a particular date:

qc_path = r"C:\Microscopes\QualityControl\251217"

- Run the code. It will save plots to the qc_path (single time point measures) as well as to the directory above (comparing measurements over time).

- Functions in QC_functions.py can also be run individually for processing of specific images, rather than full QC datasets.

Open MATLAB:

- Add the folder Microscope_QC/MATLAB_PSF to the Path in your MATLAB instance.

- Run the following line of code for .lif PSF files:

psfevalleica('file_with_one_PSF_stack.lif',1, 1, 1.7, 1, 1, 0, 1);

% psfevalleica(filename, showPlot, zoomIn, redThreshold, showMicroscopeEvalPlot, showZProjection, beadSize, extractMetaDataFromPath)

- Or this line, for other file types: 

psfeval('file_with_one_PSF_stack.nd2',1,1, 1.7, 1, 1, 0, 1);

% or .tif, .czi, etc.

### Requirements
Tested with:

Python 3.11–3.12

MATLAB R2023a

Dependencies for Python code are listed in requirements.txt.

Dependencies for MATLAB code are contained in the MATLAB_PSF folder.
