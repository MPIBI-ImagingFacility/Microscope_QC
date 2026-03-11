# Microscope_QC

Protocols and Python scripts for quality control analysis at MPI-BI microscopes. 

For questions, comments, and suggestions, please reach out to: imaging@bi.mpg.de

## Protocols

The protocols in this repository are written specifically for the microscopes at MPI-BI, but may be adapted for use elsewhere.
Currently, protocols are for:

Leica: STELLARIS 5, SP8 confocals

Zeiss: LSM900 confocal

Nikon: CSU-W1 SoRa spinning disk

## Analysis

### Installation
Clone the repository and install dependencies:

git clone https://github.com/MPIBI-ImagingFacility/Microscope_QC.git

cd microscope_qc

pip install -r requirements.txt

### Usage

Run the script containing functions for analysis:
python QC_functions.py

Open the "QC_plots_xxx" file corresponding to the protocol used (currently: Leica, Zeiss, or Nikon)

Edit the line below to match the path to a specific date's data:

#Input the path to QC folder for a particular date:

qc_path = r"C:\Microscopes\QualityControl\251217"

Run the code. It will save plots to the qc_path (single time point measures) as well as to the directory above (comparing measurements over time).

Functions in QC_functions.py can also be run individually for processing of specific images, rather than full QC datasets.

### Requirements
Tested with:

Python 3.11–3.12

Dependencies are listed in requirements.txt.
