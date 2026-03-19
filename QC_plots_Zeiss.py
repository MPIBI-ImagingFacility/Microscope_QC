# -*- coding: utf-8 -*-
"""
Created on Fri Jul 25 16:52:20 2025

For processing QC data acquired at the Zeiss LSM900

@author: dpaynter
"""

import os
import QC_functions as qcf
import tkinter as tk
from tkinter.filedialog import askopenfilename
tk.Tk().withdraw() # part of the import if you are not using other tkinter functions
from pathlib import Path

############################################################
############################################################

# Input the path to QC folder for a particular date:
qc_path = r"J:\Equipment\Microscopes\Imaging Facility\Zeiss_LSM900_O104\QualityControl\260303"

qc_folder = os.path.dirname(qc_path)
date_ = qc_path.split("\\")[-1]
microscope = qc_path.split("\\")[-3]

#assumed settings
pmt_gains= [400, 500, 600, 700, 800, 900, 1000]

print(f"Running QC analysis and plotting for {microscope} on {date_}.")

############################################################
############################################################
## Look at  today's laser power 
# find file
laserpower_files = [file for file in os.listdir(qc_path) if "laserpower.xl" in file]
if len(laserpower_files) > 1:
    print("Multiple laser power files found. Manually input the correct one:")
    filepath = askopenfilename()
elif len(laserpower_files) == 1:
    filepath = os.path.join(qc_path, laserpower_files[0])

# run function 
qcf.plot_laserpower(filepath,  microscope=microscope, date=date_)
############################################################
## Look at laser power over time

qcf.plot_laserpower_overtime(qc_folder)

############################################################
## look at detector_gain_linearity (does the measured intensity increase linearly as gain increases linearly?)
# Find the file
gainlin_files = [file for file in os.listdir(qc_path) if "gain_lin" in file and ".tif" in file]
if len(gainlin_files) != 1:
    print("Can't find gain linearity file. Manually input the correct one:")
    filepath = askopenfilename()
elif len(gainlin_files) == 1:
    filepath = os.path.join(qc_path, gainlin_files[0])

# run function
qcf.plot_gain_lin_tif(gain_tif_path=filepath, detector_gains=pmt_gains, microscope=microscope, date=date_)

############################################################
## look at detector_darknoise 
# Find the file
darknoise_files = [file for file in os.listdir(qc_path) if "darknoise" in file]
if len(darknoise_files) != 1:
    print("Can't find detector darknoise. Manually input the correct one:")
    filepath = askopenfilename()
elif len(darknoise_files) == 1:
    filepath = os.path.join(qc_path, darknoise_files[0])

# run function 
qcf.plot_detector_darknoise_leica(filepath)

############################################################
## look at field_illumination
fieldillumfiles = [file for file in os.listdir(qc_path) if "field_illum" in file and ".tif" in file]
if len(fieldillumfiles) == 0:
    print("Didn't find field_illum TIFs. Check that name contains 'field_illum' and '.tif.'")
# make a plot for each field illumination acquisition 
for file in fieldillumfiles:
    imagepath = os.path.join(qc_path, file)
    qcf.plot_field_illum(imagepath, date=date_, microscope=microscope)


