# -*- coding: utf-8 -*-
"""
Created on Fri Jul 25 16:52:20 2025

For processing QC data acquired at the Nikon Spinning Disk

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
qc_path = r"J:\Equipment\Microscopes\Imaging Facility\Nikon_SpinningDisk_N135\QualityControl\260310"

qc_folder = os.path.dirname(qc_path)
date_ = qc_path.split("\\")[-1]
microscope = qc_path.split("\\")[-3]

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
## Camera dark counts
# find the darkcounts files for each camera
darkcountfiles = [file for file in os.listdir(qc_path) if "darkcounts" in file]

# make a plot for each camera using function in 'QC_functions.py'
for file in darkcountfiles:
    qcf.plot_camera_darknoise(filepath=os.path.join(qc_path, file), save_csv=True, microscope=microscope, date=date_)


############################################################
## Field Illumination
fieldillumfiles = [file for file in os.listdir(qc_path) if "field_illum" in file and ".nd2" in file]
if len(fieldillumfiles) != 0:
    # make a plot for each field illumination acquisition 
    for file in fieldillumfiles:
        imagepath = os.path.join(qc_path, file)
        qcf.plot_field_illum(imagepath, date=date_, microscope=microscope)
else:
    print("Can't find field_illum files.n\ Files should be named 'field_illum_10x...'")


############################################################
## XY stage calibration
pos1path = os.path.join(qc_path, 'pos1.tif')
pos2path = os.path.join(qc_path, 'pos2.tif')
if os.path.isfile(pos1path) and os.path.isfile(pos2path):
    qcf.plot_xy_repos(pos1path, pos2path)
else:
    print("Can't find stage repositioning files.n\ Files should be named 'pos1.tif' and 'pos2.tif'.")


