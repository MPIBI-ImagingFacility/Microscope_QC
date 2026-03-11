# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 09:13:01 2025

After following the Leica_Confocal_QC_Protocol, this script can be used to make
plots summarizing the collected data. 
Requires "QC_functions.py"

@author: dpaynter
"""

import os
import QC_functions as qcf
import tkinter as tk
from tkinter.filedialog import askopenfilename
tk.Tk().withdraw() # part of the import if you are not using other tkinter functions

############################################################
############################################################

# Input the path to QC folder for a particular date:
qc_path = r"J:\Equipment\Microscopes\Imaging Facility\Leica_SP8_Matrix_N132\QualityControl\250305"

qc_folder = os.path.dirname(qc_path)
date_ = qc_path.split("\\")[-1]
microscope = qc_path.split("\\")[-3]

#assumed settings
hyd_gains = [20, 40, 60, 80, 100]
pmt_gains= [200, 400, 600, 800, 1000]
pmt_gains_new = [200, 300, 400, 500, 600, 700, 800, 900, 1000]

print(f"Running QC analysis and plotting for {microscope} on {date_}.")


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
gainlin_files = [file for file in os.listdir(qc_path) if "gain_lin" in file and ".lif" in file]
if len(gainlin_files) != 1:
    print("Can't find gain linearity file. Manually input the correct one:")
    filepath = askopenfilename()
elif len(gainlin_files) == 1:
    filepath = os.path.join(qc_path, gainlin_files[0])

# run function
try:
    qcf.plot_gain_lin_leica(gain_lif_path=filepath, hyd_gains=hyd_gains, pmt_gains=pmt_gains, microscope=microscope, date=date_)
except:
    qcf.plot_gain_lin_leica(gain_lif_path=filepath, hyd_gains=hyd_gains, pmt_gains=pmt_gains_new, microscope=microscope, date=date_)

############################################################
## Look at detector gain over time

qcf.plot_gain_lin_overtime(qc_folder)

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

# make a plot for each field illumination acquisition 
for file in fieldillumfiles:
    imagepath = os.path.join(qc_path, file)
    qcf.plot_field_illum(imagepath, date=date_, microscope=microscope)

############################################################
## look at stage x-y repositioning
pos1path = os.path.join(qc_path, 'pos1.tif')
pos2path = os.path.join(qc_path, 'pos2.tif')
if os.path.isfile(pos1path) and os.path.isfile(pos2path):
    qcf.plot_xy_repos(pos1path, pos2path)
else:
    print("Can't find stage repositioning files.n\ Files should be named 'pos1.tif' and 'pos2.tif'.")






