# -*- coding: utf-8 -*-
"""
Created on Fri Aug 22 11:59:30 2025

Contains functions useful for analyzing data collected for quality control purposes.

@author: dpaynter
"""

### Imports
import os
import pandas as pd
import numpy as np
from skimage.feature import peak_local_max
from readlif.reader import LifFile
import matplotlib.pyplot as plt
import seaborn as sns
import cv2
from nd2reader import ND2Reader
import csv
from scipy.stats import linregress
import multipagetiff as mtif
import matplotlib
import tkinter as tk
from tkinter.filedialog import askopenfilename
tk.Tk().withdraw() # part of the import if you are not using other tkinter functions
import tifffile
from pathlib import Path
import matplotlib.colors as mcolors
import glob



### Presets and hardcoded things:
# palette is used for assigning colors to laser power plots. If using a laser not listed here, simply add it to the list and copy/paste a color name to use for it.
palette = {'405': 'black', '442': 'deepskyblue', '445': 'deepskyblue','448': 'deepskyblue', 'Argon_488': 'limegreen', 'OPSL_488': 'limegreen', 'WLL_485':'limegreen','488': 'limegreen', 'Argon_514':'gold','Argon_561':'pink', '633':'darkred', '561':'red', 'Argon_458':'aquamarine',
           'WLL_488':'limegreen', 'WLL_514':'gold', '514':'gold','515':'gold', 'WLL_561':'red', 'WLL_594':'lightsalmon', 'WLL_633':'palevioletred', '638':'palevioletred','WLL_640':'firebrick','640':'firebrick', 'WLL_670':'maroon', 'WLL_685':'maroon', '730':'maroon'}



# settings for retaining pdf font
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['font.size'] = 11
matplotlib.rcParams['patch.edgecolor'] = 'none'
# seaborn color context (all to black)
color_context = { 'axes.edgecolor': '#000000',
'axes.labelcolor': '#000000',
'boxplot.capprops.color': '#000000',
'boxplot.flierprops.markeredgecolor': '#000000',
'grid.color': '#000000',
'patch.edgecolor': '#000000',
'text.color': '#000000',
'xtick.color': '#000000',
'ytick.color': '#000000'}
marker_style = dict(linestyle=':', color='0.8', markersize=10,
                    markerfacecolor="None", markeredgecolor="#000000")


# FUNCTIONS

def init_figure_axes( fig_size=(10,7), dpi=80, facecolor="w", edgecolor="w" ):
# Convert fig size to inches (default is inches, fig_size argument is supposed to be in cm)
    inch2cm = 2.54
    fig_size = fig_size[0]/inch2cm,fig_size[1]/inch2cm
    with sns.axes_style(style="ticks",rc=color_context):
        fig,ax = plt.subplots(num=None, figsize=fig_size, dpi=dpi,
                              facecolor=facecolor, edgecolor=edgecolor)
        ax.spines['bottom'].set_linewidth(.8)
        ax.spines['left'].set_linewidth(.8)
        ax.tick_params(width=.8)
        return fig,ax

def find_bead_position(image, threshold_rel=0.5):
    """
    Return the (row, column) coordinates of the brightest bead in an image.

    The image is first Gaussian-smoothed to reduce noise, then local maxima
    are detected. The brightest detected peak is returned, assuming the
    image contains a single bright bead.

    Parameters
    ----------
    image : ndarray
        2D grayscale image (e.g. loaded with OpenCV) containing one bright bead.
    threshold_rel : float, optional
        Relative intensity threshold for peak detection, expressed as a
        fraction of the maximum image intensity. Default is 0.5.

    Returns
    -------
    tuple or None
        (row, column) coordinates of the brightest detected bead.
        Returns None if no peaks are found above the threshold.
    """
    # Smooth to reduce noise
    blurred = cv2.GaussianBlur(image, (3, 3), 0)

    # Get local maxima
    coordinates = peak_local_max(blurred, min_distance=5, threshold_rel=threshold_rel)
    
    if len(coordinates) == 0:
        return None
    return coordinates[0]  # Return the brightest peak (assumes 1 bead)


def track_beads(image_file):

    """
    Track the position of a bright bead across frames in an image stack.

    The image stack is loaded from disk, and for each frame the position
    of the brightest bead is detected using `find_bead_position`.

    Parameters
    ----------
    image_file : str or Path
        Path to a multi-frame image file (e.g. TIFF stack) containing a
        single bright bead per frame.

    Returns
    -------
    ndarray
        Array of bead positions with shape (n_frames, 2), where each entry
        is (row, column). Frames in which no bead is detected contain None.
    """
    positions = []

    img = mtif.read_stack(image_file)
    
    for frame in img:    
        pos = find_bead_position(frame)
        positions.append(pos)
    return np.array(positions)

def compute_displacements(positions):
    """
   Compute the Euclidean displacement of each position from the mean reference position.

   The function calculates the average (centroid) of a set of positions and then
   computes the Euclidean distance of each individual position from this reference
   point.

   Parameters
   ----------
   positions : numpy.ndarray
       Array of shape (N, D) containing N positions in D-dimensional space
       (commonly D=2 or D=3).

   Returns
   -------
   numpy.ndarray
       A 1D array of length N containing the Euclidean distance of each position
       from the mean reference position.
   """
    
    reference = np.mean(positions, axis=0)
    displacements = np.linalg.norm(positions - reference, axis=1)
    return displacements


def plot_xy_repos(pos1_path=None, pos2_path=None, pixel_size_um=0.1):
    
    """
    Analyze and plot stage XY repositioning error using tracked bead positions.
    
    This function loads two datasets of bead positions, computes their displacement
    from the mean position for each repetition, converts the displacement from pixels
    to micrometers, and plots the repositioning error across repetitions. Summary
    statistics (mean ± standard deviation) are printed and displayed on the plot.
    The resulting figure is saved in the same directory as the first input file.
    
    Parameters
    ----------
    pos1_path : str, optional
        Path to the first bead position file. If None, a file selection dialog
        will open to allow the user to choose the file.
    
    pos2_path : str, optional
        Path to the second bead position file. If None, a file selection dialog
        will open to allow the user to choose the file.
    
    pixel_size_um : float, optional
        Pixel size in micrometers used to convert displacement values from
        pixels to micrometers. Default is 0.1 µm per pixel.
    
    Returns
    -------
    None
        Displays a plot of displacement vs. repetition for both beads and saves
        the figure as "stage_repositioning.png" in the same folder as `pos1_path`.
    
    Notes
    -----
    - Requires the functions `track_beads()` and `compute_displacements()`.
    - Assumes each input file contains repeated measurements of the same bead
      after stage repositioning.
    """
        
    if pos1_path == None:
        pos1_path = askopenfilename()
    if pos2_path == None:
        pos2_path = askopenfilename()
    
    folder = os.path.dirname(pos1_path)
    
    positions_b1 = track_beads(pos1_path)
    positions_b2 = track_beads(pos2_path)
    disp_b1 = compute_displacements(positions_b1)
    disp_b2 = compute_displacements(positions_b2)
    disp_b1_um = disp_b1 * pixel_size_um
    disp_b2_um = disp_b2 * pixel_size_um
    
    mean1 = np.mean(disp_b1_um)
    mean2 =  np.mean(disp_b2_um)
    std1 = np.std(disp_b1_um)
    std2 = np.std(disp_b2_um)
    print("Bead 1 mean ± std displacement (um):", np.mean(disp_b1_um), "±", np.std(disp_b1_um))
    print("Bead 2 mean ± std displacement (um):", np.mean(disp_b2_um), "±", np.std(disp_b2_um))

    fig, ax = plt.subplots(figsize=(6,4), dpi=96, facecolor="w", edgecolor="w")
    plt.plot(disp_b1_um, label='Bead 1')
    plt.plot(disp_b2_um, label='Bead 2')
    plt.ylabel('Displacement from mean (um)')
    plt.xlabel('Repetition')
    plt.legend()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.title('Stage Repositioning Error')
    # Add text to the plot
    textstr = f'Bead 1: {mean1:.3f} ± {std1:.3f} µm\nBead 2: {mean2:.3f} ± {std2:.3f} µm'
    plt.gcf().text(0.15, 0.75, textstr, fontsize=10, bbox=dict(facecolor='white', alpha=0.7))
    plt.savefig(os.path.join(folder, "stage_repositioning.png"))
    plt.show()
    plt.close()
    
    
def lif_to_stack(filename, series_index=0, channel=0, timepoint=0):
    """
    Load a (Z, Y, X) stack from a LIF file for a given series, channel, and timepoint.
    Function returns a numpy array with dimensions Z, Y, X. Used in case LIF file cannot be managed by other Python packages.
    "filename" should be a full path to a LIF file. Example: r"J:\\Equipment\\Microscopes\\Imaging Facility\\Leica_SP8_N131\\QualityControl\\250425\\field_illumination.lif"
    Using r"text" notation allows for reading through back slashes without reading it as unicode.
    
    """
    reader = LifFile(filename)
    im = reader.get_image(0)
    nx, ny, nz, nc, nt = im.dims

    # Allocate empty array
    stack = np.empty((nz, ny, nx), dtype=np.uint16)  # adjust dtype if needed

    # Iterate over z planes
    for z, plane in enumerate(im.get_iter_z(t=timepoint, c=channel)):
        stack[z] = plane

    return stack


def plot_laserpower(filepath=None, date=None, microscope=None):
    """Function to plot the laser power measured on a given day, with percentage of maximum power on x-axis, 
        and power in microwatts on y-axis.
        laserpower_file_path should point to an Excel sheet (.xlsx) formatted with three columns, named:
        laser, percentage, power"""
    
    if filepath == None:
        filepath = askopenfilename()
    if microscope == None:
        microscope = input("Enter name of microscope:  ")
    if date==None:
        date = input("Enter date of acquisition (YYYY-MM-DD):  ")
   

    laserpower_data = pd.read_excel(filepath)
    laserpower_data['percentage'] = pd.to_numeric(laserpower_data['percentage'], errors='coerce')
    
    
    # Force laser names to strings
    laserpower_data['laser'] = laserpower_data['laser'].astype(str)
    unique_lasers = laserpower_data['laser'].unique()
    
    style_dashes = {laser: (5, 1) if "WLL" in laser else '' for laser in unique_lasers}
    palette_small = {laser: palette[laser] for laser in unique_lasers}
    
    fig, ax = plt.subplots(figsize=(6,5), dpi=96, facecolor="w", edgecolor="w")
    ax = sns.lineplot(data=laserpower_data, x='percentage', y='power', hue='laser', palette=palette_small, style='laser', dashes=style_dashes)
    plt.xlabel("Percentage")
    plt.ylabel("Power (uW)")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylim(bottom=0)
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
    plt.title(f"Laser power at {microscope} on {date}")
    
    folder = os.path.dirname(filepath)
    savepath = os.path.join(folder, "laser_power.png")
    plt.savefig(savepath, bbox_inches='tight')
    plt.show()
    plt.close()



def plot_laserpower_overtime(qc_folder_path):
    """Function to plot the laser power for a given system measured at several points in time, with date on x-axis, 
        and power in microwatts on y-axis. Only laser power at 100% of maximum power is plotted.
        qc_folder_path should point to a folder specific to the microscope. It looks for subfolders labelled as dates ("250416","250813", for example).
        Function finds one .xlsx file with "laser" in the file name, for each date.
        Excel sheets should be formatted with three columns: laser, percentage, power"""
        
    power_over_time = pd.DataFrame(columns=["date", "laser", "power"])
    
    # find all the laser power files
    dates = [date for date in os.listdir(qc_folder_path) if 'laser' not in date]
    
    for date in dates:
        try:
            laserpower_file = [file for file in os.listdir(os.path.join(qc_folder_path, date)) if 'laser' in file and '.xlsx' in file][0]
            
            df = pd.read_excel(os.path.join(qc_folder_path, date, laserpower_file))
            df_100 = df[df['percentage'] == 100]
            df_100['laser'] = df_100['laser'].astype(str)
    
            # Now build a small dataframe and append
            for idx, row in df_100.iterrows():
                power_over_time = pd.concat([
                    power_over_time,
                    pd.DataFrame({
                        'date': [date],
                        'laser': [row['laser']],
                        'power': [row['power']]
                    })
                ], ignore_index=True)
    
        except Exception as e:
            print(fr"No laser power file on {date}: {e}")
    
    
    #power_over_time['date'] = pd.to_datetime(power_over_time['date'], errors='coerce')
    power_over_time = power_over_time.sort_values(by='date')
    power_over_time['power'] = pd.to_numeric(power_over_time['power'], errors='coerce')
    
    unique_lasers = power_over_time['laser'].unique()
    # Make the WLL lines dashed instead of solid
    style_dashes = {laser: (5, 1) if "WLL" in laser else '' for laser in unique_lasers}
    
    # Get a subset of the global "palette" variable, containing only lasers measured in current dataset
    palette_small = {laser: palette[laser] for laser in unique_lasers}
    
    # Make a figure
    fig, ax = plt.subplots(figsize=(6,4), dpi=96, facecolor="w", edgecolor="w")
    #Plot the data
    ax = sns.lineplot(data=power_over_time, x='date', y='power', hue='laser', style='laser', palette=palette_small, dashes=style_dashes)
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_title("Laser power over time")
    ax.set_ylabel("Power (uW)")
    ax.set_xlabel("Date")
    ax.set_ylim(bottom=0)
    plt.xticks(rotation=45, ha='right')
    
    plt.tight_layout()
    # Save the figure in the same place where the date folders are found
    plt.savefig(os.path.join(qc_folder_path, "laserpower_overtime.png"))


def plot_detector_darknoise_leica(filepath=None):
    
    """Function is used to plot "darknoise" images from PMTs and HyDs.
    Currently very hard-coded to take files from Leica confocals, where each 
    channel is a different detector, and the order of acquisition matters (i.e., 
    detector names are not read in from metadata, but are assumed).
    Filepath should be something like: r"J:\\Equipment\\Microscopes\\Imaging Facility\\Leica_SP8_N131\\QualityControl\\250321\\detector_darknoise.lif" 
    Function saves a PNG plotting the dark noise level for each detector."""
    if filepath == None:
        filepath = askopenfilename()
    
    microscope = input("Enter name of microscope:  ")
    date = input("Enter date of acquisition:  ")
    if ".lif" in filepath:
        darknoise_lif  = LifFile(filepath)
        
        # Extract the single image
        single_image = darknoise_lif.get_image(0)  # Assuming the first image contains all detectors
        
        if "STELLARIS" in filepath:
            # Stellaris scopes at MPI-BI don't have PMTs
            detector_names = ['HyD1', 'HyD2', 'HyD3']
        elif "SP8" in filepath:
            # This detector configuration reflects SP8 and SP*_Matrix at MPI-BI
            detector_names = ['PMT1', 'HyD2', 'PMT3', 'HyD4', 'HyD5']
        else:
            # If the scope is neither SP8 nor Stellaris, the detectors corresponding to each channel of the image can be manually entered
            detector_names = input("Enter the detector names as a list of strings, in order of acquisition (i.e., ['PMT1', 'HyD2']")
        detector_means = []
        
        for chan_count in range(single_image.channels):
            # One channel is expected to be from one detector
            chan_im = np.asarray(single_image.get_frame(c=chan_count))
            
            # Define middle half of image to exclude edges
            h, w = chan_im.shape
            h_start, h_end = h // 4, 3 * h // 4
            w_start, w_end = w // 4, 3 * w // 4
            mean_value = np.mean(chan_im[h_start:h_end, w_start:w_end])
            
            detector_means.append(mean_value)
        
        # Plot results
        fig, ax = plt.subplots(figsize=(6,4), dpi=96, facecolor="w", edgecolor="w")
        
        ax = sns.barplot(x=detector_names, y=detector_means)
        plt.xlabel("Detector")
        plt.ylabel("Mean Intensity")
        plt.title(f"Detector Dark Noise for {microscope} on {date}")
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        # Annotate bars with values
        for bar, value in zip(ax.patches, detector_means):
            height = bar.get_height()
            ax.text(
                bar.get_x() + bar.get_width()/2, 
                height + max(detector_means) * 0.01,  # small padding above bar
                f"{value:.1f}", 
                ha='center', va='bottom', fontsize=9, color='black'
            )
            
        folder = os.path.dirname(filepath)
        savepath = os.path.join(folder, 'detector_darknoise.png')
        plt.savefig(savepath)
        plt.show()
        plt.close()
        
def plot_detector_darknoise_tif(filepath=None):
    
    """Function is used to plot "darknoise" images from PMTs and HyDs.
    Currently very hard-coded to take files from Leica confocals, where each 
    channel is a different detector, and the order of acquisition matters (i.e., 
    detector names are not read in from metadata, but are assumed).
    Filepath should be something like: r"J:\\Equipment\\Microscopes\\Imaging Facility\\Leica_SP8_N131\\QualityControl\\250321\\detector_darknoise.lif" 
    Function saves a PNG plotting the dark noise level for each detector."""
    if filepath == None:
        filepath = askopenfilename()
    
    microscope = input("Enter name of microscope:  ")
    date = input("Enter date of acquisition:  ")
    if ".tif" in filepath:
        darknoise_tif  = tifffile.imread(filepath)

        if "LSM900" in filepath:
            detector_names = ['PMT1', 'PMT2']
        else:
            # If the scope is neither SP8 nor Stellaris, the detectors corresponding to each channel of the image can be manually entered
            detector_names = input("Enter the detector names as a list of strings, in order of acquisition (i.e., ['PMT1', 'HyD2']")
        detector_means = []
        
        for chan_count in range(darknoise_tif.shape[0]):
            # One channel is expected to be from one detector
            chan_im = np.asarray(darknoise_tif[chan_count])
            
            # Define middle half of image to exclude edges
            h, w = chan_im.shape
            h_start, h_end = h // 4, 3 * h // 4
            w_start, w_end = w // 4, 3 * w // 4
            mean_value = np.mean(chan_im[h_start:h_end, w_start:w_end])
            
            detector_means.append(mean_value)
        
        # Plot results
        fig, ax = plt.subplots(figsize=(6,4), dpi=96, facecolor="w", edgecolor="w")
        
        ax = sns.barplot(x=detector_names, y=detector_means)
        plt.xlabel("Detector")
        plt.ylabel("Mean Intensity")
        plt.title(f"Detector Dark Noise for {microscope} on {date}")
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        # Annotate bars with values
        for bar, value in zip(ax.patches, detector_means):
            height = bar.get_height()
            ax.text(
                bar.get_x() + bar.get_width()/2, 
                height + max(detector_means) * 0.01,  # small padding above bar
                f"{value:.1f}", 
                ha='center', va='bottom', fontsize=9, color='black'
            )
            
        folder = os.path.dirname(filepath)
        savepath = os.path.join(folder, 'detector_darknoise.png')
        plt.savefig(savepath)
        plt.show()
        plt.close()

def plot_camera_darknoise(filepath, exposure_times=None, save_csv=True, microscope=None, date=None):
    """
    Function analyzes dark noise data acquired with cameras at Nikon systems (.nd2) files.
    Input should be a path to a (multichannel, single plane) ND2 file. The name of the "channel"
    will be taken from the Optical Configuration (OC). MPI-BI Spinning Disk has OCs corresponding to 
    various exposure times used to make this QC measurement.
    """
    
    if microscope == None:
        microscope = input("Enter name of microscope:  ")
    if date==None:
        date = input("Enter date of acquisition (YYYY-MM-DD):  ")
    if ".nd2" in filepath:
        
        # file should be named "darkcounts_cam4009.nd2"; index out the camera name:
        cam = filepath[-11:-4]
        
        # Load ND2 file
        with ND2Reader(filepath) as images:
            images.bundle_axes = 'cyx'   # (channels, y, x)
            single_image = images[0]     # first frame
    
            # Extract exposure times from metadata if available
            if exposure_times is None:
                exposure_times = []
                try:            
                    for c in range(len(images.metadata['channels'])):
                        exposure_times.append(images.metadata['channels'][c])
                        
                except Exception:
                    exposure_times = [f"Ch{c+1}" for c in range(single_image.shape[0])]
    
        num_channels = single_image.shape[0]
        exposure_means = []
    
        for chan_count in range(num_channels):
            chan_im = single_image[chan_count, :, :]
    
            # Middle half crop
            h, w = chan_im.shape
            h_start, h_end = h // 4, 3 * h // 4
            w_start, w_end = w // 4, 3 * w // 4
            mean_value = np.mean(chan_im[h_start:h_end, w_start:w_end])
            exposure_means.append(mean_value)
    
        # ---- Plot ----
        fig, ax = plt.subplots(figsize=(6, 4), dpi=96, facecolor="w", edgecolor="w")
    
        ax = sns.barplot(x=exposure_times, y=exposure_means)
        plt.xlabel("Exposure time")
        plt.ylabel("Mean Intensity")
        plt.title(f"Dark Counts for {cam} on {date}")
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    
        # Annotate bars
        for bar, value in zip(ax.patches, exposure_means):
            height = bar.get_height()
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                height + max(exposure_means) * 0.01,
                f"{value:.1f}",
                ha='center', va='bottom', fontsize=9, color='black'
            )
    
        folder = os.path.dirname(filepath)
        savepath = os.path.join(folder, f'darknoise_{cam}.png')
        plt.savefig(savepath)
        plt.show()
        plt.close()
    
        # ---- Save results to CSV ----
        # This allows one to plot the values over time, without having to import and process every ND2 file repeatedly
        if save_csv:
            csv_path = os.path.join(folder, "darknoise_values.csv")
            file_exists = os.path.isfile(csv_path)
        
            # Load existing rows to avoid duplicates
            existing_rows = set()
            if file_exists:
                with open(csv_path, "r", newline="") as f:
                    reader = csv.reader(f)
                    for row in reader:
                        existing_rows.add(tuple(row))
        
            with open(csv_path, "a", newline="") as csvfile:
                writer = csv.writer(csvfile)
        
                # Write header if new file
                if not file_exists:
                    header = ["Date", "Microscope", "File", "ExposureTime", "MeanIntensity"]
                    writer.writerow(header)
                    existing_rows.add(tuple(header))
        
                # Write data
                for exp_time, mean_val in zip(exposure_times, exposure_means):
                    row = (date, microscope, os.path.basename(filepath), exp_time, mean_val)
        
                    if row not in existing_rows:
                        writer.writerow(row)
                        existing_rows.add(row)
        
            print(f"Darknoise for {cam} saved")
        
        
    else:
        print(r"File type not recognized; this function only takes ND2s for now.")


def plot_camera_darknoise_overtime(parent_folder):
    """
    Collect all exposure_darknoise_values.csv files under parent_folder
    and plot dark noise evolution across acquisition dates.
    """

    root_path = Path(parent_folder)
    csv_files = list(root_path.rglob("*darknoise_values.csv"))
    print(f"Found {len(csv_files)} CSV files")

    # Read and combine all CSVs
    df_list = [pd.read_csv(file, sep=None, engine='python') for file in csv_files]
    df = pd.concat(df_list, ignore_index=True)

    # --- Clean up ---
    df["Date"] = pd.to_datetime(df["Date"], format="%y%m%d", errors='coerce')
    df["Camera"] = df["File"].str.extract(r"(cam\d+)")
    # Extract numeric value and unit (ms or s)
    exp = df["ExposureTime"].str.extract(r"(\d+)(ms|s)")
    exp[0] = exp[0].astype(float)

    # Convert everything to milliseconds
    df["Exposure_ms"] = exp.apply(lambda x: x[0] * 1000 if x[1] == "s" else x[0], axis=1)    
    df = df.sort_values("Date")

    # --- Plotting ---
    plt.figure(figsize=(12, 6))

    base_colors = {"cam4009": "blue", "cam4018": "red"}  # base color for each camera

    for camera, cam_group in df.groupby("Camera"):
        exposures = sorted(cam_group["Exposure_ms"].unique())
        n_exp = len(exposures)
        
        # create shades from light to dark based on exposure rank
        for i, exposure in enumerate(exposures):
            group = cam_group[cam_group["Exposure_ms"] == exposure]
            
            # Adjust darkness: darker for higher exposures
            factor = 0.4 + 0.6 * (i / (n_exp - 1 if n_exp > 1 else 1))  # range 0.4–1
            color = mcolors.to_rgba(base_colors[camera])
            color = (*color[:3], factor)  # adjust alpha to simulate darkness
                      
            plt.plot(group["Date"], group["MeanIntensity"],
                     marker='o',
                     linestyle='-',
                     color=color,
                     label=f"{camera} - {int(exposure)} ms")

    plt.xlabel("Date")
    plt.ylabel("Mean Dark Intensity")
    plt.title("Dark Noise Over Time")
    plt.xticks(rotation=45)
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    savepath = os.path.join(parent_folder, 'detector_darknoise_overtime.png')
    plt.savefig(savepath)
    plt.show()

def plot_gain_lin_leica(gain_lif_path=None, hyd_gains=None, pmt_gains=None, microscope=None, date=None):
    """
    Plot gain linearity curves for detectors in a LIF file.
    
    Parameters
    ----------
    gain_lif_path : str
        Path to the .lif file with gain series
    hyd_gains : list
        List of gain values (e.g. %) used for HyD detectors
    pmt_gains : list
        List of gain values (e.g. volts) used for PMT detectors
    microscope : str
        Microscope name
    date : str
        Acquisition date
    qc_path : str
        Output directory for saving the figure
    """
    results = []
    
    if gain_lif_path == None:
        gain_lif_path = askopenfilename()
        
    qc_path = os.path.dirname(gain_lif_path)

    gain_lif = LifFile(gain_lif_path)

    # Create a figure with subplots
    num_detectors = gain_lif.num_images
    fig, axes = plt.subplots(nrows=1, ncols=num_detectors, figsize=(4 * num_detectors, 4), dpi=80)

    if num_detectors == 1:
        axes = [axes]  # Ensure axes is iterable

    for detector_count, ax in enumerate(axes):
        detector_im = gain_lif.get_image(detector_count)
        detector_name = gain_lif.image_list[detector_count]['name']
        chan_means = []

        # Select gain series depending on detector type
        is_hyd = "HyD" in detector_name
        xs = hyd_gains if is_hyd else pmt_gains

        for chan_count in range(detector_im.channels):
            chan_im = np.asarray(detector_im.get_frame(c=chan_count))
            blur_chan_im = cv2.GaussianBlur(chan_im, (7, 7), 0)

            h, w = blur_chan_im.shape
            h_start, h_end = h // 4, 3 * h // 4
            w_start, w_end = w // 4, 3 * w // 4
            mean_val = np.mean(blur_chan_im[h_start:h_end, w_start:w_end])

            chan_means.append(mean_val)

        chan_means = np.array(chan_means)

        # Log transform PMTs only
        if not is_hyd:
            chan_means = np.log(chan_means)

        # Perform linear regression
        slope, intercept, r_value, p_value, std_err = linregress(xs, chan_means)
        r_squared = r_value ** 2
        results.append({
            "date": date,
            "microscope": microscope,
            "detector": detector_name,
            "slope": slope,
            "r_squared": r_squared
        })
        # Plot regression
        ax.set_title(f"{detector_name}")
        sns.regplot(x=xs, y=chan_means, ax=ax, scatter=True, ci=None, line_kws={"color": "red"})
        ax.text(0.05, 0.95, f"Slope = {slope:.3f}\n$R^2$ = {r_squared:.3f}",
            transform=ax.transAxes,
            verticalalignment='top', fontsize=10,
            bbox=dict(facecolor='white', alpha=0.6))
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_ylim(bottom=0)

    fig.suptitle(f"Gain linearity for {microscope} on {date}")
    fig.supxlabel("Gain (PMTs: V; HyDs: %)")
    fig.supylabel("Mean Intensity (HyD)\nlog(Mean Intensity) (PMT)", fontsize=14, x=0.07)
    #plt.tight_layout()
    savepath = os.path.join(qc_path, 'detector_gain_linearity.png')
    plt.savefig(savepath)
    plt.show()
    plt.close()
    
    results_df = pd.DataFrame(results)

    save_csv = os.path.join(qc_path, "detector_gain_linearity_slopes.csv")
    
    # Append if file exists so you can track over time
    if os.path.exists(save_csv):
        results_df.to_csv(save_csv, mode="a", header=False, index=False)
    else:
        results_df.to_csv(save_csv, index=False)
    
def plot_gain_lin_tif(gain_tif_path=None, detector_gains=None, microscope=None, date=None):
    """
    Plot gain linearity curves for detectors, taken from a multi-channel TIF, 
    where each channel has a different gain associated.
    
    Parameters
    ----------
    gain_if_path : str
        Path to the .tif file with gain series
    detector_gains : list
        List of gain values (e.g. % or volts) used for detectors, in same order of channels in the TIF
    microscope : str
        Microscope name
    date : str
        Acquisition date
    qc_path : str
        Output directory for saving the figure
    """
    results = []

    
    if gain_tif_path == None:
        gain_tif_path = askopenfilename()
        
    qc_path = os.path.dirname(gain_tif_path)

    gain_tif = tifffile.imread(gain_tif_path)

    # Create a figure with subplots
    num_detectors = 1
    fig, axes = plt.subplots(nrows=1, ncols=num_detectors, figsize=(4 * num_detectors, 4), dpi=80)

    if num_detectors == 1:
        axes = [axes]  # Ensure axes is iterable
    

    
    for detector_count, ax in enumerate(axes):
        detector_im = gain_tif
        chan_means = []


        for chan_count in range(detector_im.shape[0]):
            chan_im = np.asarray(detector_im[chan_count])
            blur_chan_im = cv2.GaussianBlur(chan_im, (7, 7), 0)

            # Middle half (instead of hardcoded 128:384 for 512x512)
            h, w = blur_chan_im.shape
            h_start, h_end = h // 4, 3 * h // 4
            w_start, w_end = w // 4, 3 * w // 4
            mean_val = np.mean(blur_chan_im[h_start:h_end, w_start:w_end])

            chan_means.append(mean_val)
        
        chan_means = np.array(chan_means)
        
        # Log transform PMTs (assumes all detectors from TIF images are PMTs)
        chan_means = np.log(chan_means)
        
        # Perform linear regression for R²
        slope, intercept, r_value, p_value, std_err = linregress(detector_gains, chan_means)
        r_squared = r_value ** 2
        results.append({
            "date": date,
            "microscope": microscope,
            "detector": detector_count,
            "slope": slope,
            "r_squared": r_squared
        })
        # Plot regression
       # ax.set_title(f"{detector_name}")
        sns.regplot(x=detector_gains, y=chan_means, ax=ax, scatter=True, ci=None, line_kws={"color": "red"})
        ax.text(0.05, 0.95, f"$R^2$ = {r_squared:.3f}", transform=ax.transAxes,
                verticalalignment='top', fontsize=10, bbox=dict(facecolor='white', alpha=0.6))
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_ylim(bottom=0)

    fig.suptitle(f"Gain linearity for {microscope} on {date}")
    fig.supxlabel("Gain (PMTs: V; HyDs: %)")
    fig.text(x=0, y=0.5, s="Mean Intensity", fontsize=14, rotation=90, va="center")
    plt.tight_layout()
    savepath = os.path.join(qc_path, 'detector_gain_linearity.png')
    plt.savefig(savepath)
    plt.show()
    plt.close()
    
    results_df = pd.DataFrame(results)

    save_csv = os.path.join(qc_path, "detector_gain_linearity_slopes.csv")
    
    # Append if file exists so you can track over time
    if os.path.exists(save_csv):
        results_df.to_csv(save_csv, mode="a", header=False, index=False)
    else:
        results_df.to_csv(save_csv, index=False)


def plot_gain_lin_overtime(qc_folder_path):
    """
    Collect all detector_gain_linearity_slopes.csv files under qc_folder_path
    and plot slopes for each detector across acquisition dates.
    """
    
    csv_files = glob.glob(os.path.join(qc_folder_path, "*", "detector_gain_linearity_slopes.csv"))


    if len(csv_files) == 0:
        print("No slope files found.")
        return None

    dfs = []

    for f in csv_files:
        df = pd.read_csv(f)
        dfs.append(df)

    all_data = pd.concat(dfs, ignore_index=True)

    # Convert date format
    all_data["date"] = pd.to_datetime(all_data["date"].astype(str), format="%y%m%d")

    microscopes = all_data["microscope"].unique()

    for scope in microscopes:

        scope_data = all_data[all_data["microscope"] == scope]

        g = sns.FacetGrid(
            scope_data,
            col="detector",
            col_wrap=3,
            sharey=False,
            height=4
        )

        g.map_dataframe(
            sns.lineplot,
            x="date",
            y="slope",
            marker="o"
        )

        # Set subplot titles to detector name only
        g.set_titles("{col_name}")

        g.set_axis_labels("Date", "Gain linearity slope")

        for ax in g.axes:
            ax.tick_params(axis='x', rotation=45)

        # Big microscope title
        g.fig.suptitle(scope, fontsize=16)

        plt.tight_layout()
        plt.subplots_adjust(top=0.88)

        plt.show()

    return all_data   

    

    
def plot_field_illum(filepath=None, microscope=None, date=None):
    """
    Generate and save field illumination uniformity plots from a multi-channel image.

    This function loads a field illumination image stack (e.g., from a flat-field
    fluorescence slide), applies heavy Gaussian smoothing to estimate the illumination
    profile, and overlays contour lines representing normalized intensity levels
    (50–95%). Each laser/channel is displayed in a separate subplot to visualize
    illumination uniformity across the field of view.

    Metadata such as objective, camera, and laser names are parsed from the filename
    if it follows the expected naming format (e.g., "field_illum_<objective>[_cam]_...").
    The resulting figure is saved in the same directory as the input file.

    Parameters
    ----------
    filepath : str, optional
        Path to the image file (.tif or .nd2). If None, a file selection dialog
        will open for the user to choose a file.

    microscope : str, optional
        Name of the microscope used for acquisition. If None, the user will be
        prompted to enter it.

    date : str, optional
        Date of acquisition in YYYY-MM-DD format. If None, the user will be
        prompted to enter it.

    Returns
    -------
    None
        Saves a figure showing field illumination uniformity with contour overlays
        for each laser/channel.

    Notes
    -----
    - Supports `.tif` files read with `tifffile` and `.nd2` files read with `ND2Reader`.
    - Assumes channels correspond to different laser lines and are stored along
      the first dimension of the image array.
    - Gaussian smoothing (kernel size 201×201) is used to approximate the global
      illumination profile.
    - Output filename includes the microscope name and objective (and camera if present).
    """

    if filepath == None:
        filepath = askopenfilename()
    if microscope == None:
        microscope = input("Enter name of microscope:  ")
    if date==None:
        date = input("Enter date of acquisition (YYYY-MM-DD):  ")
        
    dir_name = os.path.dirname(filepath)

    im_name = filepath.split(r"field_illum_")[-1].replace(".tif","")
    im_name_parts = im_name.split("_")
    
    obj = im_name_parts[0]
    if "_cam" in im_name:
        cam = im_name_parts[1]
        lasers = im_name_parts[2:]
    else:
        lasers = im_name_parts[1:]
        print(lasers)
    num_lasers = len(lasers)
    
    if ".tif" in filepath:
        images = tifffile.imread(filepath)
        
    elif ".nd2" in filepath:
        with ND2Reader(filepath) as images_test:
            images_test.bundle_axes = 'cyx'   # (channels, y, x)
            images = images_test[0] 
    else: 
        print("Unknown file type")
        
    fig = plt.figure(figsize=(8, 3))
    if "_cam" in im_name:
        plt.suptitle(f"Field illum at {microscope}, {obj}, {cam}, on {date}")
        savepath = os.path.join(dir_name, f'field_illum_{microscope}_{obj}_{cam}.png')
    else:
        plt.suptitle(f"Field illum at {microscope}, {obj}, on {date}")
        savepath = os.path.join(dir_name, f'field_illum_{microscope}_{obj}.png')

    
    for las_it, laser in enumerate(lasers):
        print("  Laser:", laser)
        # if las_it >= images.shape[0]:
        #     continue

        im = images[las_it]
        blur_im = cv2.GaussianBlur(im, (201, 201), 0)
        blur_im_norm = blur_im / blur_im.max()

        levels = [0.5, 0.6, 0.7, 0.8, 0.9, 0.95]
        ax = plt.subplot2grid((1, num_lasers), (0, las_it))
        ax.imshow(blur_im)
        
        contour = ax.contour(blur_im_norm, levels=levels, colors='white', alpha=0.5)
        ax.clabel(contour, inline=True, fontsize=10, fmt=lambda x: f"{int(x * 100)}%")

        ax.axis('off')
        ax.set_title(laser, fontsize=14)
        if las_it == 0:
            ax.set_ylabel(obj, fontsize=14)
            
    plt.tight_layout()
    plt.savefig(savepath)
    plt.close()    



    
    
    
    
    
    