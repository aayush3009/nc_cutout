#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 10:51:52 2024

@author: aayushsaxena
"""

### JWST NIRCam cutout creator

import numpy as np
import os
from pathlib import Path

import glob

from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u

import montage_wrapper as montage

import pandas as pd

### Parallel
from joblib import Parallel, delayed
import multiprocessing


### Define function that iterates over images
def parallel_cutouts(imdata, filtname, src_ID, src_RA, src_Dec, wcs):
    """

    Parameters
    ----------
    im : FITS image
        2D FITS file
    src_ID : list
        List of source IDs.
    src_RA : list
        list of source RAs.
    src_Dec : list
        list of source Declinations.

    Returns
    -------
    None.

    """
    
    ### Define the cutout size in arcsec
    cutout_size = 2. * 0.0002777777777777778  ### Multiply with factor to convert arcseconds to degrees

    print(f"Creating cutouts for source {src_ID} for filter {filtname}...")
        
    ### Try making the directory to store gallery images
    try:
        os.makedirs(f"{filtname}_cutouts", exist_ok=True)
    except OSError as error:
        print(error)
 
    out_image = f"./{filtname}_cutouts/{src_ID}_cutout_{filtname}.fits"
    
    ### Run montage
    # montage.mSubimage(im, out_image, src_RA, src_Dec, cutout_size, hdu=1)

    # Check if cutout already exists. If it does, skip cutout creation
    if Path(out_image).exists():
        print("Cutout exists, skipping...")
    else:    
        position = SkyCoord(ra=src_RA, dec=src_Dec, unit="deg")
        size = cutout_size * u.deg
        
        # Cutout2D only reads the necessary pixels from disk
        try:
            cutout = Cutout2D(imdata, position, size, wcs=wcs)
            # Save the cutout
            hdu = fits.PrimaryHDU(data=cutout.data, header=cutout.wcs.to_header())
            hdu.writeto(out_image, overwrite=True)
        except ValueError:
            print("ValueError encountered, skipping cutout creation...")
    
    return None

def main():
    ### Source list
    src_list = pd.read_csv("/Users/aayushsaxena/Desktop/Oxford/JADES/beta_slopes/v5.1_measurements/beta_slopes_v5.1_NC_photometry.csv")
    
    ### Subset based on GOODS-S or GOODS-N
    src_list_gs = src_list[(src_list["SURVEY"]=="goods-s-deephst") |  (src_list["SURVEY"]=="goods-s-deepjwst") | (src_list["SURVEY"]=="goods-s-ultradeep") |
                           (src_list["SURVEY"]=="goods-s-mediumhst") | (src_list["SURVEY"]=="goods-s-mediumjwst")]
    src_list_gn = src_list[(src_list["SURVEY"]=="goods-n-mediumhst") |  (src_list["SURVEY"]=="goods-n-mediumjwst")]
    
    ### Convert pandas subset columns into numpy
    gs_src_ID = src_list_gs["ID_1"].to_numpy()
    gs_src_RA = src_list_gs["RA"].to_numpy()
    gs_src_Dec = src_list_gs["DEC"].to_numpy()
    
    gn_src_ID = src_list_gn["ID_1"].to_numpy()
    gn_src_RA = src_list_gn["RA"].to_numpy()
    gn_src_Dec = src_list_gn["DEC"].to_numpy()
    
    ### Image directories
    im_dir_gs = "/Volumes/data5tb/JADES-NIRCam/GOODS-S-NC-v0.9/"
    im_dir_gn = "/Volumes/data5tb/JADES-NIRCam/GOODS-N-NC-v0.7/"

    ### Num of cores to be used
    num_cores = multiprocessing.cpu_count() - 2
    
    ### Create cutouts for all sources in GOODS_S
    os.chdir(im_dir_gs)
    gs_imlist = glob.glob("*.fits")
    for gs_im in gs_imlist:
        # Extract the filter name
        filtname = gs_im.split('_')[1]
        print(f"Creating cutouts from the image: {gs_im}")

        with fits.open(gs_im, memmap=True) as hdul:
            # Extract the image data
            imdata = hdul['SCI'].data
            # Extract the wcs data
            wcs = WCS(hdul['SCI'].header)
            # Start the cutout creation
            print("Starting cutout creation...")

            # Not parallelised, go through the source list one at a time. Most stable!
            for i in range(len(gs_src_ID)):
                results_gs = parallel_cutouts(imdata=imdata, filtname=filtname, src_ID=gs_src_ID[i], src_RA=gs_src_RA[i], src_Dec=gs_src_Dec[i], wcs=wcs)

            # Parallelized, feeds in multiple sources to multiple threads, doesn't seem to work for some reason...
            # results_gs = Parallel(n_jobs=num_cores)(delayed(parallel_cutouts)(imdata=imdata, filtname=filtname, src_ID=gs_src_ID[i], src_RA=gs_src_RA[i], src_Dec=gs_src_Dec[i], wcs=wcs) for i in range(len(gs_src_ID)))
            
            # Close the fits file to free up memory
            hdul.close()
            

    ### Create cutouts for all sources in GOODS_N
    os.chdir(im_dir_gn)
    gn_imlist = glob.glob("*.fits")
    for gn_im in gn_imlist:
    # Extract the filter name
        filtname = gn_im.split('_')[1]
        print(f"Creating cutouts from the image: {gn_im}")

        with fits.open(gn_im, memmap=True) as hdul:
            # Extract the image data
            imdata = hdul['SCI'].data
            # Extract the wcs data
            wcs = WCS(hdul['SCI'].header)
            # Start the cutout creation
            print("Starting cutout creation...")

            # Not parallelised, go through the source list one at a time. Most stable!
            for i in range(len(gn_src_ID)):
                results_gn = parallel_cutouts(imdata=imdata, filtname=filtname, src_ID=gn_src_ID[i], src_RA=gn_src_RA[i], src_Dec=gn_src_Dec[i], wcs=wcs)

            # Parallelized, feeds in multiple sources to multiple threads, doesn't seem to work for some reason...
            # results_gs = Parallel(n_jobs=num_cores)(delayed(parallel_cutouts)(imdata=imdata, filtname=filtname, src_ID=gs_src_ID[i], src_RA=gs_src_RA[i], src_Dec=gs_src_Dec[i], wcs=wcs) for i in range(len(gs_src_ID)))
            
            # Close the fits file to free up memory
            hdul.close()

    print("All done!")
    
    
if __name__ == '__main__':
    main()
    
    
    