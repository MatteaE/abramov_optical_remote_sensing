#!/usr/bin/env python
# This script does DEM coregistration using xDEM, with unstable terrain mask.
# NOTE: this version does not include cropping of the reference DEM,
# we assume it is already small enough to have a fast processing.

import warnings
warnings.filterwarnings("ignore", message=".*The 'nopython' keyword.*") # Filter out numba warnings at startup.

import geoutils as gu
import numpy as np
import argparse

import xdem
from xdem import coreg

# print("Performing coregistration...")

parser = argparse.ArgumentParser(description='Coregister two DEMs using xdem.')
parser.add_argument('ref_path', metavar='ref_path', type=str, nargs=1,
                    help='Path to the reference DEM')
parser.add_argument('tba_path', metavar='tba_path', type=str, nargs=1,
                    help='Path to the DEM that should be coregistered')
parser.add_argument('unstable_path', metavar='unstable_path', type=str, nargs=1,
                    help='Path to a shapefile of unstable ground')
parser.add_argument('out_path', metavar='out_path', type=str, nargs=1,
                    help='Output file path')
args = parser.parse_args()

# Load the data using xdem and geoutils (could be with rasterio and geopandas instead)
# print("Loading input data...")

# Load the reference DEM
ref_dem = xdem.DEM(args.ref_path[0])
# Load the DEM that we have to align
tba_dem = xdem.DEM(args.tba_path[0])
# Load glacier outlines from RGI 6.0. This will act as the unstable ground.
unstable_outlines = gu.Vector(args.unstable_path[0])


# This is a boolean numpy 2D array. Note the bitwise not (~) symbol
# print("Creating inlier mask...")
inlier_mask = ~unstable_outlines.create_mask(ref_dem)

# Coregistration pipeline as suggested by https://xdem.readthedocs.io/en/latest/coregistration.html, plus subtraction of residual median bias.
pipeline = coreg.NuthKaab() + coreg.Tilt() + coreg.VerticalShift(vshift_func=np.nanmedian)

# Fit pipeline to data
# print("Fitting coregistration pipeline...")
pipeline.fit(ref_dem, tba_dem, inlier_mask=inlier_mask)


# Apply the transformation matrix to the data (or any other data)
# print("Applying coregistration pipeline...")
aligned_dem = pipeline.apply(tba_dem)

# Write aligned DEM.
# print("Writing output...")
aligned_dem.save(args.out_path[0])


# print("All done.")
