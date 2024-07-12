#!/usr/bin/env python

"""
# Author : Saif Aati
# Contact: SAIF AATI  <saif@caltech.edu> <saifaati@gmail.com>
# Copyright (C) 2022
"""

# Reworked in 2024 by Enrico Mattea to correlate a given image with
# a set of images instead of a single one (faster than calling Python
# every single time).

import sys
sys.path.insert(0, "/PATH/TO/Geospatial-COSICorr3D")

from geoCosiCorr3D.geoImageCorrelation.correlate import Correlate
from geoCosiCorr3D.geoCore.constants import *
from geoCosiCorr3D.georoutines.geo_utils import cRasterInfo
import argparse


parser = argparse.ArgumentParser(description='Call Cosi-CORR correlation on two orthoimages.')
parser.add_argument('image1', metavar='image1', type=str, nargs=1,
                    help='Path to the first image to correlate')
parser.add_argument('image2_listpath', metavar='image2_listpath', type=str, nargs=1,
                    help='Path to a text file with the paths of all the images to correlate with the first image')
parser.add_argument('outdir', metavar='outdir', type=str, nargs=1,
                    help='Path to the directory where the output correlation file will be stored')
args = parser.parse_args()


img1 = args.image1[0]
img2_listpath = args.image2_listpath[0]
outdir = args.outdir[0]


img2_listfile = open(img2_listpath, "r")
img2_list = img2_listfile.readlines()
img2_n = len(img2_list)


# windowSize:   [x_start, y_start, x_end, y_end], integer, in px. Power of two. Sizes of the sliding window. If 4 equal numbers, sliding window has constant size. Else, it is refined into a smaller one.
# steps:        [x, y], integer, in px. X and Y offsets to move the sliding window.
# maskTh:       numeric. Used to mask frequencies. Value close to 1 recommended.
# resampling:   True/False. Patches to correlate are relocated from sinc resampling. Removes sub-pixel biases, but increases processing time.
# nbIterations:  integer. Number of times per measurement the frequency mask should be adaptively re-computed. Values 2-4 recommended.
# grid:         True/False.  Force top-left corner coordinates to be integer multiple of the ground resolution. Useful for mosaics / multiple processing.
test_corr_config = {"correlator_name": "frequency",
               "correlator_params": {
                   "window_size": [64,64,32,32],
                     "step": 2 * [4],
                     "mask_th": 0.9,
                     "resampling": False,
                     "nb_iters": 4,
                     "grid": False
                     }
               }


for img_id in range(img2_n):

    corr_obj = Correlate(base_image_path=img1,
                        target_image_path=img2_list[img_id].rstrip(), # Remove trailing whitespace from path, if any.
                        output_corr_path=outdir,
                        corr_config=test_corr_config,
                        corr_show=False)


img2_listfile.close()
