# NOTE: this code is largely derived from the scripts and examples of geoCosiCorr3D.
# See its LICENSE file for more info.
# This script orthorectifies one SPOT scene using the specified geometric model and files.
import os.path

import warnings
from pathlib import Path

from geoCosiCorr3D.geoCore.constants import *
from geoCosiCorr3D.georoutines.file_cmd_routines import get_files_based_on_extension
from geoCosiCorr3D.geoCore.geoCosiCorrBaseCfg.BaseReadConfig import ConfigReader
from geoCosiCorr3D.geoOrthoResampling.geoOrtho import RSMOrtho
from geoCosiCorr3D.geoTiePoints.misc import parse_opt_report

from lib_parse_match_file import * # process() function for the binary matches file of ipmatch

import argparse

parser = argparse.ArgumentParser(description='Orthorectify a SPOT scene.')
parser.add_argument('sensor', metavar='sensor', type=str, nargs=1,
                    help='Which sensor are we going to process? Options: Spot1, Spot2, Spot3, Spot4, Spot5')
parser.add_argument('gsd', metavar='gsd', type=int, nargs=1,
                    help='Which ground resolution? Spot 1-4 have 10 (20) m in PAN (MS), Spot 5 has 5 (10) m in PAN (MS).')
parser.add_argument('raw_dir', metavar='raw_dir', type=str, nargs=1,
                    help='Path to the directory containing the raw SPOT scene and its metadata, excluding SCENE01/')
parser.add_argument('proc_subdir', metavar='proc_subdir', type=str, nargs=1,
                    help='Name of the directory to process within the processing dir')
parser.add_argument('dem_file', metavar='dem_file', type=str, nargs=1,
                   help='Path to the DEM file to be used for orthorectification')
parser.add_argument('ipfind_ipmatch_params', metavar = 'ipfind_ipmatch_params', type=str, nargs=1,
                    help='Selected ipfind and ipmatch parameter combination. For example 500k_50')
parser.add_argument('cosicorr_iter_best', metavar = 'cosicorr_iter_best', type=str, nargs=1,
                    help='Selected CosiCorr iteration (Python numbering, from 0). For example 3')
parser.add_argument('-o', dest='output_path', action='store', metavar = 'output_path', nargs=1, required=True,
                    help='Full path to the output file, including .tif extension')
args = parser.parse_args()


raw_img_dir = os.path.join(os.path.abspath(args.raw_dir[0]), "SCENE01")
ipfind_ipmatch_params = args.ipfind_ipmatch_params[0]
cosicorr_iter_best = args.cosicorr_iter_best[0]
dem_file = args.dem_file[0]
sensor = args.sensor[0]
ortho_gsd = int(args.gsd[0])
output_path = args.output_path[0]
proc_subdir = args.proc_subdir[0]

workspace_dir = os.path.abspath("./")


config_file = os.path.join(workspace_dir, "geo_ortho_config.yaml")
config = ConfigReader(config_file=config_file).get_config

raw_img_path = os.path.abspath(os.path.join(raw_img_dir, "IMAGERY.TIF"))

rsm_folder = os.path.join(proc_subdir, ipfind_ipmatch_params, "RSMs")
rsm_refinement_folder = os.path.join(proc_subdir, ipfind_ipmatch_params, "RSM_refinement")

dmp_file = os.path.abspath(os.path.join(rsm_folder, "METADATA.DIM"))



# BEGIN of function orthorectify
output_trans_path = "./transf_cur.tif"

config['ortho_params']['method']['metadata'] = dmp_file
config['ortho_params']['method']['sensor'] = sensor
config['ortho_params']['method']['corr_model'] = os.path.join(os.path.abspath(rsm_refinement_folder), "tps_proc_GCP_optloop_" + cosicorr_iter_best + "_correction.txt")
config['ortho_params']['GSD'] = ortho_gsd

RSMOrtho(input_l1a_path=raw_img_path,
        ortho_params=config['ortho_params'],
        output_ortho_path=output_path,
        output_trans_path=output_trans_path,
        dem_path=dem_file)
# END of function orthorectify

# Remove transformation raster.
Path(output_trans_path).unlink()
