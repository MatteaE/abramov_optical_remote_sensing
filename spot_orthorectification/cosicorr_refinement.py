# NOTE: this code is largely derived from the scripts and examples of geoCosiCorr3D.
# See its LICENSE file for more info.
# This script runs CosiCorr RSM refinement from
# the given control points file (and image files).
# It is called several times by 1-compute_refine_rsm.R
# to find the parameter combination yielding the
# best orthophoto.

from pathlib import Path
import shutil

from geoCosiCorr3D.geoCore.constants import *
from geoCosiCorr3D.geoCore.geoCosiCorrBaseCfg.BaseReadConfig import ConfigReader
from geoCosiCorr3D.geoOptimization.gcpOptimization import cGCPOptimization
from geoCosiCorr3D.georoutines.file_cmd_routines import get_files_based_on_extension
from geoCosiCorr3D.geoCore.core_RSM import RSM
from geoCosiCorr3D.geoRSM.geoRSM_generation import geoRSM_generation
from geoCosiCorr3D.geoTiePoints.Tp2GCPs import TPsTOGCPS

import argparse


# BEGIN parse args
parser = argparse.ArgumentParser(description='Generate and refine an RSM with CosiCorr.')
parser.add_argument('sensor', metavar='sensor', type=str, nargs=1,
                    help='Which sensor are we going to process? Options: Spot1, Spot2, Spot3, Spot4, Spot5')
parser.add_argument('gsd', metavar='gsd', type=int, nargs=1,
                    help='Which ground resolution? Spot 1-4 have 10 (20) m in PAN (MS), Spot 5 has 5 (10) m in PAN (MS).')
parser.add_argument('proc_dir', metavar='proc_dir', type=str, nargs=1,
                    help='Path to the processing directory where RSM refinement will be stored')
parser.add_argument('raw_img_dir', metavar='raw_img_dir', type=str, nargs=1,
                    help='Path to the directory containing the raw SPOT scene and its metadata, up to and including SCENE01/ (the raw scene inside must be called IMAGERY.TIF)')
parser.add_argument('ref_ortho', metavar='ref_ortho', type=str, nargs=1,
                    help='Path to the file of the reference orthorectified scene')
parser.add_argument('dem_file', metavar='dem_file', type=str, nargs=1,
                    help='Path to the DEM file to be used for orthorectification')
parser.add_argument('tps_path', metavar='tps_path', type=str, nargs=1,
                    help='Path to the plain-text file with the processed tie points.')
args = parser.parse_args()


sensor      = args.sensor[0]
ortho_gsd   = int(args.gsd[0])
proc_dir    = os.path.abspath(args.proc_dir[0])
raw_img_dir = os.path.abspath(args.raw_img_dir[0])
dem_file    = os.path.abspath(args.dem_file[0])
ref_ortho   = os.path.abspath(args.ref_ortho[0])
tps_path    = args.tps_path[0]
# END parse args


# BEGIN configuration
workspace_dir = os.path.abspath("./")

# Set sensor-specific constants.
sensor_rsm = {"Spot1": SENSOR.SPOT1,
              "Spot2": SENSOR.SPOT2,
              "Spot3": SENSOR.SPOT3,
              "Spot4": SENSOR.SPOT4,
              "Spot5": SENSOR.SPOT5}[sensor]

# Load config file.
config_file = os.path.join(workspace_dir, "geo_ortho_config.yaml")
config = ConfigReader(config_file=config_file).get_config

# Path to the raw image.
raw_img_path = os.path.join(raw_img_dir, "IMAGERY.TIF")

# Simple adaptive correlator window for GCPs optimization.
window_size = [64, 64, 64, 64]
if ortho_gsd == 20:
    window_size = [32, 32, 32, 32]
config['opt_corr_config']['correlator_params']['window_size'] = window_size
# END configuration


# BEGIN processing
dmp_file = get_files_based_on_extension(raw_img_dir, "*.DIM")[0]
os.chdir(workspace_dir)

rsm_folder = os.path.join(proc_dir, "RSMs")
os.mkdir(rsm_folder)

shutil.copy(dmp_file, os.path.join(rsm_folder, os.path.basename(dmp_file)))
dmp_file = os.path.join(rsm_folder, os.path.basename(dmp_file))

# Generate RSM model.
rsm_model = geoRSM_generation(sensor_rsm, metadata_file=dmp_file, debug=True)

# Compute footprint.
_, _, _, gdf_fp = RSM.compute_rsm_footprint(rsm_model=rsm_model, dem_file=dem_file)

# Convert TPs into GCPs.
gcp = TPsTOGCPS(in_tp_file=tps_path,
                base_img_path=raw_img_path,
                ref_img_path=ref_ortho,
                dem_path=dem_file,
                debug=True)

sat_model_params = {'sat_model': SATELLITE_MODELS.RSM, 'metadata': dmp_file, 'sensor': sensor}
rsm_refinement_dir = os.path.join(proc_dir, "RSM_refinement")
Path(rsm_refinement_dir).mkdir(parents=True, exist_ok=True)

# NOTE: this can fail if the DEM contains gaps.
opt = cGCPOptimization(gcp_file_path=gcp.output_gcp_path,
                    raw_img_path=raw_img_path,
                    ref_ortho_path=ref_ortho,
                    sat_model_params=sat_model_params,
                    dem_path=dem_file,
                    opt_params=config['opt_params'],
                    opt_gcp_file_path=os.path.join(rsm_refinement_dir, Path(gcp.output_gcp_path).stem + "_opt.pts"),
                    corr_config=config['opt_corr_config'],
                    debug=True,
                    svg_patches=True)
# END of function rsm_refinement
