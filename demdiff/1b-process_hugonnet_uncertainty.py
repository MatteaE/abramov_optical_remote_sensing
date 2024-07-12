# This script performs the full Hugonnet uncertainty workflow
# (heteroscedasticity and error correlation).
# Inputs:
  # dh map
  # reference elevation map
  # vector file with polygons of regions over which the dh is averaged and the aggregated error is computed
  # indices of the polygons to select from the previous file (format: e.g. 1-3-6-7; field name is id)
  # polygons of unstable terrain mask
  # path to output directory (will be created and filled with standard-named files)
# Outputs:
  # dh error map computed from the slope and max curvature of the reference elevation map
  # plots detailing the dh error map calculation (dh error vs slope, aspect, curvature, etc.)
  # plots detailing the variogram calculations
  # map plots with the selected region polygon and its (average dh ± uncertainty) calculated according to several different variogram models
  # csv file with the average dh and uncertainties according to the different variogram models


print(". Importing modules...")
import warnings
warnings.filterwarnings("ignore", message=".*The 'nopython' keyword.*") # Filter out stupid numba warnings at startup.

import argparse
import geoutils as gu
import matplotlib.pyplot as plt
import numpy as np
np.seterr(divide='ignore', invalid='ignore') # Filter out numpy divide by zero warnings.

import os
import xdem
from xdem.spatialstats import nmad


# Command line arguments.
parser = argparse.ArgumentParser(description='Estimate spatially variable DEMdiff uncertainties depending on slope and max curvature, following Hugonnet et al. (2022)')
parser.add_argument('dh_path', metavar='dh_path', type=str, nargs=1,
                    help='Path to the DEMdiff grid')
parser.add_argument('ref_path', metavar='ref_path', type=str, nargs=1,
                    help='Path to the reference DEM from which we compute slope, aspect and curvatures')
parser.add_argument('poly_all_path', metavar='poly_all_path', type=str, nargs=1,
                    help='Path to a shapefile/geopackage with one or more polygons over which the averaged dh error is computed. Polygons to consider are selected over attribute id, using argument poly_sel_ids.')
parser.add_argument('poly_sel_ids', metavar='poly_sel_ids', type=str, nargs=1,
                    help='Indices of the polygons over which averaged dh error should be computed (individually for each polygon). They are expressed as e.g. 1-3-7-24.')
parser.add_argument('unstable_path', metavar='unstable_path', type=str, nargs=1,
                    help='Path to a shapefile/geopackage of unstable ground, to be used as mask for the estimations on stable terrain')
parser.add_argument('outdir_path', metavar='outdir_path', type=str, nargs=1,
                    help='Path to the base directory where everything is stored')
args = parser.parse_args()



# Load data, resample reference DEM to same grid as dh.
print(". Loading data...")
dh = xdem.DEM(args.dh_path[0])
ref_dem_raw = xdem.DEM(args.ref_path[0])
ref_dem = ref_dem_raw.reproject(dst_ref = dh, resampling = "bilinear")
unstable_outlines = gu.Vector(args.unstable_path[0])
mask_glacier = unstable_outlines.create_mask(dh)
poly_all_shp = gu.Vector(args.poly_all_path[0])
poly_ids = [int(x) for x in args.poly_sel_ids[0].split("-")] # We will use these indices when we select the areas on which to compute the averaged error.


outdir_path = args.outdir_path[0]
if not os.path.isdir(outdir_path):
    os.mkdir(outdir_path)

heterosc_dirpath = os.path.join(outdir_path, "heteroscedasticity_plots")
if not os.path.isdir(heterosc_dirpath):
    os.mkdir(heterosc_dirpath)

variogram_dirpath = os.path.join(outdir_path, "variograms_plots")
if not os.path.isdir(variogram_dirpath):
    os.mkdir(variogram_dirpath)

poly_plots_dirpath = os.path.join(outdir_path, "integrated_error_poly_plots")
if not os.path.isdir(poly_plots_dirpath):
    os.mkdir(poly_plots_dirpath)



# BEGIN of the heteroscedasticity modeling
# Compute topographic parameters from reference DEM.
print(". Computing topographic parameters...")
slope, aspect, planc, profc = xdem.terrain.get_terrain_attribute(
    dem=ref_dem, attribute=["slope", "aspect", "planform_curvature", "profile_curvature"]
)


# Convert to arrays, masking out unstable terrain.
dh_arr = dh[~mask_glacier].filled(np.nan)
slope_arr = slope[~mask_glacier].filled(np.nan)
aspect_arr = aspect[~mask_glacier].filled(np.nan)
planc_arr = planc[~mask_glacier].filled(np.nan)
profc_arr = profc[~mask_glacier].filled(np.nan)


# Compute binning of dh over topographic parameters.
print(". Binning 1/3...")
df = xdem.spatialstats.nd_binning(
    values=dh_arr,
    list_var=[slope_arr, aspect_arr, planc_arr, profc_arr],
    list_var_names=["slope", "aspect", "planc", "profc"],
    statistics=["count", xdem.spatialstats.nmad],
    list_var_bins=30,
)


# Plots: NMAD of dh versus slope, aspect, curvatures.
print(". Plots...")
xdem.spatialstats.plot_1d_binning(df, var_name="slope", statistic_name="nmad", label_var="Slope (degrees)", label_statistic="NMAD of dh (m)", out_fname = os.path.join(heterosc_dirpath, "dh_nmad_slope.png"))
xdem.spatialstats.plot_1d_binning(df, "aspect", "nmad", "Aspect (degrees)", "NMAD of dh (m)", out_fname = os.path.join(heterosc_dirpath, "dh_nmad_aspect.png"))

# Now recompute binning (for maximum curvature) using more sensible
# bin limits (uniform number of samples per bin, using quantiles).
print(". Binning 2/3...")
maxc_arr = np.maximum(np.abs(planc_arr), np.abs(profc_arr))
df = xdem.spatialstats.nd_binning(
    values=dh_arr,
    list_var=[maxc_arr],
    list_var_names=["maxc"],
    statistics=["count", np.nanmedian, xdem.spatialstats.nmad],
    list_var_bins=[np.nanquantile(maxc_arr, np.linspace(0, 1, 1000))],
)
print(". Plot...")
xdem.spatialstats.plot_1d_binning(df, "maxc", "nmad", "Maximum absolute curvature (100 m$^{-1}$)", "NMAD of dh (m)", out_fname = os.path.join(heterosc_dirpath, "dh_nmad_maxc.png"))


# Finally recompute binning for 2D plot of slope and maximum curvature.
print(". Binning 3/3...")
custom_bin_slope = np.unique(
    np.concatenate(
        [
            np.nanquantile(slope_arr, np.linspace(0, 0.95, 20)),
            np.nanquantile(slope_arr, np.linspace(0.96, 0.99, 5)),
            np.nanquantile(slope_arr, np.linspace(0.991, 1, 10)),
        ]
    )
)

custom_bin_curvature = np.unique(
    np.concatenate(
        [
            np.nanquantile(maxc_arr, np.linspace(0, 0.95, 20)),
            np.nanquantile(maxc_arr, np.linspace(0.96, 0.99, 5)),
            np.nanquantile(maxc_arr, np.linspace(0.991, 1, 10)),
        ]
    )
)

df = xdem.spatialstats.nd_binning(
    values=dh_arr,
    list_var=[slope_arr, maxc_arr],
    list_var_names=["slope", "maxc"],
    statistics=["count", np.nanmedian, xdem.spatialstats.nmad],
    list_var_bins=[custom_bin_slope, custom_bin_curvature],
)
print(". Plot (2D)...")
xdem.spatialstats.plot_2d_binning(
    df,
    "slope",
    "maxc",
    "nmad",
    "Slope (degrees)",
    "Maximum absolute curvature (100 m$^{-1}$)",
    "NMAD of dh (m)",
    scale_var_2="log",
    vmin=2,
    vmax=10,
    out_fname = os.path.join(heterosc_dirpath, "dh_nmad_slope_maxc.png")
)


# Now compute the numerical approximation for the binning results, and rescale
# it (two_step_standardization) to have exactly the same dispersion as the data.
print(". Error function...")
unscaled_dh_err_fun = xdem.spatialstats.interp_nd_binning(
    df, list_var_names=["slope", "maxc"], statistic="nmad", min_count=30
)
dh_err_stable = unscaled_dh_err_fun((slope_arr, maxc_arr))

print(
    ". The spread of elevation difference is {:.2f} "
    "compared to a mean predicted elevation error of {:.2f}.".format(
        xdem.spatialstats.nmad(dh_arr), np.nanmean(dh_err_stable)
    )
)

# zscores has the standardized dh scores, only on stable terrain, as 1D array.
zscores, dh_err_fun = xdem.spatialstats.two_step_standardization(
    dh_arr, list_var=[slope_arr, maxc_arr], unscaled_error_fun=unscaled_dh_err_fun
)


# # Print some examples of modeled dh errors.
# for s, c in [(0.0, 0.1), (10.0, 0.1), (20.0, 0.1), (30.0, 0.1), (40.0, 0.1), (50.0, 0.1), (60.0, 0.1), (70.0, 0.1),
#              (0.0, 1.0), (10.0, 1.0), (20.0, 1.0), (30.0, 1.0), (40.0, 1.0), (50.0, 1.0), (60.0, 1.0), (70.0, 1.0),
#              (0.0, 5.0), (10.0, 5.0), (20.0, 5.0), (30.0, 5.0), (40.0, 5.0), (50.0, 5.0), (60.0, 5.0), (70.0, 5.0),
#              (0.0, 20.0), (10.0, 20.0), (20.0, 20.0), (30.0, 20.0), (40.0, 20.0), (50.0, 20.0), (60.0, 20.0), (70.0, 20.0)]:
#     print("Elevation measurement error for slope of {:.0f} degrees, "
#         "curvature of {:.2f} m-1: {:.1f}".format(s, c / 100, dh_err_fun((s, c))) + " meters.")


# Compute and save the elevation error map.
print(". Error map...")
maxc = np.maximum(np.abs(profc), np.abs(planc))
dh_err = dh.copy(new_array=dh_err_fun((slope.data, maxc.data)))
dh_err.save(os.path.join(outdir_path, "dh_error.tif"))


# Save topographic parameters (the two predictors of dh error) for inspection.
slope.save(os.path.join(heterosc_dirpath, "slope.tif"))
maxc.save(os.path.join(heterosc_dirpath, "maxc.tif"))


plt.close("all") # Close all pyplot windows (release memory from figures plotting).

# END of the heteroscedasticity modeling


# BEGIN of the variogram modeling

# Compute standardized dh, over the full grid.
z_dh = dh.data / dh_err

# Filter out unstable terrain and large outliers, print statistics along the way.
print(f". STD -- NMAD, before unstable terrain masking: {np.nanstd(z_dh.data*dh_err):.2f} -- {xdem.spatialstats.nmad(z_dh.data*dh_err):.2f} m.")
z_dh.data[mask_glacier.data] = np.nan
print(f". STD -- NMAD, after masking and before 4-NMAD filtering: {np.nanstd(z_dh.data*dh_err):.2f} -- {xdem.spatialstats.nmad(z_dh.data*dh_err):.2f} m.")
z_dh.data[np.abs(z_dh.data) > 4] = np.nan
print(f". STD -- NMAD, after 4-NMAD filtering: {np.nanstd(z_dh.data*dh_err):.2f} -- {xdem.spatialstats.nmad(z_dh.data*dh_err):.2f} m.")


# Scale-correction for the standardization, to ensure that the standard deviation of the data is exactly 1.
# NOTE: the original tutorial here uses nmad() to compute the "standard deviation",
# so that the scale factor is actually the NMAD and the empirical variogram later converges to 1.48 instead of 1.
# If instead we later use "dowd" estimator in the empirical variogram sampling, it converges to 2.
# We have discussed with Romain Hugonnet (issue on GitHub) and he confirmed that we are still good to go.
print(f". Standard deviation before scale-correction: {nmad(z_dh.data):.2f}")
scale_fac_std = nmad(z_dh.data)
z_dh = z_dh / scale_fac_std
print(f". Standard deviation after scale-correction: {nmad(z_dh.data):.2f}")


# Save and plot standardized elevation differences
z_dh_to_save = dh.copy(new_array=z_dh)
z_dh_to_save.save(os.path.join(outdir_path, "dh_error_standardized.tif"))

plt.figure(figsize=(8, 5))
plt_extent = [
    ref_dem.bounds.left,
    ref_dem.bounds.right,
    ref_dem.bounds.bottom,
    ref_dem.bounds.top,
]
ax = plt.gca()
unstable_outlines.ds.plot(ax=ax, fc="none", ec="tab:gray")
ax.plot([], [], color="tab:gray", label="Unstable outlines")
plt.imshow(z_dh.squeeze(), cmap="RdYlBu", vmin=-3, vmax=3, extent=plt_extent)
cbar = plt.colorbar()
cbar.set_label("Standardized elevation differences")
plt.legend(loc="lower right")
plt.savefig(os.path.join(outdir_path, "dh_map_standardized.png"))
plt.close()



# Sample empirical variogram and plot it. We use estimator="dowd" as instructed by Romain Hugonnet.
print(". Sampling empirical variogram...")
df_vgm = xdem.spatialstats.sample_empirical_variogram(
    values=z_dh.data.squeeze(), estimator = "dowd", gsd=dh.res[0], subsample=300, n_variograms=10, random_state=42
)
print(". Plotting empirical variogram...")
xdem.spatialstats.plot_variogram(df_vgm, out_fname = os.path.join(variogram_dirpath, "variogram_empirical_linscale.png"))
xdem.spatialstats.plot_variogram(df_vgm, xscale="log", out_fname = os.path.join(variogram_dirpath, "variogram_empirical_logscale.png"))
xdem.spatialstats.plot_variogram(df_vgm, xscale_range_split=[100, 1000, 10000], out_fname = os.path.join(variogram_dirpath, "variogram_empirical_splitscale.png"))


# Try to fit various functions (both single and double range) to the empirical variogram.
# Plot results.
print(". Fitting and plotting variogram model...")

maxfev = 10000
name_vgm = "Double-range exponential+exponential model"


func_sum_vgm, params_vgm = xdem.spatialstats.fit_sum_model_variogram(
    # list_models=["Gaussian", "Spherical"], empirical_variogram=df_vgm, maxfev=maxfev
    list_models=["Exponential", "Exponential"], empirical_variogram=df_vgm, maxfev=maxfev
)

xdem.spatialstats.plot_variogram(
    df_vgm,
    list_fit_fun=[func_sum_vgm],
    list_fit_fun_label=[name_vgm],
    xscale_range_split=[100, 1000, 10000],
    out_fname = os.path.join(variogram_dirpath, "variogram_modeled.png")
)

plt.close("all")




# Now integrate the variogram model over various surface areas,
# to see how it compares to a Monte-Carlo sampling (patches method).
# We also introduce a component of fully correlated variance (0 to 10 %),
# automatically estimating which fraction leads to the best match with
# the Monte-Carlo method.
# NOTE: compared to the tutorial scripts, we do all this
# using the standardized dh values.

fullycorr_pct_max   = 10
fullycorr_pct_n     = 21
fullycorr_pcts_list = np.linspace(0, fullycorr_pct_max, fullycorr_pct_n)

# Precalculate fully correlated error, it does not depend on the integration area.
stderr_fullycorr_list = []
fullycorr_pct_labels_list = []

for fullycorr_pct_id in range(fullycorr_pct_n):
    fullycorr_pct_cur = fullycorr_pcts_list[fullycorr_pct_id]
    stderr_fullycorr_list.append(np.sqrt(fullycorr_pct_cur*0.01))     # 0.01: percent to fraction.
    fullycorr_pct_labels_list.append(f"vgm + {fullycorr_pct_cur} %")  # Labels, used in plot.


# List with one list of errors for each percent of fully correlated variance, from 0 to 20.
# Each sub-list has one item per sampled area.
stderr_tot_list = [[] for fullycorr_pct_cur in fullycorr_pcts_list]

# Now integrate variogram over the areas and combine with fully correlated error.
print(". Integrating variogram over various surface areas...")
# areas = np.linspace(20, 10000, 50) ** 2
# areas = np.linspace(100, 4000, 20) ** 2 # Smaller max area, since the polygons we consider are not huge.
areas = [4000 * 2 ** (i) for i in range(2,13,1)]
for area in areas:

    # print(f"Area: {area:.0f} m²")
    # Number of effective samples integrated over the area.
    # This is the contribution from the variogram.
    neff_vgm = xdem.spatialstats.number_effective_samples(area, params_vgm)

    for fullycorr_pct_id in range(fullycorr_pct_n):

        fullycorr_pct_cur = fullycorr_pcts_list[fullycorr_pct_id]
        stderr_vgm        = np.sqrt(1-(fullycorr_pct_cur*0.01)) / np.sqrt(neff_vgm)
        stderr_tot        = stderr_fullycorr_list[fullycorr_pct_id] + stderr_vgm

        stderr_tot_list[fullycorr_pct_id].append(stderr_tot)





# Compute Monte-Carlo sampling of the integrated error, to be
# plotted alongside the error estimated with the variogram models.
# We do this using the standardized error, to be comparable with
# the variogram results.
# areas_emp = [4000 * 2 ** (i) for i in range(10)]
df_patches = xdem.spatialstats.patches_method(z_dh, gsd=dh.res[0], areas=areas)




# Plot standardized uncertainty by averaging area and fully correlated variance.
fig, ax = plt.subplots()

for fullycorr_pct_id in [0,1,2,5,10,20]:
    plt.plot(np.asarray(areas) / 1000000, stderr_tot_list[fullycorr_pct_id], label=fullycorr_pct_labels_list[fullycorr_pct_id])

#plt.plot(np.asarray(areas) / 1000000, list_stderr_vgm, label=name_vgm)
plt.scatter(
    df_patches.exact_areas.values / 1000000,
    df_patches.nmad.values,
    label="Empirical estimate",
    color="black",
    marker="x",
)
plt.xlabel("Averaging area (km²)")
plt.ylabel("Uncertainty in the mean standardized elevation difference [-]")
plt.xscale("log")
plt.yscale("log")
plt.legend(loc="lower left")
plt.savefig(os.path.join(variogram_dirpath, "std_uncertainty_by_area.png"))


plt.close("all")


# Finally select the most realistic amount of fully correlated variance based on the comparison of the patches method against the integrated variogram.
# We use RMS as metric.
fullycorr_rms_list = []
for fullycorr_pct_id in range(fullycorr_pct_n):
    fullycorr_rms_list.append(np.sqrt(np.mean(np.subtract(stderr_tot_list[fullycorr_pct_id], df_patches.nmad.values)**2)))
fullycorr_pct_sel_id = np.argmin(fullycorr_rms_list)
fullycorr_pct_sel    = fullycorr_pcts_list[fullycorr_pct_sel_id]

print(f". Best estimate for fraction of fully correlated variance: {fullycorr_pct_sel} %.")
# END variogram modeling.

# END part run just once.


# BEGIN loop on the selected polygons.

# Here we estimate the final integrated uncertainty for each selected polygon, using the exponential+exponential variogram model and the selected best estimate for the fully correlated variance.

for poly_cur_id in poly_ids:

    poly_cur = poly_all_shp[poly_all_shp.ds.id == poly_cur_id]
    poly_mask = poly_cur.create_mask(dh)

    # If the Numpy array is fully masked, it means we have no valid data within
    # our polygon, thus we skip the rest of processing for the current polygon.
    if (dh.data[poly_mask.data]).mask.all():
        print(". . Polygon " + str(poly_cur_id) + " has no valid dh values!")
    else:

        # Compute mean elevation change over the current polygon, it will be be added to the plots.
        poly_dh = np.nanmean(dh.data[poly_mask.data])

        poly_cur_plots_dirpath = os.path.join(poly_plots_dirpath, str(poly_cur_id))
        if not os.path.isdir(poly_cur_plots_dirpath):
            os.mkdir(poly_cur_plots_dirpath)

        f_out = open(os.path.join(poly_cur_plots_dirpath, "dh_error_integrated.csv"), "w")
        f_out.write("poly_dh_m,vgm_model,poly_neff,poly_dh_err_m,fullycorr_pct\n")


        # Compute number of effective samples for the polygon according to the current variogram model.
        # NOTE: we have verified that neff_circular_approx_numerical() gives the same result as number_effective_samples().
        poly_neff = xdem.spatialstats.neff_circular_approx_numerical(
            area=poly_cur.ds.area.values[0], params_variogram_model=params_vgm
        )
        # print(f"Number of effective samples of selected polygon: {poly_neff:.1f}")
        poly_z_err_vgm       = np.sqrt(1-(fullycorr_pct_sel*0.01)) / np.sqrt(poly_neff)
        poly_z_err_fullycorr = np.sqrt(fullycorr_pct_sel*0.01)
        poly_z_err_tot       = poly_z_err_fullycorr + poly_z_err_vgm
        # print(f"Standardized integrated error of selected polygon: {poly_z_err_tot:.2f}")

        # Destandardize the spatially integrated uncertainty based on the measurement error dependent on slope and
        # maximum curvature. This yields the uncertainty into the mean elevation change for the selected polygon.
        fac_poly_dh_err = scale_fac_std * np.nanmean(dh_err[poly_mask.data])
        poly_dh_err = fac_poly_dh_err * poly_z_err_tot

        f_out.write(f"{poly_dh:.5f},{name_vgm},{poly_neff:.1f},{poly_dh_err:.5f},{fullycorr_pct_sel:.1f}\n",)

        # Plot the result.
        # print("Plotting final figure...")
        plt.figure(figsize=(8, 5))
        ax = plt.gca()
        plt.imshow(dh.data, cmap="RdYlBu", vmin=-10, vmax=10, extent=plt_extent)
        cbar = plt.colorbar(ax=ax)
        cbar.set_label("Elevation differences (m)")
        poly_cur.ds.plot(ax=ax, fc="none", ec="tab:gray", lw=2)
        plt.plot([], [], color="tab:gray", label="Selected polygon")
        ax.text(
            poly_cur.ds.centroid.x.values[0],
            poly_cur.ds.centroid.y.values[0] - 1500,
            f"{poly_dh:.2f} \n$\\pm$ {poly_dh_err:.2f}",
            color="tab:gray",
            fontweight="bold",
            va="top",
            ha="center",
            fontsize=12,
        )
        plt.legend(loc="lower left")
        plt.savefig(os.path.join(poly_cur_plots_dirpath, "integrated_dh_err.png"))
        plt.close()

        f_out.close()
