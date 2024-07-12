## Code associated with the manuscript "Five decades of Abramov glacier dynamics reconstructed with multi-sensor optical remote sensing", under review (July 2024) for The Cryosphere.

- **Under box_method/:**               code to analyze glacier length changes from a vector file of terminus positions, using the rectilinear box method.

- **Under dem_spot_hrs/:**             pipeline to generate and postprocess DEMs from SPOT 5 HRS stereo pairs, including: stereo processing; filtering of saturated pixels; coregistration to a reference; and optional mosaicking of split scenes.

- **Under demdiff/:**                  pipeline to analyze a set of DEM differences, including: local coregistration; analysis and correction of DEM biases; local hypsometric filtering; aggregation of changes into polygons of interest; and uncertainty estimation.

- **Under spot_orthorectification/:**  pipeline to orthorectify a set of SPOT 1 to 5 scenes (HRV, HRVIR, HRG; both PAN and MS) from a reference orthophoto, including: GCP detection using ASP tools, automated selection of ipfind/ipmatch parameters, and distance-based GCP selection; RSM refinement with geoCosiCorr3D; orthorectification.

- **Under velocity_composites/:**      pipeline to derive annual composites of surface velocity from Sentinel-2, including: discovery of cloud-free scenes; selection of the scenes to correlate; image preprocessing; correlation; postprocessing of the displacements; extraction and plotting of longitudinal velocity profiles.

- **Under velocity_single_pairs/:**    pipeline to derive surface velocity from single pairs of orthoimages, including: image preprocessing; correlation; postprocessing of the displacements; extraction and plotting of longitudinal velocity profiles.


License notices of the third-party tools used in the manuscript are included in the respective subdirectories. All other code is licensed with GPL v3 as applicable.
