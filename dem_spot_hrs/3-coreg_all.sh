# This script xdem-coregisters SPOT 5 HRS DEMs to NASADEM reference, with RGI6.0 unstable mask.
target_dir="../dem_v2/"
output_dir="../dem_v3/"

ref_dem="/PATH/TO/REFERENCE/DEM/FOR/COREGISTRATION/"
unstable_gpkg="/PATH/TO/UNSTABLE/OUTLINES/VECTOR/FILE/"


demlist=`ls $target_dir`

for demfile in $demlist
do
	
	echo "Processing $demfile..."
	python coreg_with_crop.py $ref_dem $target_dir/$demfile $unstable_gpkg $output_dir/"${demfile%.*}_coreg.tif"
done
