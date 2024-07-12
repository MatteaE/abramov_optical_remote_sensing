# This definition is only sourced by 1-prepreprocess_set.R
# Code derived from GIV (Van Wyk de Vries and Wickert, 2021).
matlab_filter_script_code_clahe <- c(
  "[A,R]=geotiffread(fullfile(datafolder,fn));",
  "% CLAHE filter",
  "in = rescale(A); %Rescale to 0-1, necessary for CLAHE",
  "numberoftiles1=round(size(in,1)/10);",
  "numberoftiles2=round(size(in,2)/10);",
  "if numberoftiles1 < 2",
  "numberoftiles1=2;",
  "end",
  "if numberoftiles2 < 2",
  "numberoftiles2=2;",
  "end",
  "in=adapthisteq(in, 'NumTiles',[numberoftiles1 numberoftiles2], 'ClipLimit', 0.01, 'NBins', 256, 'Range', 'full', 'Distribution', 'uniform');",
  "geotiffwrite(fullfile(outfolder,fn), in, R, 'CoordRefSysCode', 32642)"
) 
