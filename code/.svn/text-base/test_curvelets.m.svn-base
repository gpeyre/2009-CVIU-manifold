% test for the curvelet transform
path(path, 'curvelab/');
path(path, 'curvelab/curvelab-toolbox/fdct_wrapping_matlab/');

name = 'barb';
n = 256;
M = load_image(name);
M = rescale(crop(M,n));

finest = 2;
nbscales = log2(n)-2;
nbangles_coarse = 16;
C = fdct_wrapping(M, 1, finest, nbscales, nbangles_coarse);

M1 = ifdct_wrapping(C, 1);

[X_rows, X_cols, F_rows, F_cols, N_rows, N_cols] = fdct_wrapping_param(C,size(M,1),size(M,2));