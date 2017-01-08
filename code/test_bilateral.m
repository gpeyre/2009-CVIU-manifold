% test for bilateral filter

n = 256;
name = 'lena';
M = load_image(name);
M = rescale(crop(M,n));

sigmaSpatial = 1/16;
sigmaRange = 1/10;

samplingSpatial = 2*sigmaSpatial;
samplingRange = 2*sigmaRange;

tic;
M1 = perform_bilateral_filter( M, M, ...
                           mmin(M), mmax(M), ...
                           sigmaSpatial, sigmaRange, ...
                           samplingSpatial, samplingRange );
toc;

imageplot({M M1});