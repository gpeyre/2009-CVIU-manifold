function y = callback_curvelets(x, dir, options)

path(path, 'curvelab/');


if dir==1
    n = size(x,1);
    finest = 2;
    nbscales = log2(n)-2;
    nbangles_coarse = 16;
    y = fdct_wrapping(x, 1, finest, nbscales, nbangles_coarse);
else
    y = ifdct_wrapping(x, 1);
end