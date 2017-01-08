function M = perform_zoom_operator(M, dir, options)

h = ones(3);
options.bound = 'sym';

s = getoptions(options, 'subsampling', 2);

filter = getoptions(options, 'filter', 'gaussian');

w = s+1; % width of filtering
n = size(M,1);
switch filter
    case 'box'
        h = ones(w);
    case 'triangle'
        h = conv2(ones(w), ones(w));
    case 'gaussian'
        h = compute_gaussian_filter([21 21], s/(2*n), [n n]);
end


if dir==1
    M = perform_convolution(M, h, options); 
    M = M(1:s:end,1:s:end);
else
    A = zeros(size(M)*s);
    A(1:s:end,1:s:end) = M;
    M = perform_convolution(A, h, options); 
end