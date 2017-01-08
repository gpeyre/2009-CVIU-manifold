function H = perform_lowdim_embedding_1d(x,options)

options.null = 0;
if isfield(options, 'bound')
    bound = options.bound;
else
    bound = 'per';
end
if isfield(options, 'k')
    k = options.k;
else
    k = 4;
end

x = x(:);
n = size(x,1);

[Y,X] = meshgrid(1:n,-k:k);
t = X+Y;
if strcmp(bound,'per')
    t = mod(t-1,n)+1;
else
    t(t<1) = 1-t(t<1);
    t(t>n) = 2*n+1-t(t>n);
end
H = x(t);