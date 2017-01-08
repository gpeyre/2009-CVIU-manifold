function M = display_packed_edges(Param,n, options)

options.null = 0;


Theta = Param.Theta;
Delta = Param.Delta;
m = size(Delta,1);
M = zeros(n);

w0 = n/m;
w = w0+1;

Delta = Delta*n/m;

for i=1:m
    for j=1:m
        selx = (i-1)*w0+1:i*w0;
        sely = (j-1)*w0+1:j*w0;
        E = compute_edge(Theta(i,j),Delta(i,j),w, options);
        M(selx,sely) = E(1:end-1,1:end-1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function E = compute_edge(t,d,w, options)


if isfield(options, 'manifold_type')
    manifold_type = options.manifold_type;
else
    manifold_type = 'edges';
end

if isfield(options, 'sigma')
    sigma = options.sigma;
else
    sigma = 0.05;
end
sigma = sigma*w;

x = linspace( -(w-1)/2,(w-1)/2, w );
[Y,X] = meshgrid(x,x);
y = cos(t)*X + sin(t)*Y - d;
if strcmp(manifold_type, 'edges')
    E = tanh(y/sigma);
else
    E = 2*( 1 - exp( -y.^2 / (2*sigma^2) ) ) - 1;
end