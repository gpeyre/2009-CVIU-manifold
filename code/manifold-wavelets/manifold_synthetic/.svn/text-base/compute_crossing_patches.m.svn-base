function [M, delta1_list,theta1_list, delta2_list,theta2_list] = compute_crossing_patches(w, options)

% [M, delta1_list,theta1_list, delta2_list,theta2_list] = compute_crossing_patches(w, options);

options.null = 0;
if isfield(options, 'n_delta')
    n_delta = options.n_delta;
else
    n_delta = 11; 
end
if isfield(options, 'n_theta')
    n_theta = options.n_theta;
else
    n_theta = 12;
end
if isfield(options, 'sigma')
    sigma = options.sigma;
else
    sigma = 0.1;
end
sigma = sigma*w;
if isfield(options, 'rescale')
    rescale = options.rescale;
else
    rescale = 1;
end
% the window to weight the distance
if isfield(options, 'bell')
    bell = options.bell;
else
    bell = 'constant';
end

switch lower(bell)
    case 'constant'
        B = ones(w);
    case 'sine'
        x = linspace(0,1,w);
        x = (1-cos(2*pi*x))/2;
        B = x'*x;
    otherwise
        error('Unknown bell shape');
end

dmax = sqrt(2)*w/2+2*sigma;
eta = 2;
t = linspace(0,pi-pi/n_theta,n_theta); t = t(:);
d = linspace(-1,1,n_delta);
d = sign(d).*abs(d).^eta; 
d = d(:)*dmax;

[delta1_list,theta1_list,delta2_list,theta2_list] = ndgrid( d, t, d, t );
theta1_list = theta1_list(:);
delta1_list = delta1_list(:);
theta2_list = theta2_list(:);
delta2_list = delta2_list(:);
p = length(delta1_list);

x = linspace( -(w-1)/2,(w-1)/2, w );
[Y,X] = meshgrid(x,x);

c1 = repmat( reshape(cos(theta1_list),[1 1 p]), [w w 1] );
s1 = repmat( reshape(sin(theta1_list),[1 1 p]), [w w 1] );
c2 = repmat( reshape(cos(theta2_list),[1 1 p]), [w w 1] );
s2 = repmat( reshape(sin(theta2_list),[1 1 p]), [w w 1] );
d1 = repmat( reshape(delta1_list,[1 1 p]), [w w 1] );
d2 = repmat( reshape(delta2_list,[1 1 p]), [w w 1] );
X = repmat( X, [1 1 p] );
Y = repmat( Y, [1 1 p] );

A1 = c1.*X + s1.*Y - d1;
A1 = 2*( 1 - exp( -A1.^2 / (2*sigma^2) ) ) - 1;
A2 = c2.*X + s2.*Y - d2;
A2 = 2*( 1 - exp( -A2.^2 / (2*sigma^2) ) ) - 1;
M = min(A1,A2);


% normalize
if rescale
    B = repmat(B,[1,1,p]);
    s = sqrt( sum( sum(B .* (M.^2),1), 2) );
    M = M ./ repmat(s,[w,w,1]);
end

return;

% compute images
for k=1:p
    t = theta_list(k);
    d = delta_list(k);
    y = cos(t)*X + sin(t)*Y - d;
    if strcmp(manifold_type, 'edges')
        A = tanh(y/sigma);
    else
        A = 2*( 1 - exp( -y.^2 / (2*sigma^2) ) ) - 1;
    end
    M(:,:,k) = A;
end
