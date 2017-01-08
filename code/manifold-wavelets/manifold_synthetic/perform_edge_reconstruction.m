function M = perform_edge_reconstruction(Param, options)

% M = perform_edge_reconstruction(Param, options);

options.null = 0;
if isfield(options, 'w')
    w = options.w;
else
    w = 11;
end

if isfield(options, 'manifold_type')
    manifold_type = options.manifold_type;
else
    manifold_type = 'edges';
end

if isfield(options, 'sigma')
    sigma = options.sigma;
else
    sigma = 0.1;
end
sigma = sigma*w;

Delta = Param.Delta;

if strcmp(manifold_type, 'edges')
    M = tanh(Delta/sigma);
else
    M = 2*( 1 - exp( -Delta.^2 / (2*sigma^2) ) ) - 1;
end


return;

Theta = Param.Theta;
if isfield(Param, 'Lambda')
    Lambda = Param.Lambda; 
else
    Lambda = Delta*0+1;
end
if isfield(Param, 'Mu')
    Mu = Param.Mu;
else
    Mu = Delta*0;
end

n = size(Theta,1);

if isfield(options, 'n_full')
    n_full = options.n_delta;
else
    n_delta = 16; 
end
if isfield(options, 'n_theta')
    n_theta = options.n_theta;
else
    n_theta = 24;
end
s = w^2;

options.rescale = 0;
options.add_flat = 1;
[H,delta_list,theta_list] = compute_edge_patches(w,options);


M = zeros(n);
for i=1:n
    for j=1:n
        % distance on circle
        t = Theta(i,j);
        d = Delta(i,j);
        if d==Inf
            M(i,j) = -1;
        elseif d==-Inf
            M(i,j) = 1;
        else
            Et = min([abs(theta_list-t), abs(theta_list-t+2*pi), abs(theta_list-t-2*pi)], [], 2);
            E = abs(delta_list-d) + Et;
            [tmp,I] = min(E);
            x = H(:,:,I);
            x = x((w+1)/2,(w+1)/2);
            M(i,j) = x*Lambda(i,j)+Mu(i,j);
        end
    end
end