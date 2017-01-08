function [M,errcg] = compute_operator_projection(M, y, pbm, options)


% compute_operator_projection - compute projection on constraint
%
%   If options.niter_cg=1, compute Phi'*(y-Phi*M)
%   otherwise compute projection on {M \ y=Phi*M} using gradient descent.
%
%   Copyright (c) 2008 Gabriel Peyre

niter_cg = getoptions(options, 'niter_cg', 1, 1);
lambda = getoptions(options, 'lambda', .8);

if strcmp(pbm, 'inpainting')
    niter_cg = 1;
end

if niter_cg>1
    options.niter_cg = 1;
    errcg = [];
    for k=1:niter_cg
        h = compute_operator_projection(M, y, pbm, options);
        M = M + lambda*h;
        errcg(end+1) = norm(h, 'fro');
    end
    return;
end


switch pbm
    case 'inpainting'
        U = getoptions(options, 'U', [], 1);
        M(U~=Inf) = y(U~=Inf);
    case 'cs'
        % M = reshape( M(:) + Phi' * ( (Phi*Phi') \ (y-Phi*M(:)) ), n, n);
        n = size(M,1);
        M = callback_sensing_rand( y - callback_sensing_rand(M(:),+1,options), -1,options );
        M = reshape(M,n,n);
    case 'zoom'
        M = perform_zoom_operator( y - perform_zoom_operator(M,+1,options), -1,options );
    otherwise
        error('Unknown method');
end
errcg = [];