function M = perform_local_manifold_projection(M, D, options)

% perform_local_manifold_projection - project patch on nearest manifold point
%
%   M = perform_local_manifold_projection(M, D, options);
%
%   set options.dr to use PCA projection to speed up distance computing.
%
%   Copyright (c) 2007 Gabriel Peyre


d = size(D,1);
w = getoptions(options, 'w', sqrt(d));
d = w^2;

% for inpainting
U = getoptions(options, 'U', []);
% for affine matching
niterdec = getoptions(options, 'niterdec', 3);

% extract patches
H0 = compute_all_patch(M,w, options);
if isempty(U)
    H = H0;
else
    % only patches crossing the inpainted region
    HU = compute_all_patch(U,w, options);
    sel = find( sum(sum(HU==Inf,1),2)>0 );
    H = H0(:,:,:,sel);
end
q = size(H,4);
Y = reshape( H, d,q );

if isstr(D) && strcmp(D, 'regular')
    % special case of regular manifold
    % look for a decomposition P=a*x+b*y+c*z
    a = ones(1,q);
    b = zeros(1,q);
    c = zeros(1,q);
    [x,y] = meshgrid(1:w,1:w); x = x(:); y = y(:); z = x*0+1;
    x = repmat(x, [1 q]);
    y = repmat(y, [1 q]);
    z = repmat(z, [1 q]);
    err = [];
    for k=1:10
        Y1 = Y - repmat(b,[d 1]).*y - repmat(c,[d 1]).*z;
        a = sum(Y1.*x)./sum(x.*x);
        Y1 = Y - repmat(a,[d 1]).*x - repmat(c,[d 1]).*z;
        b = sum(Y1.*y)./sum(y.*y);
        Y1 = Y - repmat(a,[d 1]).*x - repmat(b,[d 1]).*y;
        c = sum(Y1.*z)./sum(z.*z);
        err(end+1) = sum( sum( (Y - repmat(a,[d 1]).*x - repmat(b,[d 1]).*y - repmat(c,[d 1]).*z ).^2  ) );
    end
    if err(end)>err(1)
        warning('Convergence problem');
    end
    H0 = reshape( repmat(a,[d 1]).*x - repmat(b,[d 1]).*y - repmat(c,[d 1]).*z, w, w, 1, q );

else
    % look for decomposition y=a*D+b
    a = ones(1,q);
    b = zeros(1,q);
    for k=1:niterdec
        Y1 = (Y - repmat(b,[d 1])) ./ repmat(a,[d 1]);
        dist = compute_distance_to_points(Y1, D, options);
        [tmp,I] = min(dist);
        X = D(:,I);
        if niterdec>1
            % perform linear regression for all the patches
            err = [];
            for l=1:30
                b = mean( Y-repmat(a,[d 1]).*X);
                b = clamp(b,0,1);
                e = sum( X.^2 ); e(e<1e-9) = 1;
                a = sum( X.*( Y-repmat(b,[d 1]) ) ) ./ e;
                a = clamp(a,.1,1);
                err(end+1) = sum( sum( (X.*repmat(a,[d 1]) + repmat(b,[d 1])-Y).^2 ));
            end
            if err(end)>err(1)
                warning('Convergence problem');
            end
        end
    end
    % reconstruct by averaging
    H1 = reshape( D(:,I).*repmat(a,[d 1]) + repmat(b,[d 1]), w, w, 1, q );
end

if not(isempty(U))
    % add back the patch not crossing inpainted region
    H0(:,:,:,sel) = H1; H1 = H0;
end
% average to reconstruct back
M = compute_all_patch(H1, size(M,1), options);
