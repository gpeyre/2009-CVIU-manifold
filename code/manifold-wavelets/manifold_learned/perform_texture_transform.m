function [F,Fint] = perform_texture_transform(M,F,Jmin,dir,options)

% perform_texture_transform - peform non linear wavelet transform
%
% for forward transform
%   D = perform_texture_transform(M,UV,Jmin,+1,options);
% or for backward transform
%   UV = perform_texture_transform(M,D,Jmin,+1,options);
%
%   Copyright (c) 2006 Gabriel Peyré

options.null = 0;
options.sampling = 'uniform';

if ~isfield(options, 'wmax')
    options.wmax = 8;
end

% number of color channels
s = size(M,3);
n = size(F,1);
if dir==1
    % dimensionality of the details
    if isfield(options, 'dim')
        d = options.dim;
    else
        d = 5;
    end
    %%% Forward transform %%%
    UV = F; D = zeros(n,n,d);
    Jmax = log2(n)-1;
    selx = 1:n;
    sely = 1:n;
    % perform multiscale transform
    for j=Jmax:-1:Jmin
        w = (n/2^j)/2;
        options = compute_manifold(M,w,d,options);
        Xs = options.Xs;
        transfo_axis = [-1,1] * (2*(mod(j,2)==0)-1);
        for tau=2:-1:1
            t = transfo_axis(tau);
            if t>0
                % X axis transform
                selxd = selx(end)/2+1:selx(end); selyd = sely;
                selx = 1:selx(end)/2;
            else
                % Y axis transform
                selyd = sely(end)/2+1:sely(end); selxd = selx;
                sely = 1:sely(end)/2;
                UV =  transpose_nd(UV);
            end
            UV1 = UV(1:2:end,:,:);
            UV2 = UV(2:2:end,:,:);
            % perform interpolation
            UVi = perform_texture_interpolation(M,UV1,w,t, options);
            % compute difference in the tangent plane
            IndC = sub2ind( size(M),UVi(:,:,1),UVi(:,:,2) ); % index of center
            IndD = sub2ind( size(M),UV2(:,:,1),UV2(:,:,2) ); % index of point to predict
            Dj = zeros( [size(UV1,1) size(UV1,2) d] );
            for i1=1:size(UV1,1)
                for i2=1:size(UV1,2)
                    ic = IndC(i1,i2);
                    id = IndD(i1,i2);
                    Delta = Xs(:,ic) - Xs(:,id);
                    Dj(i1,i2,:) = reshape( Delta, [1 1 d] );
                end
            end
            Dj = reshape(Dj, [size(UV1(:,:,1)),d]);
            if t<0
                UV1 =  transpose_nd(UV1);
                Dj =  transpose_nd(Dj);
            end
            % assign results
            D(selxd,selyd,:) = Dj;
            UV = UV1;
        end
    end
    % assign remaining coarse
    D(selx,sely,1:2) = UV;
    F = D;
else
    %%% Backward transform %%%
    D = F; % details
    UV = D(1:2^Jmin,1:2^Jmin,1:2);
    d = size(D,3); % dimensionnality
    Jmax = log2(n)-1;
    Fint = {};
    selx = 1:2^Jmin;
    sely = 1:2^Jmin;
    for j=Jmin:Jmax
        w = (n/2^j)/2;
        options = compute_manifold(M,w,d,options);
        Xs = options.Xs;
        atria = options.atria;
        transfo_axis = [-1,1] * (2*(mod(j,2)==0)-1);
        for tau=1:2
            t = transfo_axis(tau);
            if t>0
                % X axis transform
                selxd = selx(end)+1:selx(end)*2; selyd = sely;
                selx = 1:selx(end)*2;
            else
                % Y axis transform
                selyd = sely(end)+1:sely(end)*2; selxd = selx;
                sely = 1:sely(end)*2;
                UV =  transpose_nd(UV);
            end
            % perform interpolation
            UVi = perform_texture_interpolation(M,UV,w,t, options);
            UV1 = UV;       % even point, 
            UV2 = UV1*0;    % odd points, not yet available
            % add details to UVi
            Dj = D(selxd,selyd,:);
            if t<0
                Dj = transpose_nd(Dj);
            end
            IndC = sub2ind( size(M),UVi(:,:,1),UVi(:,:,2) ); % index of center
            for i1=1:size(UVi,1)
                for i2=1:size(UVi,2)
                    ic = IndC(i1,i2);
                    Delta = Dj(i1,i2,:);
                    a = Xs(:,ic) - Delta(:);
                    % find id by closest point
                    [id,distance] = nn_search(Xs', atria, a', 1, 0);
                    [UV2(i1,i2,1),UV2(i1,i2,2)] = ind2sub( size(M), id );
                end
            end
            % interleaved responses
            UV = zeros( 2*size(UV,1),size(UV,2),2); 
            UV(1:2:end,:,:) = UV1;
            UV(2:2:end,:,:) = UV2;
            if t<0
                UV = transpose_nd(UV);
            end
        end
        % save temporary results
        Fint{end+1} = UV;
    end
    F = UV;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UUVV = perform_texture_interpolation(M,UV,w,t, options)

U = UV(:,:,1); V = UV(:,:,2);
Xloc = options.Xloc;
Yloc = options.Yloc;
X = options.X;
Xs = options.Xs;
P = options.P;
Psi = options.Psi;
atria = options.atria;

s = size(M,3);
ww = size(X,1);
w = (ww-1)/2;
if t<0
    % U = U'; V = V';
    X = transpose_nd(X); % transpose the images
end
n1 = size(U,1); m1 = size(U,2);
% perform patch melting
p = n1*m1; % number of current patches
sel1 = 1:n1;        
sel2 = [2:n1,1];
U1 = U(sel1,:,:); V1 = V(sel1,:,:); I1 = sub2ind(size(M),U1(:),V1(:));
U2 = U(sel2,:,:); V2 = V(sel2,:,:); I2 = sub2ind(size(M),U2(:),V2(:));
XX = zeros( ww,ww, s, p );
XX(1:w+1,:,:,:) = X(w+1:end,:,:,I1);
XX(w+1:end,:,:,:) = X(1:w+1,:,:,I2);
if t<0
    XX = transpose_nd(XX); % transpose the images
end
% perform dimension reduction
XXs = reshape( XX, [prod(size(XX(:,:,:,1))), p] );
XXs = XXs - repmat(Psi,[1,size(XXs,2)]);
XXs = P' * XXs;
% search best fit
[I,distance] = nn_search(Xs', atria, XXs', 1, 0);
UU = reshape(Xloc(I),size(U));
VV = reshape(Yloc(I),size(U));
UUVV = zeros( [size(UU) 2] );
UUVV(:,:,1) = UU;
UUVV(:,:,2) = VV;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function options = compute_manifold(M,w,d,options)

% d is the parameters dimensionality

% load patches library
fprintf(['--> Preparing patch libary of size ' num2str(w) ' ... ']);
options.sampling = 'uniform';
[X,Xloc,Yloc] = compute_patch_library(M,w,options);
ww = size(X,1); w = (ww-1)/2;
p0 = size(X,4);
Xs = reshape( X, [prod(size(X(:,:,:,1))), p0] );

% perform sub-sampling to speed up
rand('state',0);
sel = randperm(p0);
m = min(2000,p0); sel = sel(1:m);

% compute global dimension reduction
[P,tmp,v,Psi] = pca(Xs(:,sel),d);
Xs = Xs - repmat(Psi,[1,size(Xs,2)]);
Xs = P' * Xs;

% prepare for nn search
atria = nn_prepare(Xs');

options.P = P;
options.Psi = Psi;
options.Xs = Xs;
options.Xloc = Xloc;
options.Yloc = Yloc;
options.X = X;
options.atria = atria;

fprintf('done.\n');


return;


if isfield(options, 'pca_numvecs')
    m = options.pca_numvecs;
else
    m = 15;
end


% perform dimension reduction
m = min(m,size(Xs,1));
[P,tmp,v,Psi] = pca(Xs,m);
Xs = Xs - repmat(Psi,[1,size(Xs,2)]);
Xs = P' * Xs;
% prepare for nn search
atria = nn_prepare(Xs');

% perform local pca
[options.LocP,options.LocPsi,options.Segm] = perform_local_pca(Xs, options);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = transpose_nd(X)
X = permute(X,[2 1 3 4 5]); % transpose the images