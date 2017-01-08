name = 'brick';
name = 'square';
name = 'chessboard';
name = 'grass';
name = 'lines_orientation_smooth';
name = 'locally-parallel';
name = 'text';
name = 'reptilskin';
name = 'fingerprint';
name = 'fabric2';
name = 'dead_leaf';
name = 'locally-parallel';
name = 'locparx';
name = 'grid-circles';
name = 'barb';
name = 'corral';

path(path, '../images/');
path(path, 'images/');

if not(exist('locality'))
    locality = 'nonlocal';
    locality = 'semilocal';
end
use_sparse = 0;
use_eigen = 0;
if not(exist('use_normalized'))
    use_normalized = 1;
end

%%
n = 68;
nn = n^2;
k = 4;  % half size of the window

n0 = [];
mu = 0;
switch name
    case 'fingerprint'
        n0 = 200;
    case 'chessboard'
        n0 = n;
        options.width = n/4;
        options.sigma = 1;
        mu = 0.003; % additional noise level
    case 'corral'
        n0 = 120;
        k=2;
    case 'dead_leaf'
        n0 = 350;
        k = 5;
    case 'fabric2'
        n0 = 200;
    case 'grass'
        n0 = 100;
    case 'lines_orientation_smooth'
        n0 = 80;
        k = 2;
    case 'reptilskin'
        n0 = 180;
        k = 5;
    case 'locally-parallel'
        k = 4;
    case 'text'
        k = 4;
        n0 = 300;
    case 'locparx'
        n0 = 200;
    case 'grid-circles'
        n0 = n;
        
end
m = min((2*k+1)^2 - 4, 30); % for dimreduc

rep = ['results/nonlocal-2d/' name '/']; 
if not(exist(rep))
    mkdir(rep);
end

options.k = k;
kk = (2*k+1)^2;

options.width = 0.4;
options.frequency = 15;

%% load some shape
mu = 0.12;
zeta = 0.07;
M0 = load_image(name,n0, options);
M0 = rescale( crop(M0,n), zeta,1-zeta );
randn('state',666);
M = M0+randn(size(M0))*mu;


%% save original
clf;
imageplot(clamp(M));
warning off;
imwrite(clamp(M), [rep name '-original.png'], 'png');
warning off;


if not(exist('sigma'))
    sigma = 0.1;
    sigma = 0.1;
end
%% compute the mapping
if strcmp(locality,'nonlocal')
    options.ndims = m;
    [H,Proj,Psi] = perform_lowdim_embedding(M,options);
    H = reshape(H, nn, m)';
    % compute euclidian distance
    disp('--> Computing distance matrix.');
    D = sqrt( compute_distance_matrix(H)/kk );
    D = D-diag(diag(D))+diag(Inf*ones(nn,1));
    W0 = exp( -D.^2/(2*sigma^2) );
    W0 = W0-diag(diag(W0));
else
    if not(exist('sigma_sp'))
        sigma_sp = 1.6;
    end
    if not(exist('sigma_va'))
        sigma_va = 0.3;
    end
    fprintf('--> Computing semi-local weights:');
    [tmp,W0] = compute_neighborhood_laplacian(M, sigma_sp, sigma_va); clear tmp;
    D = 1./(1+W0);
end


%% Computation of Laplacians
d = sum(W0,2);
% renormalized smoothing operator
W = W0 ./ repmat( d, [1 nn]);
% symmetric smoothing operator
% Wt = diag(d.^(-1/2))*W0*diag(d.^(-1/2));
Wt = W0 ./ (sqrt(d*d'));
if use_normalized
    % symmetric laplacian operator
    L = eye(nn)-Wt;
else
    % symmetric unnormalized laplacian
    L = diag(d)-W0;
end
% enforce symmetry (avoid <0 eigenvalues)
L = (L+L')/2;

%% compute a reduced gradient operator
disp('--> Computing gradient.');
nb_nn = 30; % number of nearest neighbors
[tmp,I] = sort(D);
I = I(2:nb_nn+1,:);
nv = n^2; ne = nv*nb_nn;
G = sparse(ne,nv);
for i=1:nv
    for j=1:nb_nn
        a = i; b = I(j,i);
        if use_normalized
            G((i-1)*nv + j,a) = +sqrt(W0(a,b))/sqrt(d(a));
            G((i-1)*nv + j,b) = -sqrt(W0(a,b))/sqrt(d(b));
        else
            % connexion between vertex i and I(i,j)
            G((i-1)*nv + j,a) = +sqrt(W0(a,b));
            G((i-1)*nv + j,b) = -sqrt(W0(a,b));
        end
    end
end
L2 = (G'*G)/2;
% L2 should be close to L
norm(L2-L, 'fro');


if use_eigen
    %% compute eigendecomposition
    fprintf('--> Computing eigendecomposition ...'); tic;
    if not(use_sparse)
        [V,delta] = eig(L);
    else
        opts.issym = 1; opts.isreal = 1; opts.disp = 0;
        [V,delta,flags] = eigs(L2,200,'SM', opts);
    end
    fprintf([' done, ' num2str(toc) 'sec.\n']);
    delta=real(diag(delta));
    V = V(:,end:-1:1); delta = delta(end:-1:1);

    % display some eigenvector
    nb = [3 5];
    eig_list = round(linspace(3,25, prod(nb)));
    for i=1:prod(nb)
        A = reshape(V(:,eig_list(i)),n,n);
        imageplot( A,  ['eigv.' num2str(eig_list(i))], nb(1), nb(2), i );
    end
end


%% test for thresholding flow
if strcmp(locality, 'nonlocal')
    tmax = 2.5;
else
    tmax = 8;
end
Tlist = linspace(0,tmax,60); % [0 0.2 0.2 0.3 0.4];
nT = length(Tlist);
Fhard = zeros(n,n,nT); Fsoft = zeros(n,n,nT); Fquad = zeros(n,n,nT); Fgauss = zeros(n,n,nT);
fprintf('--> Computing thresholding flow: ');
if use_eigen
    MV = V'*M(:);
end
for i=1:nT
    progressbar(i,nT);
    t = Tlist(i);
    if use_eigen
        Fhard(:,:,i) = reshape( V*perform_thresholding(MV, t, 'hard'), n,n);
        Fsoft(:,:,i) = reshape( V*perform_thresholding(MV, t, 'soft'), n,n);
        A = diag( 1./(1+delta*t) );
        Fquad(:,:,i) = reshape( V*(A*MV), n,n);
        A = diag( exp( -delta*t ) );
        Fgauss(:,:,i) = reshape( V*(A*MV), n,n);
    else
        Fhard = [];
        Fsoft = [];
        Fquad(:,:,i) = reshape( ( speye(n^2) + t*L2 )\M(:), n,n );
    end
end

if not(use_eigen)
    % solve the heat equation iteratively
    niter = 100;
    dT = tmax/niter;
    Fgauss = zeros(n,n,niter);
    Fgauss(:,:,1) = M;
    for i=2:niter
        f = Fgauss(:,:,i-1);
        Fgauss(:,:,i) = (1-dT) * f + dT*reshape( W*M(:), n,n );
    end
end

%% wavelet denoising
% ne = 2^(ceil(log2(n)));
% Me = perform_image_extension(M,ne);
Me = M;
nTh = 40; Jmin = 3;
ThList = linspace( 2.2,3, nTh )*mu;
MW = perform_atrou_transform(Me,Jmin);
Fwav = zeros(n,n,nTh);
fprintf('--> Computing wavelet flow: ');
for i=1:nTh
    progressbar(i,nTh);
    MWt = perform_thresholding(MW, ThList(i), 'hard');
    Fwav(:,:,i) = perform_atrou_transform(MWt,Jmin);
end

%% local denoising
ns = 40;
slist = linspace(0.1,3,ns);
for i=1:ns
    Floc(:,:,i) = perform_blurring(M,slist(i));
end


%% Non local TV Flow
options.niter = 40;
options.eta = 0.1;
options.lambda_min = 0.01; 
options.lambda_max = 0.06;
fprintf('--> Computing Non Local TV flow: ');
[tmp,Freg] = perform_analysis_regularization(M(:), G, options);
Freg = reshape(Freg,n,n,options.niter);

warning on;
if not(isempty(Fhard))
err = sum( sum( (Fhard - repmat(M0,[1 1 nT]) ).^2, 1 ), 2); err=sqrt( err(:) );
[tmp,k] = min(err); Hhard = Fhard(:,:,k);
if k==1 || k==nT
    warning('Out of bound reached.');
end
err = sum( sum( (Fsoft - repmat(M0,[1 1 nT]) ).^2, 1 ), 2); err=sqrt( err(:) );
[tmp,k] = min(err); Hsoft = Fsoft(:,:,k);
if k==1 || k==nT
    warning('Out of bound reached.');
end
end
err = sum( sum( (Fquad - repmat(M0,[1 1 size(Fquad,3)]) ).^2, 1 ), 2); err=sqrt( err(:) );
[tmp,k] = min(err); Hquad = Fquad(:,:,k);
if k==1 || k==size(Fquad,3)
    warning('Out of bound reached.');
end
err = sum( sum( (Fgauss - repmat(M0,[1 1 size(Fgauss,3)]) ).^2, 1 ), 2); err=sqrt( err(:) );
[tmp,k] = min(err); Hgauss = Fgauss(:,:,k);
if k==1 || k==size(Fgauss,3)
    warning('Out of bound reached.');
end
err = sum( sum( (Freg - repmat(M0,[1 1 options.niter]) ).^2, 1 ), 2); err=sqrt( err(:) );
[tmp,k] = min(err); Hreg = Freg(:,:,k);
if k==1 || k==options.niter
    warning('Out of bound reached.');
end
err = sum( sum( (Fwav - repmat(M0,[1 1 size(Fwav,3)]) ).^2, 1 ), 2); err=sqrt( err(:) );
[tmp,k] = min(err); Hwav = Fwav(:,:,k);
if k==1 || k==size(Fwav,3)
    warning('Out of bound reached.');
end
err = sum( sum( (Floc - repmat(M0,[1 1 size(Floc,3)]) ).^2, 1 ), 2); err=sqrt( err(:) );
[tmp,k] = min(err); Hloc = Floc(:,:,k);
if k==1 || k==size(Floc,3)
    warning('Out of bound reached.');
end

pnoisy = psnr(M,M0,1);
if not(isempty(Fhard))
    phard = psnr(Hhard,M0,1);
    psoft = psnr(Hsoft,M0,1);
else
    phard = -1; psoft = -1;
    Hhard = zeros(n); Hsoft = zeros(n);
end
pquad = psnr(Hquad,M0,1);
pgauss = psnr(Hgauss,M0,1);
preg = psnr(Hreg,M0,1);
pwav = psnr(Hwav,M0,1);
ploc = psnr(Hloc,M0,1);

clf;
% imageplot(M0, 'Original', 2,3,1);
imageplot(clamp(M), ['Noisy, psnr=' num2str(pnoisy,3)], 2,3,1);
imageplot(clamp(Hwav), ['Wav, psnr=' num2str(pwav,3)], 2,3,2);
imageplot(clamp(Hreg), ['Loc, psnr=' num2str(ploc,3)], 2,3,6);
% imageplot(clamp(Hhard), ['Hard, psnr=' num2str(phard,3)], 2,3,2);
imageplot(clamp(Hsoft), ['Soft, psnr=' num2str(psoft,3)], 2,3,3);
imageplot(clamp(Hquad), ['Quad, psnr=' num2str(pquad,3)], 2,3,4);
imageplot(clamp(Hgauss), ['Gauss, psnr=' num2str(pgauss,3)], 2,3,5);
imageplot(clamp(Hreg), ['Reg, psnr=' num2str(preg,3)], 2,3,6);
saveas(gcf, [rep name '-results.png'], 'png');

filename = [rep name '-results.txt'];
fid = fopen(filename, 'a');
if strcmp(locality, 'nonlocal')
    fprintf(fid, '--> denoising, n=%d, locality=%s, sigma=%.3f, normalized=%d\n', n, locality, sigma, use_normalized);
else
    fprintf(fid, '--> denoising, n=%d, locality=%s, sigma_sp=%.3f, sigma_va=%.3f, normalized=%d\n', n, locality, sigma_sp, sigma_va, use_normalized);
end
    
fprintf(fid, 'Noisy: PSNR=%.3f\n', pnoisy);
fprintf(fid, 'Hard:  PSNR=%.3f\n', phard);
fprintf(fid, 'Soft:  PSNR=%.3f\n', psoft);
fprintf(fid, 'Quad:  PSNR=%.3f\n', pquad);
fprintf(fid, 'Gauss: PSNR=%.3f\n', pgauss);
fprintf(fid, 'TV:    PSNR=%.3f\n', preg);
fprintf(fid, 'Wav:   PSNR=%.3f\n', pwav);
fprintf(fid, 'Loc:   PSNR=%.3f\n', ploc);
fclose(fid);

addstr = '';
if use_normalized==0
    addstr = '-unorm';
end
if strcmp(locality, 'semilocal')
    addstr = [addstr '-semiloc'];
else
    addstr = [addstr '-nonloc'];
end

warning off;
imwrite(clamp(M0), [rep name '-original.png'], 'png');
imwrite(clamp(M), [rep name '-noisy.png'], 'png');
imwrite(clamp(Hhard), [rep name '-denoise' addstr '-hard.png'], 'png');
imwrite(clamp(Hsoft), [rep name '-denoise' addstr '-soft.png'], 'png');
imwrite(clamp(Hquad), [rep name '-denoise' addstr '-quad.png'], 'png');
imwrite(clamp(Hgauss), [rep name '-denoise' addstr '-gauss.png'], 'png');
imwrite(clamp(Hreg), [rep name '-denoise' addstr '-tv.png'], 'png');
imwrite(clamp(Hwav), [rep name '-denoise-wav.png'], 'png');
imwrite(clamp(Hloc), [rep name '-denoise-loc.png'], 'png');
warning on;
