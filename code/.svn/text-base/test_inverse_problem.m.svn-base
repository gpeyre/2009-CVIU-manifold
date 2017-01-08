% test for inverse problem using manifold model

path(path, '../images/');

if not(exist('pbm'))
    pbm = 'zoom';
    pbm = 'cs';
    pbm = 'inpainting';
end

if not(exist('name'))
    name = 'fnoise';
    name = 'piece-regular';
    name = 'disk';
    name = 'geometrical';
    name = 'fingerprint';
    name = 'reptilskin';
    name = 'barb';
    name = 'corral';
    name = 'boat';
    name = 'lena';
end

if not(exist('method'))
    method = 'sparsity-haar';
    method = 'sparsity-curvelets';
    method = 'sparsity-dct-orth';
    method = 'sparsity-dct';
    method = 'manifold';
    method = 'sparsity-wavelets';
    method = 'nl-means';
    method = 'tv';
end


rep = ['results/' pbm '/' name '/'];
repiter = [rep method '/iter/'];
if not(exist(rep))
    mkdir(rep);
end
if not(exist(repiter))
    mkdir(repiter);
end

disp(['---- Processing ' name ', method=' method ' ----']);

% kind of manifold
nbdims = 2;
switch name
    case {'geometrical', 'disk'}
        manifold = 'edges';
    case 'fnoise'
        manifold = 'regular';
    case 'piece-regular';
        nbdims = 1;
    case 'fingerprint'
        manifold = 'oscillations';
    case {'reptilskin' 'corral'}
        manifold = 'exemplar';
    otherwise
        manifold = 'exemplar';
%        error('Unknown name.');
end

if nbdims==2
    n = 256;
    N = n^2;
else
    n = 2048;
    N = n;
end

n0 = [];
switch name
    case 'fingerprint'
        n0 = 180;
    case 'reptilskin'
        n0 = 180;
    case 'coral'
        n0 = 140;
end

options.Jgeometrical = 2;
options.alpha = 2.8;
options.nbdims = nbdims;

if nbdims==2
    M0 = load_image(name, n0, options);
    M0 = rescale( crop(sum(M0,3),n) );
    if strcmp(name, 'geometrical') || strcmp(name, 'disk')
        M0 = perform_blurring(M0,2);
    end
else
    M0 = rescale(load_signal(name, n0, options) );
end
options.M0 = M0;


%% parameters for the methods
switch method
    case 'sparsity-haar'
        options.wavelet_type = 'daubechies'; options.wavelet_vm = 1;    
        mca_callback = @callback_atrou;
    case 'sparsity-wavelets'
        if nbdims==2
            options.wavelet_type = 'biorthogonal'; 
        else
            options.wavelet_type = 'daubechies';
        end
        options.wavelet_vm = 3;
        mca_callback = @callback_atrou;
    case 'sparsity-curvelets'
        mca_callback = @callback_curvelets;
    case 'sparsity-dct'
        mca_callback = @callback_localdct;
        options.w = 16;
        options.n = n;
        options.dct_type = 'redundant';
        options.remove_lowfreq = 0;
    case 'sparsity-dct-orth'
        mca_callback = @callback_localdct;
        options.w = 16;
        options.n = n;
        options.dct_type = 'orthogonal';
        options.remove_lowfreq = 0;
    case 'nl-means'
        options.k = 3;          % half size for the windows
        options.T = 0.08;       % width of the gaussian, relative to max(M(:))  (=1 here)
        options.max_dist = 8;  % search width, the smaller the faster the algorithm will be
        options.ndims = 25;     % number of dimension used for distance computation (PCA dim.reduc. to speed up)
        options.do_patchwise = 0;
    case 'tv'
        options.order = 2;
        options.bound = 'sym';
end


%% perform measurement
switch pbm
    case 'cs'
        %% compressive sampling
        if not(exist('p'))
            p = round(N/4);
        end
        if 0
            Phi = compute_compressed_sensing_matrix(p,N, 'random', 'normalize');
            y = Phi*M0(:);
        else
            options.n = N; options.p = p;
            y = callback_sensing_rand(M0(:), +1, options);
        end
        U = []; % no mask
    case 'inpainting'
        %% inpainting
        inpaint_mode = 'mask';
        inpaint_mode = 'missing';
        if nbdims==2
            switch inpaint_mode
                case 'mask'
                    if not(exist('U'))
                        options.r = 3;
                        options.mode = 'line';
                        U = grab_inpainting_mask(M0, options);
                    end
                case 'missing'
                    if not(exist('nmiss'))
                        nmiss = .8; % ratio of missing pixels
                    end
                    rand('state', 12345);
                    S = rand(n); [tmp,I] = sort(S(:));
                    U = zeros(n);
                    U(I(1:round(n^2*nmiss))) = Inf;
            end
        else
            x = linspace(0,1,n);
            c = [.2 .5 .7];
            c = 0.05:.06:.95;
            r = .01*.8;
            U = zeros(n,1);
            for i=1:length(c)
                U( abs(x-c(i))<r ) = Inf;
            end
        end
        ninpaint = sum(U(:)==Inf);
        % erase pixels
        y = M0; y(U==Inf)=rand(ninpaint,1);
    case 'zoom'
        if not(exist('s'))
            s = 2;
        end
        options.subsampling = s; % sub-sampling
        y = perform_zoom_operator(M0, +1, options);
        U = []; % no mask
        %% add some noise
        noise = .0;
        randn('state', 123456);
        y = y + noise*randn(size(y));
end
options.U = U;


%% save original input data y
if strcmp(pbm, 'inpainting')
    M1 = M0; M1(U==Inf) = 1;
    M2 = M0; M2(U==Inf) = 0;
    warning off;
    imwrite(cat(3,M1,M2,M2), [rep name '-mask-' num2str( round(nmiss*100) ) '.png'], 'png');
    warning on;
elseif strcmp(pbm, 'zoom')
    warning off;
    imwrite(rescale(y), [rep name '-input-s' num2str(s) '.png'], 'png');
    warning on;
end

% load manifold model
if strcmp(method, 'manifold')
    w = 10;
    d = w^2; % dimensionality
    options.w = w;
    options.n_delta = 12;
    options.n_theta = 22;
    options.sigma = 0.02;
    options.sigma = 0.005;
    options.mode = manifold;
    if strcmp(manifold, 'regular')
        D = 'regular';
    elseif strcmp(manifold, 'exemplar')
        MM = load_image(name, n0, options);
        MM = rescale( MM(1:64,1:64) );
        H = compute_all_patch(MM,w, options);        
        q = d*8; 
        sel = randperm( size(H,4) ); sel = sel(1:q);
        D = reshape( H(:,:,:,sel), [d q] );
    else
        D = compute_synthetic_dictionary(manifold, w, options);
    end
    % dimensionality reduction
    options.dr = min(25,d);
end

%% initial guess
switch pbm
    case 'cs'
        if exist('Phi')
            M = Phi'*y;
        else
            M = callback_sensing_rand(y, -1, options); M = reshape(M,n,n);
        end
        if nbdims==2
            M = reshape(M,n,n);
        end
        M = clamp(M-mean(M(:))+mean(M0(:)));
        
        sigma = linspace(6,1,2);
        for i=1:length(sigma)
            progressbar(i,length(sigma));
            M = perform_blurring(M,sigma(i));
            options.niter_cg = 10;
            M = compute_operator_projection(M, y, pbm, options);
            M = clamp(M-mean(M(:))+mean(M0(:)));
        end
    case 'inpainting'
        M = y;
        % initialize with linear reconstruction
        for i=1:40
            M = perform_blurring(M,3);
            M(U~=Inf) = M0(U~=Inf);
        end
    case 'zoom'
        %% bilinear interpolation
        xi = 1:1/s:(n/s+1-1/s);
        xi(xi>n/s) = 2*n/s-xi(xi>n/s);
        [Xi,Yi] = meshgrid(xi,xi);
        [X,Y] = meshgrid(1:n/s,1:n/s);
        M = interp2(X,Y,y, Xi,Yi);
        %% compute projection of the first iteration
        options.niter_cg = 60;
        M = compute_operator_projection(M, y, pbm, options);
end



%% save initial guess
if nbdims==2
    warning off;
    imwrite(rescale(M), [repiter name '-iter-000.png'], 'png');
    warning on;
    warning off;
    imwrite(rescale(M0), [rep name '-original.png'], 'png');
    warning on;
else
    lw = 2;        
    clf; 
    hold on;
    h = plot(x, M.*(1+U)); set(h, 'LineWidth', lw);
    V = 1-(U==Inf); V(V==1) = Inf;
    h = plot(x, M0.*(1+V), 'r'); set(h, 'LineWidth', 1);
    axis tight; box on;
    hold off;
	saveas(gcf, [repiter name '-iter-000.png'], 'png');
end

% set to 1 for no affinity transform
options.niterdec = 3;
if strcmp(name, 'disk')
    options.niterdec = 1;
elseif strcmp(manifold, 'exemplar')
    options.niterdec = 1;
end

options.bound = 'sym';
% options.bound = 'per';

%% optimization parameters schedule
if strcmp(method, 'nl-means')
    niter = 60;
    tlist = linspace(.08,.015,niter); % width of filtering
    % it is required that tau*mu<1
    if strcmp(pbm, 'cs')  
        mumax = 1;
        mumax = .1;     
        tlist = linspace(.05,.015,niter);
    elseif strcmp(pbm, 'zoom')
        mumax = 1;
    else
        mumax = .01;
    end
    mulist = linspace(mumax, 0, niter); 
    tau = .1; % discretization time steps
elseif strcmp(method, 'tv')
    %% total variation 
    niter = 200;
    tau = 1; % gradient descent step
    tau = .1;
    if strcmp(pbm, 'cs')       
        mumax = .1;
    elseif strcmp(pbm, 'inpainting')
        tau = .05;
        mumax = .1;
    else
        mumax = .01;
    end
    mulist = linspace(mumax, 0, niter);
else
    %% sparsity
    thresh = 'soft';
    thresh = 'hard';
    niter = 200;
    tlist = linspace(1,0,niter).^4;
end

Mtemp = [];
if strcmp(pbm, 'cs') && strcmp(method, 'nl-means')
    Mtemp = load_image( ['cs-templates/' name '-sparsity-wavelets-p' num2str( round(N/p) )] );
    Mtemp = rescale(crop(Mtemp,n));
    % options.Ma = rescale(crop(Mtemp,n));
    % options.Ma = M0;
    M = Mtemp;
end
    
Msvg = [];
err = [];
for i=1:niter
    progressbar(i,niter);
    % perform local nearest neighbor projection
    if strcmp(method, 'manifold')
        M = perform_local_manifold_projection(M, D, options);
    elseif strcmp(method, 'nl-means')
        %% Nlmeans
        mu = mulist(i); % regularization weight by NLmeans
        options.T = tlist(i);
        [M1,Wx,Wy] = perform_nl_means(M, options);
       
        if 1
            
            %% compute projection
            options.niter_cg = 30;
            Ma = compute_operator_projection(M1, y, pbm, options);
            Ma = clamp(Ma-mean(Ma(:))+mean(M0(:)));            
            M = Ma;

%            options.niter_cg = 1;
%            Ma = M + .5*compute_operator_projection(M, y, pbm, options);
%            Ma = clamp(Ma-mean(Ma(:))+mean(M0(:)));
            
            % M = (1-mu)*Ma + mu*M1;
            
                        
        else
            g = M1-M;

            %% compute projection
            options.niter_cg = 1;
            h = compute_operator_projection(M, y, pbm, options);

            M = M + tau * ( h+mu*g );

            %% post projection
            options.niter_cg = 10;
            M = compute_operator_projection(M, y, pbm, options);
        end
        
    elseif strcmp(method, 'tv')
        
        % perform projection
        M = compute_operator_projection(M, y, pbm, options);
        % compute prox for regularization mu
        mu = mulist(i);
        errtv = []; nitertv = 50; Ma = M;
        tautv = .02;
        for itv = 1:nitertv
            options.order = 1;
            mu = mulist(i); % regularization weight by TV
            g = grad(M, options);
            d = sqrt(sum(g.^2, 3));
            d = max(d, 1e-6);
            g = g ./ repmat(d, [1 1 2]);
            g = div(g, options);
            errtv(end+1) = .5*norm( Ma-M, 'fro' )^2 + mu*sum(d(:));    
            M = M + tautv * ( (Ma-M) + mu*g );         
        end
        if errtv(1)<errtv(end)
            warning('Pbm with TV sub-routine');
        end
        
        if 0
        
        %% compute projection
        if strcmp(pbm, 'inpainting')
%            M1 = compute_operator_projection(M, y, pbm, options);
            h = zeros(n);
            h(U~=Inf) = M(U~=Inf) - M0(U~=Inf);
            M = M + tau * ( -h + mu*g ); 
%            M = compute_operator_projection(M, y, pbm, options);
        else
            options.niter_cg = 1;
            h = compute_operator_projection(M, y, pbm, options);
            M = M + tau * ( h+mu*g );
            options.niter_cg = 10;
            M = compute_operator_projection(M, y, pbm, options);
        end
        
        end
    else
        %% sparsity
        MW = feval(mca_callback, M,+1,options);
        % thresholding
        MW = perform_thresholding(MW,tlist(i),thresh);
        % transform back
        M = feval(mca_callback, MW,-1,options);
        %% compute projection
        options.niter_cg = 30;
        M = compute_operator_projection(M, y, pbm, options);
    end
    
    if exist('mu') && mu==0
        %% compute final projection
        options.niter_cg = 40;
        M = compute_operator_projection(M, y, pbm, options);
    end
    
    M = clamp(M-mean(M(:))+mean(M0(:)));
    
    % save the result
    if nbdims==2
        warning off;
        % imwrite(rescale(M), [repiter name '-iter-' num2string_fixeddigit(i,3) '.png'], 'png');
        warning on;
        % Msvg(:,:,end+1) = M;
        clf; imageplot(clamp(M)); drawnow;
    else
        clf; h = plot(x, M); axis tight; set(h, 'LineWidth', lw);
        saveas(gcf, [repiter name '-iter-' num2string_fixeddigit(i,3) '.png'], 'png');
    end
    
    M = clamp(M-mean(M(:))+mean(M0(:)));
    err(i) = norm(M-M0, 'fro');
    if i==1 || err(i)<=min(err(1:i-1))
        Mbest = M;
    end
end

addstr = [];
if strcmp(pbm, 'zoom')
    addstr = ['-s' num2str(s)];
elseif strcmp(pbm, 'cs')
    addstr = ['-p' num2str( round(N/p) )];
elseif strcmp(pbm, 'inpainting')
    addstr = ['-r' num2str( round(nmiss*100) )];
end
warning off;
imwrite(rescale(Mbest), [rep name '-' method addstr '.png'], 'png');
warning on;
% compute_movie_file(Msvg, [rep name '-' method addstr '-iter.gif']);

disp(['PSNR=' num2str(psnr(M0,Mbest)) '.']);

% record in a file the parameters
if strcmp(pbm, 'cs') || strcmp(pbm, 'zoom') || strcmp(pbm, 'inpainting')
    fid = fopen([rep name '-results.txt'], 'at'); % '  p=' num2str(p)
    fprintf(fid, ['----> ' method addstr '\t\tpsnr=' num2str(psnr(M0,Mbest)) '\n']);
    if strcmp(method, 'sparsity-wavelets')  
        fprintf(fid, ['thresh=' thresh '\n']);
    end
    if strcmp(method, 'nl-means')  
        fprintf(fid, ['  tmin =' num2str(tlist(end)) '  tmax =' num2str(tlist(1)) '\n']);
        fprintf(fid, ['  mumin=' num2str(mulist(end)) '  mumax=' num2str(mulist(1)) '\n']);
    end
    fclose(fid);
end