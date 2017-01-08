% test for inpainting

n0 = 140;
n = 100;
name = 'fingerprint3';

problem = 'inpainting';
problem = 'cs';

rep = ['results/locpar-' problem '/'];
if not(exist(rep))
    mkdir(rep);
end

M = load_image(name, n0);
M = rescale( crop(M,n) );

if strcmp(problem, 'inpainting')
    if not(exist('U'))
        options.r = 3;
        options.mode = 'line';
        U = grab_inpainting_mask(M, options);
    end
    mask = (U==Inf);
else
    % number of measurements
    P = round(n^2/4);
    A = rand(n);
    k = 6;
    A(end/2-k:end/2+k,end/2-k:end/2+k) = -Inf;
    A = fftshift(A);
    v = sort(A(:));
    mask = 1-(A<=v(P));
end

% parameter of the dictionary
w = 8;
options.fmin = 3;
options.fmax = 6;
options.nfreq = 12;
options.nphase = 12;
options.ntheta = 16;
options.gamma = .35;
D = compute_locpar_dictionary(w,options);

y = M;


%% Linear inpainting
options.mask = mask;
if strcmp(problem, 'inpainting')
    U = @callback_inpainting;
else    
    U = @callback_fourier_cs;
end
y = feval(U,M,+1,options);
M1 = feval(U,y,-1,options);


% inpainting
niter = 30;
filename = {}; Msvg = {};
err = [];
for i=1:niter
    if i<10 && strcmp(problem,'inpainting')
        filename{end+1} = ['iter' num2str(i)];
        M2 = M1; M2(mask==1) = perform_histogram_equalization( M1(mask==1), y(mask==0) );
        Msvg{end+1} = M1;
    end 
    options.sub = 2;
    if i>.6*niter
        options.sub = 2;
    end
    if i>.8*niter
        options.sub = 2;
    end
    if i>=niter-7
        options.sub = 1;
    end
    progressbar(i,niter);
    M1 = perform_patch_manifold_projection(M1,D,options);
    M1 = clamp(M1);
    % impose known values
    if strcmp(problem, 'inpainting')
        M1(mask==0) = y(mask==0);
    else
        M1 = fft2(M1)/n;
        M1(mask==0) = y(mask==0);
        M1 = real(ifft2(M1))*n;
    end
    M1 = clamp(M1);
    err(end+1) = snr(M,M1);
    clf; imageplot(M1); drawnow;
end

clf;
imageplot(M1, ['SNR=' num2str(snr(M,M1))]);

if strcmp(problem,'inpainting')
    M2 = M1; M2(mask==1) = perform_histogram_equalization( M1(mask==1), y(mask==0) );
    clf;
    imageplot({M1 M2}, {['SNR=' num2str(snr(M,M1))] ['SNR=' num2str(snr(M,M2))]});
else
    M2 = M1;
end


filename{end+1} = 'manifold';
Msvg{end+1} = M2;
filename{end+1} = 'original';
Msvg{end+1} = M;


options.clamp = 1;
options.base_str = [rep name '-'];
save_image(Msvg,filename,options);

%% DCT options
options.n = n;
options.w = 16; % width of patches
options.q = options.w/4;
options.dct_type = 'orthogonal2';
options.dct_type = 'redundant';
options.dct_type = 'orthogonal4';
options.dct_type = 'gabor';
options.remove_lowfreq = 0;
options.tau = 1;

%% wavelets options
Jmax = log2(n)-1;
options.Jmin = Jmax - 3;
options.wavelet_type = 'daubechies';        % real TF
options.wavelet_type = 'biorthogonal';      % approx TF

D = @callback_atrou;
D = @callback_localdct;
        
%% sparsity inpainting
Mi = feval(U,y,-1,options);
options.x = feval(D, Mi, -1, options);
options.niter = 600;
options.thresh = 'soft';
options.Tmax = .2/5;
options.Tmin = 0;
options.M1 = Mi;
options.D = D;
options.drawiter = 1;
options.x0 = feval(D, M, -1, options); % to monitor error
[xwav,err,lun,Tlist,Err] = perform_iterative_thresholding(U,y,options);
Msol = feval(D, xwav, 1, options);

warning off;
imwrite(clamp(Msol), [rep name '-spars.png'], 'png');
warning on;

dmp.name = name;
dmp.n = n;
dmp.p = P;
dmp.snr_manifold = snr(M,M1);
dmp.snr_sparse = snr(M,Msol);
dump_struct(dmp, [rep 'results.txt'], ['---- ' problem ' ----']);