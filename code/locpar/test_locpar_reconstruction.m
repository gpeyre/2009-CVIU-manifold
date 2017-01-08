% test for manifold of locally parallel textures

n0 = 140;
n = 64*2;
n = 96;
name = 'fingerprint3';

rep = ['results/locpar/'];
if not(exist(rep))
    mkdir(rep);
end

M = load_image(name, n0);
M = rescale( crop(M,n) );

w = 12;

% parameter of the dictionary
options.fmin = 3;
options.fmax = 6;
options.nfreq = 16;
options.nphase = 14;
options.ntheta = 18;
options.gamma = .35;
D = compute_locpar_dictionary(w,options);

%% perfrom analysis-reconstruction
clf;
sublist = [12 4 2];
M1 = {}; lgd = {}; filename = {};
for i=1:length(sublist)
    options.sub = sublist(i);
    M1{end+1} = perform_patch_manifold_projection(M,D,options);
    lgd{end+1} = ['SNR=' num2str(snr(M,M1{end}))];
    filename{end+1} = ['w' num2str(sublist(i))];
end
imageplot(M1, lgd, 2,2);

M1{end+1} = M;
filename{end+1} = 'original';

options.clamp = 1;
options.base_str = [rep name '-'];
save_image(M1,filename,options);

%% perform analysis for various gamma
gammalist = .05:.05:1;
options.sub = 2;
w = 8;
M1 = {}; lgd = {}; filename = {}; err = [];
for i=1:length(gammalist)
    progressbar(i, length(gammalist));
    gamma = gammalist(i);
    options.gamma = gamma;
    D = compute_locpar_dictionary(w,options);
    M1{end+1} = perform_patch_manifold_projection(M,D,options);
    err(end+1) = snr(M,M1{end});
    lgd{end+1} = ['gamma=' num2str(gamma) ', SNR=' num2str(snr(M,M1{end}))];
    filename{end+1} = ['gam' num2str(100*gamma)];
end
imageplot(M1, lgd);
save_image(M1,filename,options);
plot(gammalist, err); axis tight;