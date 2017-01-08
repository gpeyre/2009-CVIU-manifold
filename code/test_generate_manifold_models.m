% generate sets of images

n = 150;
nb = 5;
warning off;

rep = 'results/manifold-examples/';
if not(exist(rep))
    mkdir(rep);
end

name = {'fnoise', 'geometrical', 'lic-texture', 'tv-image'};
nameprint = {'regular', 'cartoon', 'oscillating', 'tv'};


name = {'lic-texture'};
nameprint = {'oscillating'};

options.alpha = 3;
options.Jgeometrical = 2;

lic_regularity_list = linspace(5,20, nb);

for k=1:length(name)
    disp(['Processing ' nameprint{k} '.']);
    for i=1:nb
        options.lic_regularity  = lic_regularity_list(i);
        M = load_image(name{k}, n, options);
        imwrite(rescale(M), [rep nameprint{k} '-' num2str(i) '.png'], 'png');
        w = 20; 
        m = 1000;
        p = 3; %  number of patch to save
        H = compute_random_patches(M,w,m);
        H1 = reshape(H, [w^2, m]);
        s = std(H1, 0, 1); [tmp,I] = sort(s, 'descend');
        I = I(1:m/10);
        a = randperm(m/10); I = I(a(1:p));
        for j=1:p
            q = (i-1)*p + j;
            imwrite((H(:,:,I(j))), [rep nameprint{k} '-patch-' num2string_fixeddigit(q,2) '.png'], 'png');
        end
    end
end