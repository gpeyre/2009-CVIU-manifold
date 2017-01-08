% test for patch creation

name = 'oriented_grid';
name = 'edges';
name = 'disks';
options.dim = [1 1] * w;
options.smoothing = 0.02;
options.nbr = 500;
if strcmp(name, 'oriented_grid')
    options.nbr = 5*w; 
end
options.sampling = 'unif';
options.radius = 0.25;
options.spacing = 0.5;
D{Jmax} = load_images_dataset(name, options);