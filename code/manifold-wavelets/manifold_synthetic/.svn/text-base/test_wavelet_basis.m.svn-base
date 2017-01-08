n = 64;
Jmin = 2;
s = 2;

options.sigma = 0.05;
options.w = 9;
options.n_delta = 9;
options.n_theta = 12;

Theta = zeros(n);
Delta = zeros(n);

[selx,sely] = perform_subband_selection(Jmin,s);
Delta(selx(2),sely(2)) = 0;
Theta(selx(2),sely(2)) = 1;

ParamW.Delta = Delta;
ParamW.Theta = Theta;

Param = perform_orientation_transform(ParamW,1,-1,options);