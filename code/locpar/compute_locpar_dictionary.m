function D = compute_locpar_dictionary(w,options)

% parameter of the dictionary
nfreq = getoptions(options, 'nfreq', 16);
nphase = getoptions(options, 'nphase', 14);
ntheta = getoptions(options, 'ntheta', 18);
fmin = getoptions(options, 'fmin', 3);
fmax = getoptions(options, 'fmax', 6);
gamma = getoptions(options, 'gamma', 1);

freq = linspace(fmin,fmax,nfreq);
phase = linspace(0,1,nphase+1); phase(end) = [];
theta = linspace(0,pi,ntheta+1); theta(end) = [];
[X,Y,Phase,Theta,Freq] = ndgrid(1:w,1:w,phase,theta,freq);
D = cos( 2*pi./Freq .* ( X.*sin(Theta) - Y.*cos(Theta) ) - 2*pi*Phase );
D = abs(D).^gamma .* sign(D);

% normalize
D = D ./ repmat( sqrt(sum(sum(D.^2))), [w w] );
m = nfreq*nphase*ntheta;
D = reshape(D,[w w m]);
D1 = reshape(D,[w*w m]);
