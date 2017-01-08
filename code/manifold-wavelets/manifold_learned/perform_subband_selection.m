function [selx,sely] = perform_subband_selection(j,s,M)

% perform_subband_selection - select a subset of wavelet coefficients.
%
%   [selx,sely] = perform_subband_selection(j,s);
% or
%   M1 = perform_subband_selection(j,s,M);
%
%   j is the scale, s is the orientation (1 or 2).
%
%   Copyright (c) 2006 Gabriel Peyré

if s==0
    selx = 1:2^j;
    sely = 1:2^j;    
elseif s==1
    selx = 1:2^(j+1);
    sely = 2^j+1:2^(j+1);
else
    selx = 2^j+1:2^(j+1);
    sely = 1:2^j;
end

if mod(j,2)==0
    a = selx; 
    selx = sely; sely = a;
end

if nargin==3
    selx = M(selx,sely);
end