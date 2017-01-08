function B = compute_circle_mean(B1,B2,w1,w2, s)

% compute_circle_mean - average on circle
%
% B = compute_circle_mean(B1,B2,w1,w2, s);
%
%   B1,B2 and B3 are assumed to be angle, ie. to be always
%   computed modulo s, with detaulf s=2*pi.
%
%   Copyright (c) 2006 Gabriel Peyré

if nargin<3
    w1 = B1*0+1/2;
end
if nargin<4
    w2 = B1*0+1/2;
end
if nargin<5
    s = 2*pi;
end

if size(w1,1)==1 && size(w1,2)==1 && size(w2,1)==1 && size(w2,2)==1
    w1 = B1*0 + w1;
    w2 = B1*0 + w2;
end

% weighted circle mean
b1 = exp( 2i*pi / s * B1 );
b2 = exp( 2i*pi / s * B2 );
b = b1.*w1 + b2.*w2;

I = find( abs(b(:))<1e-5 );
B = B1;
B(I) = mod( w1(I).*B1(I) + w2(I).*B2(I), s );
I = find( abs(b(:))>=1e-5 );
B(I) = mod(angle(b(I)), s);
