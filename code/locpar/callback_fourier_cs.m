function y = callback_fourier_cs(x,dir,options)

mask = getoptions(options, 'mask', 1, 1);
n = size(x,1);
if dir==1
   y = fft2(x) / n;
   y(mask==1) = 0;    
else
   x(mask==1) = 0;    
   y = real(ifft2(x)) * n;
end