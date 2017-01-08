% test different parameters for denoising


for use_normalized=0:1
    locality = 'nonlocal';
    test_nonlocal_regularization_2d;
    locality = 'semilocal';
    test_nonlocal_regularization_2d;
end

return;


sigma_sp_list = linspace(2,10,10);
sigma_va_list = linspace(0.03,0.3, 10);
locality = 'semilocal';
nsp = length(sigma_sp_list);
nva = length(sigma_va_list);

inb = 0;
Err = zeros(nsp,nva);
for isp = 1:nsp
    for iva = 1:nva
        inb = inb + 1;
        sigma_sp = sigma_sp_list(isp);
        sigma_va = sigma_va_list(iva);
        disp( ['####### test ' num2str(inb) '/' num2str(nsp*nva) ' #######'] );
        test_nonlocal_regularization_2d;
        Err(isp,iva) = pgauss;
    end
end

lambda = sigma_sp_list./sigma_va_list;
clf;
imageplot(sigma_va_list, sigma_sp_list, Err); axis image;
xlabel('\si');
ylabel('\la');
