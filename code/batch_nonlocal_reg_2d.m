% test for denoising with several parameters

sigma_list = [0.01 0.02 0.05 0.1 0.15 0.2 0.3];
for it=1:length(sigma_list);
    sigma = sigma_list(it);
    test_nonlocal_regularization_2d;
end