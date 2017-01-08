% batch for non-local zooming


pbm = 'zoom';
pbm = 'inpainting';
pbm = 'cs';

% 'barb',
namelist = { 'peppers','lena' };  %'boat',
methlist = {'tv'};
methlist = {'sparsity-wavelets', 'nl-means'};
methlist = {'nl-means'};


switch pbm
    case 'cs'
        %%%%%% COMPRESSED SENSING %%%%%
        slist = [8]; %  12 4 16
        N = 256^2;

        for is = 1:length(slist)
            for iname = 1:length(namelist)
                for imeth = 1:length(methlist)
                    name = namelist{iname};
                    p = round(N/slist(is));
                    method = methlist{imeth};
                    test_inverse_problem;
                end
            end
        end

    case 'zoom'
        %%%%%% ZOOMING %%%%%
        slist = [8 4 2];
        for iname = 1:length(namelist)
            for is = 1:length(slist)
                for imeth = 1:length(methlist)
                    name = namelist{iname};
                    s = slist(is);
                    method = methlist{imeth};
                    test_inverse_problem;
                end
            end
        end

    case 'inpainting'
        nmiss_list = [.7 .8 .9 .95];
        for imiss = 1:length(nmiss_list)
            for iname = 1:length(namelist)
                for imeth = 1:length(methlist)
                    name = namelist{iname};
                    nmiss = nmiss_list(imiss);
                    method = methlist{imeth};
                    test_inverse_problem;
                end
            end
        end
end