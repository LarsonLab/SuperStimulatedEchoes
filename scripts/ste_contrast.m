clear all

dataset = 202;
setup_dataset

%% read in data

for n= 1:2
    
    overlay_root = sprintf('%s%sspectra_L1_TV0.0001_xfm0.0005_iter5x20/%s_phased',...
        dirpath, exppath{n}, specname{n});
    
    overlay_image{n}=read_ddf_image(overlay_root);
    
    peaks{n} = getpeaks(overlay_image{n}.img / pol(n), Imets);
end

[tags bits_used] = read_tags(sprintf('%s%sspectra_L1_TV0.0001_xfm0.0005_iter5x20/%s_sampled_tags.txt',dirpath, exppath{1}, specname{1}), ...
    prod(overlay_image{1}.ddf.npix));
Ntags = floor(log2(bits_used))+1;

%% data analysis

vox_thresh = max(peaks{1}(1,1,1,:)) * vox_thresh_scale;
if length(vox_thresh_scale) == Nmets
    Sp = size(peaks{1});
    vox_mask = peaks{1} > repmat(reshape(vox_thresh,[1 1 1 Nmets]), [Sp(1:3) 1]);
else
    vox_mask = peaks{1} > vox_thresh;
end

ste_ratio = peaks{2} ./ peaks{1};

for m = 1:Nmets

    for n=1:2
        temp{n} = peaks{n}(:,:,:,m);
        mets_all{n,m} =  temp{n}(find(vox_mask(:,:,:,m)));
    end
    
    temp_ratio = ste_ratio(:,:,:,m);
    met_ratios{m} = temp_ratio(find(vox_mask(:,:,:,m)));
    
    met_ratio_avg(m) = mean(met_ratios{m});
    met_ratio_std(m) = std(met_ratios{m});
    
    for t = 1:Ntags
        ratio_tags(m,t) = mean(temp_ratio(find(tags == 2^(t-1))));
        for n=1:2
            mets_tags(n,m,t) = mean(temp{n}(find(tags == 2^(t-1))));
        end
    end

end

