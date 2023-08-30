%dataset = 201;
dataset = 103;
setup_dataset

do_overlays = 1;
write_pngs = 1;

thresh = -1;

base_low = 0;
base_high = 1.2e4;

overlay_low = 0;
overlay_high_scale = [.7 .9];

figsize = [150 480];
sep = .01;
imname = 't2fse_cor';

alpha = 0.4;
colorbar_on = 0;        

%% read in data

if do_overlays
    bottom_image_root = [dirpath 'images/' imname];
    bottom_image=read_idf_image(bottom_image_root);
end

for n= 1:2
    
    overlay_root = sprintf('%s%sspectra_L1_TV0.0001_xfm0.0005_iter5x20/%s_phased',...
        dirpath, exppath{n}, specname{n});
    
    overlay_image{n}=read_ddf_image(overlay_root);
    
    peaks{n} = getpeaks(overlay_image{n}.img / pol(n), Imets);
end

if do_overlays
    
    old_image.idf = ddf_to_idf(overlay_image{1}.ddf);  % this function is a little bit of a hack
    old_image.idf.center = old_image.idf.center * old_image.idf.dcos * bottom_image.idf.dcos';
    old_image.idf.fov = old_image.idf.fov * old_image.idf.dcos * bottom_image.idf.dcos';
    old_image.idf.pixelsize = old_image.idf.pixelsize * old_image.idf.dcos * bottom_image.idf.dcos';
    old_image.idf.npix = old_image.idf.npix * old_image.idf.dcos * bottom_image.idf.dcos';
    old_image.idf.toplc = old_image.idf.toplc * old_image.idf.dcos * bottom_image.idf.dcos';
    
    [xtemp ytemp overlay_zlocs] = find_voxel_centers(old_image.idf);
    [xtemp ytemp bottom_zlocs] = find_voxel_centers(bottom_image.idf);
    Islices = find( (overlay_zlocs >= min(bottom_zlocs)) & (overlay_zlocs <= max(bottom_zlocs)));
    for s= 1:length(Islices)
        [temp slices(s)] = min(abs(bottom_zlocs - overlay_zlocs(Islices(s))));
    end
    
end

if exist('Idisp') && ~isempty(Idisp)
    slices = slices(Idisp);
end


%% data analysis

vox_thresh = max(peaks{1}(1,1,1,:)) * vox_thresh_scale;
if length(vox_thresh_scale) == Nmets
    Sp = size(peaks{1});
    vox_mask = peaks{1} > repmat(reshape(vox_thresh,[1 1 1 Nmets]), [Sp(1:3) 1]);
else
    vox_mask = peaks{1} > vox_thresh;
end


%%
if exist('overlay_peaks_steam') && (length(overlay_peaks_steam) == length(Imets))
    do_interp = 0;
else
    do_interp = 1;
end

    overlay_peaks_nosteam_high = max(peaks{1}(:)) *overlay_high_scale(1);
    overlay_peaks_steam_high = max(peaks{2}(:)) *overlay_high_scale(2);
for m = 1:length(Imets)
    peaks_steam = peaks{2}(:,:,:,m);
    peaks_nosteam = peaks{1}(:,:,:,m);
    
    Imask = find(vox_mask(:,:,:,m));
    
    if do_overlays
        if do_interp
            old_image.img = permute(peaks_nosteam, [1 3 2]);
            overlay_peaks_nosteam{m} = resample_m(old_image, bottom_image.idf, [], 'cubic');
            
            old_image.img = permute(peaks_steam, [1 3 2]);
            overlay_peaks_steam{m} = resample_m(old_image, bottom_image.idf, [], 'cubic');
        end
        
        x_bottom = 1:256; y_bottom = 1:256;
        
        for s = slices
            if write_pngs
                fs = figure(1);
                set(fs,'Color', 'white', 'Position', [1 1  figsize]);
                figure(1);
            else
                fs = figure((m-1)*max(slices)+s);
                set(fs,'Color', 'white', 'Position', [50*find(s==slices)   50   figsize]);
                figure((m-1)*max(slices)+s);
            end
                        
            subplot('Position', [sep, sep, 1-2*sep, (1-2*sep)/2])
            color_overlay_coords(bottom_image.img(:,:,s), overlay_peaks_steam{m}.img(:,:,s), ...
                x_bottom, y_bottom, x_bottom, y_bottom, base_low, base_high, ...
                overlay_low, overlay_peaks_steam_high, thresh, alpha, colorbar_on);
            axis(axzoom), view([90 90]);
            axis off;
            
            subplot('Position', [sep, sep+1/2, 1-2*sep, (1-2*sep)/2])
            color_overlay_coords(bottom_image.img(:,:,s), overlay_peaks_nosteam{m}.img(:,:,s), ...
                x_bottom, y_bottom, x_bottom, y_bottom, base_low, base_high, ...
                overlay_low, overlay_peaks_nosteam_high, thresh, alpha, colorbar_on);
            
            axis(axzoom), view([90 90]);
            axis off;
            
            if write_pngs
                
                print('-dpng', sprintf('%ssteamcompare_%s_cor_slice%d.png', dirpath, mets{m}, s))
                
            end
        end
    end

    
end

