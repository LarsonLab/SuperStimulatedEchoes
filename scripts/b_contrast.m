clear all

vox_thresh_scale = 3;
vox_thresh_num = 1;

mets = {'pyr', 'ala', 'lac'};
Nmets = 3;


dataset = 2;
switch dataset
    % Normal
    case 1
        dirpath = '/data/vig1/C13/animal/Normal_Mice/FVB_Normal_100421/';
        imname = 't2fse_ax';
        exppath = {'experiment2-3d/', 'experiment3-steprep_TM5ms/', 'experiment1-steprep/'};
        specname = {'3d', '3dsteprep2', '3dsteprep'};
        pol = [25.8 30.5 33.2];
        Imets = [228 148 55];
        
    case 2
        dirpath = '/data/vig1/C13/animal/TRAMP/MS223_100729/';
        imname = 't2fse_ax';
        exppath = {'experiment2-steprep_tm10ms/', 'experiment1-steprep_tm1s/'};
        specname = {'steprep2', 'steprep1'};
        pol = [20 20];
        Imets = [224 144 49];
end


Nexp = length(exppath);

%% read in data

for n= 1:Nexp
    
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

for m = 1:Nmets
    
    for n=1:Nexp
        temp{n} = peaks{n}(:,:,:,m);
        mets_all{n,m} =  temp{n}(find(vox_mask(:,:,:,m)));
        for t = 1:Ntags
            mets_tags(n,m,t) = mean(temp{n}(find(tags == 2^(t-1))));
        end
        
        mets_avg(n,m) = mean(mets_all{n,m});
        
    end
    
    
end


%%
I1 = 1;  I2 = 2;

mean(mets_all{I1,1} ./mets_all{I2,1})
mean(mets_all{I1,2} ./mets_all{I2,2})
mean(mets_all{I1,3} ./mets_all{I2,3})
mean([mets_all{I1,1};mets_all{I1,2};mets_all{I1,3}] ./ [mets_all{I2,1};mets_all{I2,2};mets_all{I2,3}])


plot([mets_all{I1,1};mets_all{I1,2};mets_all{I1,3}], [mets_all{I2,1};mets_all{I2,2};mets_all{I2,3}], 'x')
plot(peaks{1}(:), peaks{2}(:), 'x')
axis equal