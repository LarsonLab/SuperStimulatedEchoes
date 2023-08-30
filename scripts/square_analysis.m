%% setup global parameters
clear all

mets = {'lac','ala','pyr'};

thresh = 0;

extra_scale = 1;
vox_thresh_scale = 2.5;

%% choose dataset for overlay
dataset = 1;

switch dataset
    case 1
        % square encoding
        dirpath = '/data/vig1/C13/animal/UCSF_rats/20090428_1/';
        imname = 't2fse_ax';
        exppath = {'experiment2-nosteam/', 'experiment3-hardste/','experiment1-square2/'};
        specname = {'spec_nosteam', 'spec_hardste', 'steam_square2'};
        pol = [25 25 25];
        
        Imets = [99 191 21];
        I_liver = [3 4 5 3 4 5;
            5 4 4 4 4 4;
            10 10 10 11 11 11];
        I_kidney = [3 3 7;
            5 5 5;
            12 13 13];
        
    case 2
        % square encoding
        dirpath = '/data/vig1/C13/animal/UCSF_rats/20090805_1/';
        imname = 't2fse_ax';
        exppath = {'experiment2-3d/', 'experiment3-steam/','experiment1-steprep/'};
        specname = {'csi_3d', 'steamcsi', 'csi_steprep'};
        pol = [16.6 16.6 15];
        
        Imets = [61 155 237];
        %dummies
        I_liver = [3 4 5 3 4 5;
            5 4 4 4 4 4;
            10 10 10 11 11 11];
        I_kidney = [3 3 7;
            5 5 5;
            12 13 13];
        
end


%% read in data

for n= 1:3
    
    overlay_root = sprintf('%s%sspectra/%s',...
        dirpath, exppath{n}, specname{n});
    
    overlay_image{n}=read_ddf_image(overlay_root);
    
    peaks{n} = getpeaks(overlay_image{n}.img / pol(n), Imets);
end

I1_liver = (I_liver(3,:)-1)*size(peaks{1},1)*size(peaks{1},2) + ...
    (I_liver(2,:)-1)*size(peaks{1},1) + I_liver(1,:);
I1_kidney = (I_kidney(3,:)-1)*size(peaks{1},1)*size(peaks{1},2) + ...
    (I_kidney(2,:)-1)*size(peaks{1},1) + I_kidney(1,:);


%% data analysis

vox_thresh = max(peaks{1}(1,1,1,:)) * vox_thresh_scale;
if length(vox_thresh_scale) ==3
    Sp = size(peaks{1});
    vox_mask = peaks{1} > repmat(reshape(vox_thresh,[1 1 1 3]), [Sp(1:3) 1]) & ...
        peaks{2} > repmat(reshape(vox_thresh,[1 1 1 3]), [Sp(1:3) 1]) & ...
        peaks{3} > repmat(reshape(vox_thresh,[1 1 1 3]), [Sp(1:3) 1]);
else
    vox_mask = peaks{1} > vox_thresh & peaks{2} > vox_thresh & peaks{3} > vox_thresh;
end

for n = 1:2
    
    square_ratio = peaks{n+1} ./ peaks{1}  .* vox_mask;
    
    for m = 1:3
        temp = square_ratio(:,:,:,m);
        met_ratios{n,m} = temp(find(vox_mask(:,:,:,m)));
        met_ratio_avg(n,m) = mean(met_ratios{n,m});
        met_ratio_std(n,m) = std(met_ratios{n,m});
        ratio_liver(n,m) = mean(temp(I1_liver));
        ratio_kidney(n,m) = mean(temp(I1_kidney));
    end
    ratio_avg(n) = mean([met_ratios{n,1}; met_ratios{n,2}; met_ratios{n,3}]);
    
end

for n=1:3
    for m=1:3
        temp = peaks{n}(:,:,:,m);
        mets_all{n,m} =  temp(find(vox_mask(:,:,:,m)));
    end
end

%plot(peaks{2}(:), peaks{3}(:),'x'), xlabel('STE'), ylabel('sSTE')
figure(1)
plot([mets_all{2,1};mets_all{2,2};mets_all{2,3}], [mets_all{3,1};mets_all{3,2};mets_all{3,3}], 'x')
plot([mets_all{2,1};mets_all{2,2};mets_all{2,3}], [mets_all{3,1};mets_all{3,2};mets_all{3,3}], 'x')
figure(2)
peaks_temp1 = peaks{2}(:,:,:,[1,2 3]); peaks_temp2 = peaks{3}(:,:,:,[1,2 3]); plot(peaks_temp1(:), peaks_temp2(:), 'x')

