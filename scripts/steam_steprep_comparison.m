clear all

vox_thresh_scale = 3;
vox_thresh_num = 2;

mets = {'pyr', 'ala', 'lac'};


% dirpath = '/data/vig1/C13/animal/TRAMP/MS213_101007/';
% imname = 't2fse_ax';
% exppath = {'exp1-steam/', 'exp2-steprep/'};
% specname = {'steam', 'steprep'};
% pol = [31.4 24.3];
% Imets = [224 144 49];
% 
% dirpath = '/data/vig1/C13/animal/Normal_Mice/M036_101201/';
% imname = 'ax';
% exppath = {'', ''};
% specname = {'steam', 'steprep'};
% pol = [19.3 19.3];
% Imets = [234 88 47]; % pyr urea lac
% 
% dirpath = '/data/vig1/C13/animal/UCSF_rats/20101203_1/';
% imname = 'ax';
% exppath = {'', ''};
% specname = {'steam1', 'steprep2'};
% pol = [20 20];
% Imets = [234 149 88 47] -2; % pyr ala urea lac
% Imask = [1:7, 15:18];

dirpath = '/data/vig1/C13/animal/Normal_Mice/M003_110228/';
imname = 'ax';
exppath = {'', ''};
specname = {'steam_phased', 'step_phased'};
pol = [20 20];
Imets = [226 145 53]; % pyr ala lac
Imask = [];

Nmets = length(Imets);


Nexp = length(exppath);

%% read in data

for n= 1:Nexp
    
    overlay_root = sprintf('%s%sspectra/%s',...
        dirpath, exppath{n}, specname{n});
    
    overlay_image{n}=read_ddf_image(overlay_root);
    
    peaks{n} = getpeaks(overlay_image{n}.img / pol(n), Imets);
%     [tags{n} bits_used] = read_tags(sprintf('%s%sspectra/%s_raw_tags.txt',dirpath, exppath{n}, specname{n}), ...
%         prod(overlay_image{n}.ddf.npix));
end

% Ntags = floor(log2(bits_used))+1;



%% data analysis

vox_thresh = max(peaks{vox_thresh_num}(1,1,1,:)) * vox_thresh_scale;
vox_mask = peaks{1} > vox_thresh & peaks{2} > vox_thresh;
vox_mask(:,:,Imask,:) = 0; % remove cath(?)


for m = 1:Nmets
    
    for n=1:Nexp
        temp{n} = peaks{n}(:,:,:,m);
        mets_all{n,m} =  temp{n}(find(vox_mask(:,:,:,m)));
        
        temp_full{n} = peaks{n}(:,:,:,m);
%        for t = 1:Ntags
%             mets_tags(n,m,t) = mean(temp_full{n}(find(tags{n} == 2^(t-1))));
%         end
%         
        mets_avg(n,m) = mean(mets_all{n,m});
        
    end
    
    
end


%%
I1 = 2;  I2 = 1;

mean(mets_all{I1,1} ./mets_all{I2,1})
mean(mets_all{I1,2} ./mets_all{I2,2})
mean(mets_all{I1,3} ./mets_all{I2,3})
%mean(mets_all{I1,4} ./mets_all{I2,4})
mean([mets_all{I1,1};mets_all{I1,2};mets_all{I1,3}] ./ [mets_all{I2,1};mets_all{I2,2};mets_all{I2,3}])

figure(1)
plot([mets_all{1,1};mets_all{1,2};mets_all{1,3}], [mets_all{2,1};mets_all{2,2};mets_all{2,3}], 'x')
plot([mets_all{1,1};mets_all{1,2};mets_all{1,3}], [mets_all{2,1};mets_all{2,2};mets_all{2,3}], 'x')
figure(2)
peaks_temp1 = peaks{1}(:,:,:,[1,2 3]); peaks_temp2 = peaks{2}(:,:,:,[1,2 3]); plot(peaks_temp1(:), peaks_temp2(:), 'x')
%axis equal