function [spec, data] = recon_madsteam_1Depsi(fname, outname, Iurea);
% recon phase sensitive 1D MRSI
addpath /netopt/share/lib/local/brain/matlab/

%%fb epsi
%fname = '/data/vig1/C13/animal/TRAMP/MS243_110421/experiment2-madsteam_epsi/madsteam_epsi';
Iref = Iurea +[-2:2];
%Iref = 16:20 ;

ddf_info = readddf_file(['spectra/' fname]);
epsi_params = get_epsi_info(fname);

epsi_params.Nlobes = 59;  % Best
epsi_params.dt = 1/25000;
epsi_params.Nskip  = 27;
epsi_params.Nepsi = 16;
hz_apod = 10;

[pfile_data, header, rhuser] = rawloadX(fname, [], 1, 1, 1);
Nt = size(pfile_data,2);

w = spec_apod_win(epsi_params.Nlobes, epsi_params.Nepsi*epsi_params.dt*1e-6, hz_apod,1);
% TEST PSF!!pfile_data(:) = 1;

%pfile_data = sum(pfile_data,2);


%%
Nstart = 0; % looks goodish
data = reshape(pfile_data(Nstart+(1:(epsi_params.Nskip+epsi_params.Nepsi)*epsi_params.Nlobes),:), [epsi_params.Nskip+epsi_params.Nepsi epsi_params.Nlobes Nt]);
data = data(1:epsi_params.Nepsi, :,:) .* repmat(w, [epsi_params.Nepsi 1 Nt]);

%tshift ~ -14 (epsi_params.Nlobes=58,59)
for tshift = 0
    data_kx_f = fft(data, [], 2);
    data_kx_f = rephase_epsi(data_kx_f, 1, 2, (epsi_params.Nskip+epsi_params.Nepsi)*epsi_params.Nlobes, tshift);
    
    data_reph = ifft(data_kx_f, [], 2);
    
    spec = zeros(size(data_reph));
    for t = 1:Nt
        spec(:,:,t) = fftnc(data_reph(:,:,t));
    end
    
    if 1
        for x = 1:size(spec,1)
            for t = 1:size(spec,3)
                phi0 = find_phase_corr(spec(x,Iref,t));
                spec(x,:,t) = spec(x,:,t) * exp(1i*phi0);
            end
        end
    else
        phi0= find_phase_corr(spec);
        spec = spec * exp(1i*phi0);
    end
    
    %spec(1, 1:size(spec,1), 1:size(spec,2)) = spec; spec = spec(1,1:size(spec,1),:);
    if 0
        specplot = permute(spec, [1 3 2]);
        %    cplot_voxels(spec(Nt/2+1 +[-2:2], epsi_params.Nepsi/2+1 +[-2:2],:))
        cplot_voxels(specplot,[],[],[],[min(imag(spec(:))) max(real(spec(:)))]);
        %    pause
    end
end


%% peaks
Imets = [-5 0 34]+ Iurea;
mets = {'lac', 'urea', 'pyr'};
Nwidth = 3;

Amets = zeros(length(Imets), size(spec,1), Nt);
for n =1:length(Imets)
    Amets(n,:,:) = sum(spec(:,Imets(n) + (-Nwidth:Nwidth), :), 2);
if 1
    t = [0:Nt-1];       
    figure(n)
    plot(t,real(squeeze(Amets(n,9,:))),'x--',t,imag(squeeze(Amets(n,9,:))),'o--', ...
        t, abs(squeeze(Amets(n,9,:))), ':');
    xlabel('time (s)'), legend('real', 'imag')
    title(mets{n})
    if 0
        print('-depsc', [expdir fname '_' mets{n}]);
    end
end
end


%% write

ddf_info.npix(epsi_params.epsiaxis) = epsi_params.Nepsi;
ddf_info.pixel_spacing(epsi_params.epsiaxis) = epsi_params.epsires;
ddf_info.specpoints = epsi_params.Nlobes;

ddf_info.numdim = 4;
ddf_info.root_name = fname;

ddf_info.center(epsi_params.epsiaxis) = ...
    ddf_info.center(epsi_params.epsiaxis) - 0.5 .* epsi_params.epsires;
ddf_info.toplc(epsi_params.epsiaxis)=...
    ddf_info.center(epsi_params.epsiaxis) + ((epsi_params.Nepsi-1)/2)*epsi_params.epsires;
% ddf_info.toplc(1)=ddf_info.center(1) - ((N(1)-1)/2)*res(1);
% ddf_info.toplc(2)=ddf_info.center(2) - ((N(2)-1)/2)*res(2);
% ddf_info.toplc(3)=ddf_info.center(3) + ((N(3)-1)/2)*res(3);
ddf_info.box_center = ddf_info.center;
ddf_info.box_size(epsi_params.epsiaxis) = epsi_params.Nepsi*epsi_params.epsires;
ddf_info.acq_toplc = ddf_info.toplc;
ddf_info.acq_spacing = ddf_info.pixel_spacing;
ddf_info.acq_n_data_points = ddf_info.specpoints;
ddf_info.acq_n_points = ddf_info.npix;

outdata.ddf = ddf_info;
for t = 1:Nt
    outdata.img = flipdim(permute(spec(:,:,t), [2 1]),2);
    write_ddf_image(sprintf('spectra/%s%02d', outname, t), outdata);
end

outdata.ddf.specpoints = Nt;
for n =1:length(mets)
    outdata.img = flipdim(permute(Amets(n,:,:), [3 2 1]),2);
    write_ddf_image(sprintf('spectra/%s_%s', outname, mets{n}), outdata);
end



