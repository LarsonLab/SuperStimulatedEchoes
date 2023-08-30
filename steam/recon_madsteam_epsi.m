% recon phase sensitive MRSI

%%fb epsi
fname = '/data/vig1/C13/animal/TRAMP/MS227_110202/madsteam_cor';  Iref = 15:18;
fname = '/data/vig1/C13/animal/Normal_Mice/M038_101222/madsteam1_raw';  Iref = 12:15;

Nlobes = 59;  % Best
dt = 1/25000;
Nskip  = 27;
Nread = 16;
hz_apod = 10;

[pfile_data, header, rhuser] = rawloadX(fname, [], 1, 1, 1);
Npe = size(pfile_data,2);

w = spec_apod_win(Nlobes, Nread*dt*1e-6, hz_apod,1);
% TEST PSF!!pfile_data(:) = 1;

%pfile_data = sum(pfile_data,2);


%%
Nstart = 0; % looks goodish
data = reshape(pfile_data(Nstart+(1:(Nskip+Nread)*Nlobes),:), [Nskip+Nread Nlobes Npe]);
data = data(1:Nread, :,:) .* repmat(w, [Nread 1 Npe]);

%tshift ~ -14 (nlobes=58,59)
for tshift = 0
    data_kx_f = fft(data, [], 2);
    data_kx_f = rephase_epsi(data_kx_f, 1, 2, (Nskip+Nread)*Nlobes, tshift);
    
    data_reph = ifft(data_kx_f, [], 2);
    
    spec = fftnc(data_reph);
    
    if 1
        for x = 1:size(spec,1)
            for y = 1:size(spec,3)
                phi0 = find_phase_corr(spec(x,Iref,y));
                spec(x,:,y) = spec(x,:,y) * exp(1i*phi0);
            end
        end
    else
        phi0= find_phase_corr(spec);
        spec = spec * exp(1i*phi0);
    end
    
    %spec(1, 1:size(spec,1), 1:size(spec,2)) = spec; spec = spec(1,1:size(spec,1),:);
    spec = permute(spec, [1 3 2]);
%    cplot_voxels(spec(Npe/2+1 +[-2:2], Nread/2+1 +[-2:2],:))
cplot_voxels(spec,[],[],[],[min(imag(spec(:))) max(real(spec(:)))]);
    pause
end

return

%% kpe, kepsi, kf raw data
for Nstart = 0:size(pfile_data,1)-(Nskip+Nread)*Nlobes-1
    data = reshape(pfile_data(Nstart+(1:(Nskip+Nread)*Nlobes),:), [Nskip+Nread Nlobes Npe]);
    data = data(1:Nread, :,:) .* repmat(w, [Nread 1 Npe]);
    
    data_kx_f = fft(data, [], 2);
    data_kx_f = rephase_epsi(data_kx_f, 1, 2, (Nskip+Nread)*Nlobes);
    
    data_reph = ifft(data_kx_f, [], 2);
    
    spec = fftnc(data_reph);
    
    for x = 1:size(spec,1)
        for y = 1:size(spec,3)
            phi0 = find_phase_corr(spec(x,Iref,y));
            spec(x,:,y) = spec(x,:,y) * exp(1i*phi0);
        end
    end
    
    spec = permute(spec, [3 1 2]);
    
    cplot_voxels(spec(Npe/2+1 +[-2:2], Nread/2+1 +[-4:4],:))%cplot_voxels(spec);
    real_sum(Nstart+1) = sum(real(spec(:)* exp(1i * phi0)));
    abs_sum(Nstart+1) = sum(abs(spec(:)* exp(1i * phi0)));
    pause
end
