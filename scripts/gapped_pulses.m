% demonstrate gapped pulse design

Ngap = 100;
Npulse = 10;
x= [-2:.01:2];

rf{1} = [1 1]/2 * pi;

rf{2} = adiabatic_sech(18, 1, 4 * 2*pi/1.87, 4);
rf{2} = 1.5*pi * rf{2} / sum(rf{2});



for n = 1:2
    rf_c{n} = zeros(1,length(rf{n})*Ngap);
    for r = 1:length(rf{n})
        rf_c{n}([1:Ngap] + (r-1)*Ngap) = rf{n}(r)/Ngap;
    end
   rf_c{n} = resample(rf{n},Ngap,1)/Ngap;
    
    mz = ab2inv(abr(rf_c{n},x * length(rf{n})));
    
    rf_gapped{n} = zeros(1,length(rf{n})*Ngap);
    for m= 1:Npulse
        rf_gapped{n}(m+Ngap/2:Ngap:end) = rf{n}/Npulse;
    end
    mz_gapped = ab2inv(abr(rf_gapped{n},x * length(rf{n})));
    
    figure(n)
    subplot(131)
    plot(1:length(rf_c{n}), real(rf_c{n}), 'b--', ...
        1:length(rf_c{n}), imag(rf_c{n}), 'g--', ...
        1:length(rf_gapped{n}), real(rf_gapped{n}), 'b-', ...
        1:length(rf_gapped{n}), imag(rf_gapped{n}), 'g-');
    
    subplot(132)
    plot(mz), axis tight
    subplot(133)
    plot(mz_gapped), axis tight
end
