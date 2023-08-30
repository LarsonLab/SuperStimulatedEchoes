function [data, spec, phi0, Amets] = process_madsteam(fname, Iref, phi1, mets, Imets, TR, do_plots);

if nargin < 7
    do_plots = 0;
end

data = rawloadX(fname);
[Nspec, Nt] = size(data);

% 
% %[spec_correct phi_correct spec_orig M] = phase_est_twopeaks('exp1-MADsteam_TE1/steam_te1', Imets([2,4]));
% %[spec_correct phi_correct spec_orig M] = phase_est_twopeaks('exp2-MADsteam_TE2/steam_te2', Imets([2,4]));
% 
% phi1 = .102;  %.1  OR SKIP 8 samples!
phi1mat = repmat( exp(1i * (1:256).' * phi1), [1 Nt]);

phi0width = 3;
Iphi0 = Iref + [-phi0width:phi0width];

df = 2500/Nspec;
f = [1:Nspec] *df;

win = spec_apod_win(Nspec, 1/2500, 10, 1).';

spec = fftc(data.* repmat(win, [1 Nt])) .* phi1mat;

for t = 1:Nt
    phi0(t) = -angle(sum(spec(Iphi0,t)));

    spec(:,t) = spec(:,t) * exp(1i * phi0(t));
    
    if do_plots
        figure(99)
        cplot(f, spec(:,t))
        pause(.5)
    end
end


Nwidth = 5;
t = (0:Nt-1) * TR;

Amets = zeros(length(Imets), Nt);
for n =1:length(Imets)
    Amets(n,:) = sum(spec(Imets(n) + (-Nwidth:Nwidth), :), 1);
    figure(n)
    plot(t,real(Amets(n,:)),'x--',t,imag(Amets(n,:)),'o--', ...
        t, abs(Amets(n,:)), ':');
    xlabel('time (s)'), legend('real', 'imag')
    title(mets{n})
end
