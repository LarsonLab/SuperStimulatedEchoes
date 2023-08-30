% To test STE pulse design:

% using some parameters from successful C13 experiments

GAMMA = 1070.5;

z = linspace(-5, 5, 201);  % mm
f = linspace(-500, 500, 201);  % Hz

Ivox = find(abs(z) < 2.5);

% want fshift*TE/2 = integers, for in-phase after conversion
% fshift = 385 (lac), 193 (ala)
% multiples of 5.2ms are good
TE = 31.2e-3;
dt = 10*4e-6;

% conventional
n = 1;
a_gencode = 3.5 / 10;  % G/mm
pw_gencode = 4e-3;
pw_rf = 264e-6;
peak_B1 = pi/2 / (2*pi*pw_rf * GAMMA);

Nrf = round(pw_rf/dt);
Nrf_hard = Nrf;
Ngencode = round(pw_gencode/dt);
Nwait = round((TE/2 - pw_rf - pw_gencode)/dt);

rf{n} = [pi/2*ones(1,Nrf)/Nrf, zeros(1, Ngencode+Nwait), pi/2*ones(1,Nrf)/Nrf];
g{n} = [zeros(1,Nrf), ones(1,Ngencode)*a_gencode, zeros(1,Nwait+Nrf)];

rfex{n} = [pi/2*ones(1,Nrf)/Nrf, zeros(1, Ngencode+Nwait)];
gex{n} = [zeros(1,Nrf), ones(1,Ngencode)*a_gencode, zeros(1,Nwait)];

% some linear phase accrual in frequency after final RF pulse...
% not dependant on TE, but on RF duration


% frequency selective encoding
n = n+1;
pw_rf = 10e-3;
bw = 250;
Nrf = round(pw_rf/dt);
rf_pyr = dzrf(Nrf, pw_rf*bw, 'ex');
Nwait = round((TE/2 - pw_rf - pw_gencode)/dt);

rf{n} = [rf_pyr, zeros(1, Ngencode+Nwait), rf_pyr(end:-1:1)];
g{n} = [zeros(1,Nrf), ones(1,Ngencode)*a_gencode, zeros(1,Nwait+Nrf)];

Nrefocus = round((TE/2 - pw_gencode)/dt - Nrf_hard/2);
rfex{n} = [pi/2*ones(1,Nrf_hard)/Nrf_hard, zeros(1, Ngencode+Nrefocus)];
gex{n} = [zeros(1,Nrf_hard), ones(1,Ngencode)*a_gencode, zeros(1,Nrefocus)];

fshift = 385; % 192.6; 385.6;

for n = 1:length(rf)
    mz{n} = ab2inv(abr(rf{n}, 2*pi*dt*(GAMMA*g{n} + i*ones(size(g{n}))), ...
        z, f));
    mxytemp = ab2ex(abr(rfex{n}, 2*pi*dt*(GAMMA*gex{n} + i*ones(size(gex{n}))), ...
        z, f));
    mxy{n} = -mz{n} .* mxytemp;
    
    fmod = 1; %exp(i*2*pi*fshift*[0:length(rfex{n})-1]*dt);
    mxytemp_shift = ab2ex(abr(rfex{n} .* fmod, 2*pi*dt*(GAMMA*gex{n} + i*ones(size(gex{n}))), ...
        z, f+fshift));
    mxy_shift{n} = -mz{n} .* mxytemp_shift;
   
    figure(n)
    subplot(331)
    imagesc(f, z, mz{n})
    colorbar
    subplot(332)
    plot(z, mz{n}(:, 101))
    subplot(333)
    plot(f, mz{n}(101, :))
    subplot(334)
    imagesc(f, z, imag(mxy{n}))
    colorbar
    subplot(335)
    cplot(z, mxy{n}(:, 101))
    subplot(336)
    cplot(f, mxy{n}(101,:))
    subplot(337)
    imagesc(f, z, imag(mxy_shift{n}))
    colorbar
    subplot(338)
    cplot(z, mxy_shift{n}(:, 101))
    subplot(339)
    cplot(f, mxy_shift{n}(101,:))
    
    %mxy_z{n} = sum(mxy{n}(Ivox,:), 1);
    mxy_z{n} = sum(mxy{n}, 1);
    mxy_z_shift{n} = sum(mxy_shift{n}, 1);
    
    figure(100+n)
    subplot(211)
    cplot(f, mxy_z{n})
    subplot(212)
    cplot(f, mxy_z_shift{n})
    
    t = [1:length(g{n})]*dt;
    figure(200+n)
    subplot(211)
    cplot(t,rf{n} / (2*pi*GAMMA*dt))
    subplot(212)
    plot(t,g{n})
    
%     disp(['Im = ' num2str(sum(imag(mxy{n}(:,101)))) ...
%         ', Re = ' num2str(sum(real(mxy{n}(:,101))))])
end

%rfwrite(rfscaleg(real(rf_pyr), pw_rf*1e3, GAMMA), pw_rf, pi/2, GAMMA, pw_rf/2, bw);