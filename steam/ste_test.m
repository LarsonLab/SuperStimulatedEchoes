% To test STE pulse design:

GAMMA = 1070.5;

z = linspace(-40, 40, 201);
f = linspace(-400, 400, 201);

%fmod = 

TE = 2e-3;
dt = 10e-6;
NTE = round(TE/dt);

% conventional
n = 1;
gamp = 0.1;
rf{n} = [pi/2, zeros(1, NTE), pi/2];
g{n} = [gamp, ones(1,NTE) * gamp, gamp];
rfex{n} = [pi/2, zeros(1,NTE)];
gex{n} = [gamp, ones(1,NTE)*gamp];

if 0
    % % BIR-4?
    % nope
    n = n+1;
Tinv = 4e-3;
ninv = round(Tinv/dt);
BW = 1000; mu = 6;
%NTE = 0;
rf_ad = adiabatic_sech(ninv, 1*2*pi*GAMMA*Tinv/ninv, BW*Tinv*2*pi/1.87, mu);
rf{n} = [rf_ad([ninv/2+1:ninv]), zeros(1,NTE/2), rf_ad*exp(j*3*pi/2), zeros(1,NTE/2), rf_ad(1:ninv/2)];
g{n} = [ones(1,length(rf{n})) * 0.02];
rfex{n} = [rf_ad([1:ninv/2]), zeros(1,NTE)];
gex{n} = [ones(1,length(rfex{n}))*0.02];
end

if 0
% min/max phase
% same as hard pulse
n = n+1;
gamp = 0.1;
rfm = dzrf(50, 6, 'ex', 'min');
rf{n} = [rfm, zeros(1, NTE), rfm(end:-1:1)];
g{n} = [zeros(1,50), ones(1,NTE) * gamp, zeros(1,50)];
rfex{n} = [rfm, zeros(1,NTE)];
gex{n} = [zeros(1,50), ones(1,NTE)*gamp];
end


n = n+1;
gamp = 0.4;
gamp2 = 0.4;
N = 24; % even
tbw = N/2;
% real filter:
tw = 0.05;
%finv = [0 .1 .2 .4 .6 .8 .9 1];a = [1 1 0 0 1 1 0 0];
finv = [0 .5-tw .5+tw 1];a = [1 1 0 0];
binv = firls(N, finv, a, [1 .01]);
rfinv = b2rf(binv);
Nwait = 16;
rf{n} = zeros(1, Nwait*(N+1));
rf{n}(1:Nwait:end) = rfinv;
g{n} = ones(size(rf{n}))*gamp2;
g{n}(1:Nwait:end) = gamp;
%Ngrad = 19*(N+1);
Ngrad = Nwait-1; %round(NTE*0.33);
rfex{n} = [pi/2, zeros(1,Ngrad)];
gex{n} = [gamp, gamp2*ones(1,Ngrad)];

for n = 1:length(rf)
    mz{n} = ab2inv(abr(rf{n}, 2*pi*dt*(GAMMA*g{n} + i*ones(size(g{n}))), ...
        z, f));
    mxytemp = ab2ex(abr(rfex{n}, 2*pi*dt*(GAMMA*gex{n} + i*ones(size(gex{n}))), ...
        z, f));
    mxy{n} = mz{n} .* mxytemp;
    figure(n)
    subplot(221)
    imagesc(f, z, mz{n})
    colorbar
    subplot(222)
    plot(z, mz{n}(:, 101))
    subplot(223)
    imagesc(f, z, imag(mxy{n}))
    colorbar
    subplot(224)
    cplot(z, mxy{n}(:, 101))
    disp(['Im = ' num2str(sum(imag(mxy{n}(:,101)))) ...
        ', Re = ' num2str(sum(real(mxy{n}(:,101))))])
end