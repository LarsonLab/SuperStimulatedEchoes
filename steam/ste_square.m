% To test STE pulse design:

% can this do 3D box selection?

% using some parameters from successful C13 experiments

GAMMA = 1070.5;

z = linspace(-10, 10, 201);  % mm

Ivox = find(abs(z) < 2.5);

%TE = 20e-3;
% designing based on pulse spacing
dT = 3e-3;
f = linspace(-3/dT, 3/dT, 201);  % Hz
a_gencode = 3 / 10;  % G/mm
pw_gencode = 2e-3;

dt = 40e-6;%16*4e-6;

% conventional
n = 1;
TE(n) = dT*2;
pw_rf = 264e-6;
peak_B1 = pi/2 / (2*pi*pw_rf * GAMMA);
peak_B1 = 1.5;  % G - C13 coil maxes

Nrf = round(pw_rf/dt);
Ngencode = round(pw_gencode/dt);
Nwait = round((TE(n)/2 - pw_rf - pw_gencode)/dt);

rf{n} = [pi/2*ones(1,Nrf)/Nrf, zeros(1, Ngencode+Nwait), pi/2*ones(1,Nrf)/Nrf];
g{n} = [zeros(1,Nrf), ones(1,Ngencode)*a_gencode, zeros(1,Nwait+Nrf)];

rfex{n} = [pi/2*ones(1,Nrf)/Nrf, zeros(1, Ngencode+Nwait)];
gex{n} = [zeros(1,Nrf), ones(1,Ngencode)*a_gencode, zeros(1,Nwait)];

% some linear phase accrual in frequency after final RF pulse...
% not dependant on TE, but on RF duration


% inversion encoding
n = n+1;
pulseselect = 11;

N = 6;
tbw = 3;
switch pulseselect
    case 1
        N = 12;
        d1 = .01/8; d2 = sqrt(.01/2);
        d = dinf(d1,d2);
        % real filter:
        fw = d/tbw;
        finv = [0 .5*(1-fw) .5*(1+fw) 1];a = [1 1 0 0];
        binv = firls(N-1, finv, a, [1 d1/d2]);
        rfinv = b2rf(binv);
        
        rfinv = real(rfinv);
    case 2
        % adiabatic inversion
        N = 12;
        mu = 3.;
        rfinv = adiabatic_sech(N, 1, 0.9* tbw * 2*pi/1.87, mu);
        rfinv = 1.7*pi*rfinv/sum(rfinv);
    case 3
        % adiabatic inversion
        N = 8;
        mu = 2.;
        rfinv = adiabatic_sechn(N, 1,  0.7*tbw * 2*pi/1.87, mu, 2);
        rfinv = 1.2*pi*rfinv/sum(rfinv);
    case 5
        % for SLR matched design
        N = 9;  % seems to breakdown (inversion design mostly) for N=7
        tbw = 3;
        d1 = .01/8; d2 = sqrt(.01/2);
        d = dinf(d1,d2);
        % real filter:
        fw = d/tbw;
        finv = [0 .5*(1-fw) .5*(1+fw) 1];a = [1 1 0 0];
        binv = firls(N-1, finv, a, [1 d1/d2]);
        rfinv = b2rf(binv);
        
        rfinv = real(rfinv);
    case 6
        % adiabatic inversion for matched design
        N = 11;
        mu = 2.5; tbw = 4;
        rfinv = adiabatic_sech(N, 1, 0.65* tbw * 2*pi/1.87, mu);
        rfinv = 1.5*pi*rfinv/sum(rfinv);
    case 7
        % adiabatic inversion for matched design
        N = 9;
        mu = 2.5; tbw = 3;
        rfinv = adiabatic_sech(N, 1, 0.65* tbw * 2*pi/1.87, mu);
        rfinv = 1.5*pi*rfinv/sum(rfinv);
    case 10
        % longer SLR (N~20, based on simulations in
        % square_tradeoff_dnp.m)
        N = 18;
        tbw = N/2;
        d1 = .01/8; d2 = sqrt(.01/2);
        
        if 0 % max phase
            n2 = 2*N-1;
            d = 0.5*dinf(2*d1,0.5*d2.*d2);
            w = d/tbw;
            finv = [0 .5*(1-fw) .5*(1+fw) 1];a = [1 1 0 0];
            binv2 = firpm(n2-1,finv,a,[1 2*d1/(0.5*d2*d2)]);
            binv = fmp(binv2);
        else % lin phase
            d = dinf(d1,d2);
            fw = d/tbw;
            finv = [0 .5*(1-fw) .5*(1+fw) 1];a = [1 1 0 0];
            binv = firls(N-1, finv, a, [1 d1/d2]);
        end
        rfinv = b2rf(binv);
        
        rfinv = real(rfinv);
    case 11
        % longer adiabatic inversion
        % (N~20, based on simulations in square_tradeoff_dnp.m)
        N = 18;
        mu = 4; tbw = 4;
        rfinv = adiabatic_sech(N, 1, tbw * 2*pi/1.87, mu);
        rfinv = 1.5*pi*rfinv/sum(rfinv);

end

TE(n) = dT * N;
temp_B1 = abs(rfinv) / (2*pi*GAMMA*dt);
Nw_all = ceil(temp_B1/peak_B1);
Nw = max(Nw_all);

Nwait = round((TE(n)/ dt - (N*Nw+1)) / N);

rf{n} = zeros(1, Nwait*(N-1) + N*Nw);

for w = 1:N
        rfadd = 0;%mod(w,2)*i*pi;
        rfsign = 1;%exp(i*pi*w);
    Nw_all(w) = Nw;
    rf{n}([1:Nw_all(w)]+round((Nw-Nw_all(w))/2) + (w-1)*(Nwait+Nw))= (rfinv(w)*rfsign + rfadd)/Nw_all(w);
end


if 0
    nompw = length(rf{n})*dt;
    rfwrite(rfscaleg(rf{n}, nompw*1e3,GAMMA), nompw, pi, GAMMA, nompw/2, 0);
    disp('delta_rf1_lobes = ')
    disp((Nw+Nwait)*dt)
    disp('pw_rf1_lobes = ')
    disp(Nw*dt)
    return
end


% gradient design
a_gsquare = 2.77 / 10;  % G/mm
g{n} = zeros(size(rf{n}))*a_gsquare;
if 1
    max_slew = 1.5e3; % G/mm/s
    Nramp = ceil(a_gsquare / max_slew / dt);
    Nplat = Nwait - 2*Nramp;
    globe = [ [1:Nramp]/Nramp ones(1,Nplat) [Nramp:-1:1]/Nramp];
        
    for w = 1:N-1
        gcsign = 1;%mod(w,2)*2 -1;
        g{n}([1:Nwait]+Nw+ (Nwait+Nw)*(w-1)) = globe*a_gsquare*gcsign;
    end
    
    % with gradient on during pulses, get nasty rect spatial profile
    % could fix by using spectral-spatial pulse design methods:
    % - modulate inversion samples by sincs or better
    % - adjust/ramp gradient as needed for slice thickness
    Ngrad = Nwait;
    gex{n} = [zeros(1,Nrf), a_gsquare*globe];
elseif 1
    g{n} = ones(size(rf{n}))*a_gsquare;
    Ngrad = Nwait+1;
    gex{n} = [zeros(1,Nrf), a_gsquare*ones(1,Ngrad)];
else
    % matched to 90-90
    a_gsquare = Ngencode*a_gencode / Nwait;
    
    for w = 1:N-1
        g{n}([1:Nwait]+Nw+ (Nwait+Nw)*(w-1)) = a_gsquare;
    end

    gex{n} = [zeros(1,Nrf), a_gsquare*ones(1,Nwait)];
  
end

exselect = 4;

switch exselect
    case 0
    rfex{n} = [pi/2*ones(1,Nrf)/Nrf, zeros(1,Ngrad)];
%gex{n} = [a_gsquare*ones(1,Nrf), a_gsquare*ones(1,Ngrad)];
    case 1
        foff = 1/(2*dt*(Nwait + Nw));
    % This strategy is no good
    Nex = N;
  rf90 = dzrf(Nex, 1.1*tbw, 'ex');
temp_B1 = max(abs(rf90)) / (2*pi*GAMMA*dt);
Nwex = ceil(temp_B1/peak_B1);

Nwaitex = (Nwait+Nw) - Nwex;
Nwaitend = 0;%round((Nwait+Nw) - Nwex/2);

rfex{n} = zeros(1, Nwaitex*(Nex-1) +Nwaitend + Nex*Nwex);
for w = 1:Nwex
    rfex{n}([w:(Nwaitex+Nwex):end]) = rf90/Nwex;
end
gex{n} = zeros(size(rfex{n}));
    for w = 1:Nex-1
        gex{n}([1:Nwait]+Nwex+ (Nwaitex+Nwex)*(w-1)) = globe*a_gsquare;
    end
    
    rfex{n} = (1 + exp(j*2*pi*foff*[1:length(rfex{n})]*dt)).*rfex{n};
    
    case 2
    % phase-control in wideband pulse...??  too fine.
    Nex = 1*(Nw+Nwait);
    
    ffilt = [-1e3:10:1e3];
    fenc = 1 / (2*(Nw+Nwait)*dt);
    a = ones(1,length(ffilt));
    ph = (mod(ffilt+fenc/2,fenc)-fenc/2) * (pi/fenc);
    d = .01*a;
    
    [b, status] = fir_linprog_phase(Nex, ffilt*2*dt, a, ph, d, [], 1);
if status == 'Failed'
    error('Filter design failed')
end
    rfex{n} = b2rf(sin(pi/4)*b);
    gex{n} = zeros(size(rfex{n}));

    case 3
        foff = 1/(2*dt*(Nwait + Nw));
        rfex{n} = [-rf{n}.*(1 + exp(j*2*pi*foff*[1:length(rf{n})]*dt))/(1.7*2)];
        gex{n} = [g{n}];
    case 4
        % STE as prep pulses
        rfex{n} = [rf{n} -pi/2*ones(1,Nrf)/Nrf];
        gex{n} = [g{n} zeros(1,Nrf)];
    case 5
        d1 = sqrt(.01/2); d2 = sqrt(.01/2);
        d = dinf(d1,d2);
        % real filter:
        fw = d/tbw;
        fex = [0 .5*(1-fw) .5*(1+fw) 1];a = [1 1 -1 -1];
        bex = firls(N-1, fex, a, [1 d1/d2]);
        rfex_temp = b2rf(bex*sqrt(1/2));
rfex{n} = zeros(1, Nwait*(N-1) + N*Nw);
for w = 1:Nw
    rfex{n}([w:(Nwait+Nw):end]) = rfex_temp/Nw;
end
        
      % bex = (-binv + 0.5) * 2/sqrt(2);        rfex{n} = b2rf(bex);
        gex{n} = [g{n}];
        
       % end up with signal to be refocused by spin-echo?  (refocused at
       % center of pulse?)
        Nse_wait = round(( (N-1)*Nwait + (N-2)*Nw )/2 +8);
        rfex{n} = [rfex{n} pi zeros(1,Nse_wait)];
        gex{n} = [gex{n} 0 (N-1)/2 *globe*a_gsquare zeros(1,Nse_wait-length(globe))];
        
end

for n = 1:length(rf)
    mz{n} = ab2inv(abr(rf{n}, 2*pi*dt*(GAMMA*g{n} + i*ones(size(g{n}))), ...
        z, f));
    mxytemp = ab2ex(abr(rfex{n}, 2*pi*dt*(GAMMA*gex{n} + i*ones(size(gex{n}))), ...
        z, f));
    mxy{n} = -mz{n} .* mxytemp;
    figure(n)
    subplot(331)
    imagesc(f, z, mz{n})
    colorbar
    subplot(332)
    plot(z, mz{n}(:, 101))
    subplot(333)
    plot(f, mz{n}(101, :))
    subplot(334)
    imagesc(f, z, imag(mxytemp))
    colorbar
    subplot(335)
    cplot(z, mxytemp(:, 101))
    subplot(336)
    cplot(f, mxytemp(101,:))
    subplot(337)
    imagesc(f, z, imag(mxy{n}))
    colorbar
    subplot(338)
    cplot(z, mxy{n}(:, 101))
    subplot(339)
    cplot(f, mxy{n}(101,:))
    
    %mxy_z{n} = sum(mxy{n}(Ivox,:), 1);
    mxy_z{n} = sum(mxy{n}, 1);
    
    figure(100+n)
    plot(f, real(mxy_z{n}),f, imag(mxy_z{n}),f, abs(mxy_z{n}),'--')
    
    t = [1:length(g{n})]*dt;
    figure(200+n)
    subplot(211)
    cplot(t,rf{n} / (2*pi*GAMMA*dt))
    subplot(212)
    plot(t,g{n})
    
    disp(['Im = ' num2str(sum(imag(mxy{n}(:,101)))) ...
        ', Re = ' num2str(sum(real(mxy{n}(:,101))))])
end

return
%% EPG
train = cell(2*N-1,1); trainargs =train;

for w=1:N
train{2*w-1} = 'T';
trainargs{2*w-1} = rfinv(w);

if w < N
train{2*w} = 'S';
trainargs{2*w} = round(exp(i*pi* w * 2));
end
end
figure
[F,k] = epg(train, trainargs, [], [], 1);

train = cell(2*N-1,1); trainargs =train;

% disp(sum(abs(F(2:3:end))) ) %?

%%
sc = 0.8:0.2:1.2;
for n = 1:2
    for Isc = 1:length(sc)
        mz{n, Isc} = ab2inv(abr(sc(Isc)*rf{n}, 2*pi*dt*(GAMMA*g{n} + i*ones(size(g{n}))), ...
            z, f));
        mxytemp = ab2ex(abr(sc(Isc)*rfex{n}, 2*pi*dt*(GAMMA*gex{n} + i*ones(size(gex{n}))), ...
            z, f));
        mxy{n,Isc} = -mz{n,Isc} .* mxytemp;
    end
end

figure(99)

subplot(211)
    plot(z, mz{1,1}(:, 101), 'b--', z, mz{1,2}(:, 101), 'b-', z, mz{1,3}(:, 101), 'b--', ...
        z, mz{2,1}(:, 101), 'g--', z, mz{2,2}(:, 101), 'g-', z, mz{2,3}(:, 101), 'g--')
ylabel('M_{Z}')
%     subplot(222)
%     plot(f, mz{1,1}(101, :), 'b--', f, mz{1,2}(101, :), 'b-', f, mz{1,3}(101, :), 'b--', ...
%         f, mz{2,1}(101, :), 'g--', f, mz{2,2}(101, :), 'g-', f, mz{2,3}(101, :), 'g--')
    subplot(212)
    plot(z, imag(mxy{1,1}(:, 101)), 'b--', z, imag(mxy{1,2}(:, 101)), 'b-', z, imag(mxy{1,3}(:, 101)), 'b--', ...
        z, imag(mxy{2,1}(:, 101)), 'g--', z, imag(mxy{2,2}(:, 101)), 'g-', z, imag(mxy{2,3}(:, 101)), 'g--')
xlabel('Position (mm)'), ylabel('M_{X}'), ylim([0 1])
%     subplot(224)
%     plot(f, imag(mxy{1,1}(101, :)), 'b--', f, imag(mxy{1,2}(101, :)), 'b-', f, imag(mxy{1,3}(101, :)), 'b--', ...
%         f, imag(mxy{2,1}(101, :)), 'g--', f, imag(mxy{2,2}(101, :)), 'g-', f, imag(mxy{2,3}(101, :)), 'g--')
