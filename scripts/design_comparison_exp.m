addpath ~/matlab/steam/

dT = 3e-3;
peakB1 = 1.5;  % C13 coil maxes
dt = 40e-6;
GAMMA = 1071;
B1_scale = .6:.1:1.4;

xsim = [-1:.01:1];
T1 = 30; T2 = .195;

do_plots = 0;

%% conventional 90-90
Nw(1) = ceil(pi/2 / (2*pi*GAMMA*dt*peakB1));
Nwait(1) = round(dT/dt - Nw(1));
rf{1} = zeros(1,2*Nw(1) + Nwait(1));
rf{1}(1:Nw(1)) = pi/2 / Nw(1);  rf{1}(end-Nw(1)+1:end) = pi/2 / Nw(1);


%%
N = 18;
params.dt = dt; params.GAMMA = GAMMA; params.peakB1 = peakB1;

%% SLR
params.tbw = N/2;  params.method = 'SLR';
%params.d1 = 1e-3;params.d2 = 1e-2;


[rf{2} Nw(2) Nwait(2)] = design_sSTE_inv(N, dT, params);
rf{2} = rf{2}*1.08;  % factor for B1 experimental data.

%% Sech
params.method = 'sech';
params.mu = 4; params.tbw = 4;

[rf{3} Nw(3) Nwait(3)] = design_sSTE_inv(N, dT, params);
rf{3} = rf{3} * 1;

Senc_xy = zeros(length(B1_scale), length(rf));
Senc_z = Senc_xy; Sfinal_90 =  Senc_xy; Sfinal_g90 =  Senc_xy; Sprep =  Senc_xy;
for S = B1_scale
    
    %% Menc simulation
    Menc_xy = zeros(length(xsim), length(rf)); Menc_z = zeros(length(xsim), length(rf));
    for n = 1:length(rf)
        [atemp btemp] = abr(S*rf{n}, xsim *length(rf{n}) *dt/dT);
        Menc_xy(:,n) = ab2ex(atemp,btemp);
        Menc_z(:,n) = ab2inv(atemp,btemp);
        %[mx my mz] = bloch(rfscaleg(rf{n}, dt*1e3*length(rf{n})), zeros(size(rf{n})), ...
        %    dt, T1, T2, xsim/dT, 0, 0);
    end
    
    if do_plots
        figure(1)
        plot(xsim,Menc_z)
        figure(2)
        cplot(xsim,Menc_xy)
    end
    
    %% excitations
    
    % 90
    rf_90 = [-pi/2 zeros(1, round(dT/dt-1))];
    mxy_90 = ab2ex(abr(S*rf_90, xsim * length(rf_90)*dt/dT));
    
    Mfinal_90 = zeros(length(xsim), length(rf));
    for n = 1:length(rf)
        Mfinal_90(:,n) = Menc_z(:,n) .* mxy_90;
    end
    
    if do_plots
        figure(3)
        cplot(xsim, Mfinal_90)
    end
    
    %% gapped 90
    params.tbw = (N-1)/2;
    [rf_g90 Nw_g90 Nwait_g90] = design_sSTE_ex((N-1), dT, params);
    
    % add spin-echo pulse for refocusing
    
    Nse_wait = round( length(rf_g90) /2 );
    rf_g90 = [rf_g90 pi zeros(1,Nse_wait)];
    
    mxy_g90 = ab2ex(abr(S*rf_g90, xsim * length(rf_g90)*dt/dT));
    
    Mfinal_g90 = zeros(length(xsim), length(rf));
    for n = 1:length(rf)
        Mfinal_g90(:,n) = Menc_z(:,n) .* mxy_g90;
    end
    
    if do_plots
        figure(4)
        cplot(xsim, Mfinal_g90)
        
        % figure(6), mpplot(mxy_g90)
    end
    
    %% prep
    
    Mprep = Menc_z .^2 * sin(pi/2 * S);
    
    if do_plots
        figure(5)
        plot(xsim, Mprep);
    end
    
    %%
    if do_plots
        for f = 1:5
            figure(f)
            axis([-1 1 -1 1])
        end
    end
    
    
    %%
    IS = find(S==B1_scale);
    Senc_xy(IS,:) = sum(Menc_xy,1)/length(xsim);
    Senc_z(IS,:) = sum(abs(Menc_z),1)/length(xsim);
    Sfinal_90(IS,:) = sum(Mfinal_90,1)/length(xsim);
    Sfinal_g90(IS,:) = sum(Mfinal_g90,1)/length(xsim);
    Sprep(IS,:) = sum(Mprep,1)/length(xsim);
    
end

%% b-value (inc T2?) calculations - < 1% expected difference
pw_ramp = 560e-6;  Nramp = round(pw_ramp/dt);
a_g = 0.38; D = [0 1e-3];  % water D = 2.3e-3;
Nx = 101;
    Ng = dT/dt - 2*Nramp - Nw;

for n =  []%1:length(rf)
    g{n} = zeros(size(rf{n}));
    for In = 1:floor(length(rf{n})*dt/dT)
        g{n}((1:Ng+2*Nramp) + (In-1)*round(dT/dt) + Nw(n)) = ...
            [(1:Nramp)/Nramp ones(1,Ng(n)), (Nramp:-1:1)/Nramp] * a_g;
    end
    
    
    for Ds = D
        
        Nbuff = 2*length(rf{n});
        x = [-Nbuff:Nx-1+Nbuff] * 2*pi / (Nx-1) - pi;
        M = repmat([0;0;1], [1 length(x)]);
        for puls = 1:2
            for It = 1:length(rf{n})
                M = xyrot(rf{n}(It)) * M;
                for xn = x
                    M(:, xn==x) = zrot(g{n}(It)*xn / (a_g * (Ng(n)+Nramp+1))) * M(:, xn==x);
                end
                dM = Ds * (GAMMA*a_g*(Ng(n)+Nramp+1)*dt * Nx)^2 * diff(M,2,2) * dt;
                if max(abs(dM) > 1e-2)
                    disp('warning: large dM');
                end
                M = M(:,2:end-1) + dM;  % check for too large dM??
                x = x(2:end-1);
                if 0
                    figure(99)
                    subplot(211), plot(M.'), subplot(212), plot(dM.'),drawnow
                end
            end
            M(1:2,:) = 0;
        end
        Senc(n,Ds==D) = sum(M(3,:))/Nx;
        
    end
end

% approximate:  looks ok.
% b_STE = (GAMMA * 2*pi *a_g * (Ng+Nramp+1)*dt)^2 * (2/3)*(Ng+Nramp+1)*dt
% .5*exp(-b_STE*D)


%% FINAL figures

B1_sim = B1_scale; load ../data/phantom_B1_data.mat

figure(6)
subplot(121)
plot(B1_sim, abs(Sfinal_90), '--', ...
    B1_scale, Sphantom(:, 1:3) ./ repmat(T2scale(1:3), [length(B1_scale), 1]) / Sphantom(5,1) / 2, 'x' )
axis([.6 1.4 .25 1])
legend(Iphantom(1:3))
xlabel('B_1 Scale'), ylabel('Echo Amplitude')

subplot(122)
plot(B1_sim, abs(Sprep), '--', ...
    B1_scale, Sphantom(:, 4:6) ./ repmat(T2scale(4:6), [length(B1_scale), 1]) / Sphantom(5,4), 'x' ,...
    B1_sim, sin(pi/2 * B1_sim), '--')
axis([.6 1.4 .25 1])
legend(Iphantom(4:6))
xlabel('B_1 Scale'), ylabel('Echo Amplitude')