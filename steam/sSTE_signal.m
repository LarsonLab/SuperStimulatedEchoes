function Ss = sSTE_signal(ds, a, bn, N, Sds, Sa, T1, T2, GAMMA)

if ds <= 0
    Ss = 0;
    return;
end

fsim = linspace(-1/(2*(ds+a)), 1/(2*(ds+a)), 201);
Mz_0 = 1;

params.dt = (ds+a) / 100;
if params.dt > T2/10
    params.dt = (ds+a) / ceil(10 * (ds+a) / T2);
end

params.method = 'SLR';
params.tbw = max(floor(N/2), 3);
params.GAMMA = GAMMA;
params.peakB1 = 1e3; % for short pulses
rf = design_sSTE_inv(N, ds+a, params);
rfsim = rfscaleg(rf, params.dt*length(rf)*1e3, GAMMA);


[mx, my, mz] = bloch(rfsim, zeros(size(rfsim)), params.dt, T1, T2, ...
    fsim, 0, 0, 0, 0, Mz_0, GAMMA);
Mz_store = sum(abs(mz))/length(fsim);

% bn = normalized (b/G^2)
%Tm = bn ./ (2*pi*GAMMA*ds).^2  - ((N-1)*(ds+a) - ds/3);
Tm = bn ./ (2*pi*GAMMA*ds).^2  - Sds*ds - Sa * a;
if Tm < 0
    Tm = 0;
end

Ss = Mz_store^2 * exp(-Tm/T1);
