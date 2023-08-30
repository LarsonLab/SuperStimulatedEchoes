% slice-selective sSTE

dT = 3e-3;
peak_B1 = 1.5;  % C13 coil maxes
dt = 40e-6;
GAMMA = 1071;

z = linspace(-20, 20, 201);  % mm
f = linspace(-2/dT, 2/dT, 201);  % Hz

do_plots = 0;


%%
N = 18;
params.dt = dt; params.GAMMA = GAMMA; params.peak_B1 = peak_B1;

%% Sech
params.method = 'sech';
params.mu = 4; params.tbw = 4;

rfinv = adiabatic_sech(N, 1, params.tbw * 2*pi/1.87, params.mu);
rfinv = 1.5*pi*rfinv/sum(rfinv);


%% sub pulses
tbw_x = 3;
Nw = 20;
rf_x = dzrf(Nw, tbw_x);
rf_x = rf_x/sum(rf_x);

TE = dT * N;
temp_B1 = max(abs(rfinv)) / (2*pi*params.GAMMA*params.dt);
%Nw = ceil(temp_B1/params.peak_B1);

Nwait = round((TE/ params.dt - (N*Nw+1)) / N);


gx_amp = 4/10; % G/mm
gdiff_amp = 4/10;

rf = zeros(1, Nwait*(N-1) + N*Nw);
g = ones(1, Nwait*(N-1) + N*Nw)*gdiff_amp;
for w = 1:N
    rf([1:Nw] + (w-1)*(Nwait+Nw)) = rf_x * rfinv(w);
    g([1:Nw] + (w-1)*(Nwait+Nw)) = gx_amp;
end

mz = ab2inv(abr(rf, 2*pi*dt*(GAMMA*g + 1i*ones(size(rf))), ...
    z, f));

figure(2)
imagesc(f,z,mz)

gx = ones(1, Nwait*(N-1) + N*Nw)*gdiff_amp;
gy = ones(1, Nwait*(N-1) + N*Nw)*gdiff_amp;
for w = 1:N/2
    gx([1:Nw] + (w-1)*(Nwait+Nw)) = gx_amp;
    gy([1:Nw] + (w-1)*(Nwait+Nw)) = 0;
end
for w = N/2+1:N
    gx([1:Nw] + (w-1)*(Nwait+Nw)) = 0;
    gy([1:Nw] + (w-1)*(Nwait+Nw)) = gx_amp;
end


mz = ab2inv(abr(rf, 2*pi*dt*GAMMA*(gx + 1i*gy), ...
    z, z));

figure(3)
imagesc(z,z,mz)

