%%
clear all

GAMMA = 1070.5;

dt = 5e-5;
T2 = Inf; T1 = Inf;
dT = 7e-3; % time between RF pulses
f = linspace(-5/dT, 5/dT, 501);

G = .4;% G/mm

%% square encode??
% [rf Nw Nwait rfinv] = design_sSTE_inv(N, dT, params)


%% phase spec:  add hydrate?
mets = setup_C13_mets;
f_met = ([mets.pyr.Hz mets.ala.Hz mets.lac.Hz]' - (mets.pyr.Hz+mets.lac.Hz)/2);
ph_comp = pi/2 - mod((f_met-f_met(1))*2*pi*dT, pi);
ph_comp(1) = 0;
df = 10;

%% 90-90-90 tests
% linear - gets refocused  max phase - could work
Trf1 = 3e-3 + dt; tbw =4;
Nrf1 = Trf1/dt;
rf1 = dzrf(Nrf1, tbw, 'ex', 'max');
Trf1 = .2e-3 + dt; Nrf1 = Trf1/dt; rf1 = ones(1,Nrf1) * pi/2 / Nrf1;
%phase pulse
f_pulse([1:length(f_met)]*2-1) = (f_met-df) * dt;
f_pulse([1:length(f_met)]*2) = (f_met+df) * dt;
ph_pulse([1:length(f_met)]*2-1) = ph_comp;
ph_pulse([1:length(f_met)]*2) = ph_comp;
a_pulse = ones(size(f_pulse));
if 1
  f_pulse = [-1,-.9,f_pulse,0.9,1];
  a_pulse = [0,0,a_pulse,0,0];
  ph_pulse = [0,0,ph_pulse,0,0];
end
d_pulse = ones(size(f_pulse)) * 1e-2;

de = rf_ripple(d_pulse, a_pulse, pi/2, 'ex');

% b = cfirpm(Nrf1-1, f_pulse, a_pulse.*exp(1i*ph_pulse),'both');
% freqz(b,1,256,'whole');
[b,status] = fir_linprog_phase(Nrf1,f_pulse, a_pulse, ph_pulse, de, [], 2);
%rf1 = b2rf(sqrt(1/2)*b);


Trf2 = .2e-3 + dt;
Nrf2 = Trf2/dt;
%rf2 = dzrf(Nrf2, tbw, 'ex', 'max');
rf2 = ones(1,Nrf2) * pi/2 / Nrf2;

%Trf3 = 1.1e-3 + dt; Nrf3 = Trf3/dt; rf3 = ones(1,Nrf3) * pi/2 / Nrf3;
Trf3 = Trf1; Nrf3 = Nrf1; rf3 = rf1; %dzrf(Nrf3, tbw, 'ex', 'max');

%%
Ng = dT/dt - max([Nrf1,Nrf2,Nrf3]);
z = linspace(-1/(Ng*dt*G*GAMMA), 1/(Ng*dt*G*GAMMA), 201);  % mm

rfenc = [rf1, zeros(1, dT/dt - Nrf1/2 - Nrf2/2), rf2];
genc = [zeros(1,Nrf1), ones(1, Ng), zeros(1,length(rfenc)-Ng-Nrf1)] * G;

rfex = [-rf3, zeros(1, dT/dt - Nrf3/2)];    
gex = [zeros(1,Nrf3), ones(1, Ng), zeros(1,length(rfex)-Ng-Nrf3)]* G;


%%
mz = ab2inv(abr(rfenc, 2*pi*dt*(GAMMA*genc + 1i*ones(size(genc))), z, f));
mxy = ab2ex(abr(rfex, 2*pi*dt*(GAMMA*gex + 1i*ones(size(gex))), z, f));

figure(1)
imagesc(f,z,mz, [-1 1])
colorbar
figure(2)
subplot(121), imagesc(f,z,real(mxy), [-1 1])
subplot(122), imagesc(f,z,imag(mxy), [-1 1])
colorbar

figure(3)
subplot(121), imagesc(f,z,real(mxy.*mz), [-1 1])
subplot(122), imagesc(f,z,imag(mxy.*mz), [-1 1])
colorbar

figure(4)
subplot(121),cplot(f, sum(mxy, 1))
subplot(122), cplot(f, sum(mxy.*mz, 1))

figure(5)
plot(f, mz(z==0,:))
figure(6)
mpplot(f, sum(mxy.*mz, 1))