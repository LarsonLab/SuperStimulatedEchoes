clear all

GAMMA = 1070.5;

%% phase spec:
mets = setup_C13_mets;
f_met = ([mets.pyr.Hz mets.ala.Hz mets.lac.Hz]' - (mets.pyr.Hz));

dt = 5e-5;
T2 = Inf; T1 = Inf;
dT = 5.5 / (mets.lac.Hz - mets.pyr.Hz) /2; % time between RF pulses
f = linspace(-5/dT, 5/dT, 501);

G = .4;% G/mm

%% 90-90-90 tests
% linear - gets refocused, max phase - could work

Nrf1 = 10; rf1 = exp(i*pi*.2)* ones(1,Nrf1) * pi/2 / Nrf1;

Nrf2 = Nrf1; rf2 = exp(i*pi*.1)* ones(1,Nrf2) * pi/2 / Nrf2;
Nrf3 = 10; rf3 = exp(i*pi*.4)* ones(1,Nrf3) * pi/2 / Nrf3;

%% gradient
Ng = round(dT/dt - max([Nrf1,Nrf2,Nrf3]) -1);
z = linspace(-1/(Ng*dt*G*GAMMA), 1/(Ng*dt*G*GAMMA), 201);  % mm

rfenc = [rf1, zeros(1, dT/dt - Nrf1/2 - Nrf2/2), rf2];
genc = [zeros(1,Nrf1), ones(1, Ng), zeros(1,length(rfenc)-Ng-Nrf1)] * G;

rfex = [rf3, zeros(1, dT/dt - (Nrf3)/2 +1)];
gex = [zeros(1,Nrf3), ones(1, Ng), zeros(1,length(rfex)-Ng-Nrf3)]* G;


%%
mz = ab2inv(abr(rfenc, 2*pi*dt*(GAMMA*genc + 1i*ones(size(genc))), z, f_met));
mxy = ab2ex(abr(rfex, 2*pi*dt*(GAMMA*gex + 1i*ones(size(gex))), z, f_met));

figure(98), plot(z,mz)
figure(99), cplot(z,mxy)

[a1 b1] = abr([rf1 pi zeros(1,Nrf1/2)], 2*pi*dt*ones(1,length(rf1) + Nrf1/2+1), f_met);
[a2 b2] = abr([rf2 pi zeros(1,Nrf2/2)], 2*pi*dt*ones(1,length(rf2) + Nrf2/2+1), f_met);
[a3 b3] = abr([rf3 pi zeros(1,Nrf3/2)], 2*pi*dt*ones(1,length(rf3) + Nrf3/2+1), f_met);

for m1 = 1:length(f_met)
    for m2 = 1:length(f_met)
        S(m1, m2) = sum(mz(:,m1) .* mxy(:,m2))/length(z);

        % very close to correct phase from pulses...
        Spred(m1,m2) = 2*conj(a1(m1))*b1(m1)*...
            -conj(a2(m1)) * conj(b2(m1))  *...
            2*conj(a3(m2))*b3(m2) *...
            exp(1i * 2*pi*(f_met(m2) - f_met(m1)) * dT);
    end
figure(1)
    subplot(length(f_met),1,m1)
    cplot(f_met, S(m1,:))
figure(2)
    subplot(length(f_met),1,m1)
    cplot(f_met, -conj(Spred(m1,:)))
end

angle(S)/pi
