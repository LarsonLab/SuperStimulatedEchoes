%% simulation data
m0 = [1; 0.0]; % P, PH

k1 = .0019; % 1/s
k2 = .002 *20;
A = [k1 -k2; -k1 k2];

t = 0:20; % s

% syringe data /home/plarson/data/steam/20120217_HPpyrh
% exp 2
% estimated kP->PH = 0.0019 1/s
TE = 13e-3;



dph = mod(TE * 270 / 2, 1);

for It = 1:length(t)
    m(:,It) = expm(-A*t(It)) * m0;
    mc(:, It) = expm(-A*t(It)) .* [1 exp(i*dph * 2*pi); exp(-i*dph * 2*pi) 1] * m0;

end
% 
% subplot(311)
% plot(t, m)
% subplot(312)
% cplot(t, mc(1,:))
% subplot(313)
% cplot(t, mc(2,:))
% 





A0 = abs(Amets(2,:));
close all
subplot(311)
plot(t, m,'--', t, abs(Amets(2,:))./A0, t, abs(Amets(1,:))./A0)
subplot(312)
cplot(t, mc(1,:),'--')
hold on, cplot(t, Amets(2,:)./A0)
subplot(313)
cplot(t, mc(2,:),'--')
hold on, cplot(t, Amets(1,:)./A0)
