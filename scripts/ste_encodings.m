Nperiods = 2;
Nsamp = 100;
t = [-pi*Nperiods:pi/Nsamp:pi*Nperiods];
tn = t/(2*pi);
Q = sign(cos(t));
Q(find(abs(cos(t))< 1e-3)) = 0;

% STEAM

figure(1)
plot(tn, sin(t),'k--', tn, cos(t),'k-')
xlabel('\phi / 2\pi'), ylabel('M_{enc}'), axis([-Nperiods/2 Nperiods/2, -1 1])
legend('M_Y', 'M_Z')

figure(2)
plot(tn, cos(t).^2,'k-', tn, cos(t).*sin(t),'k--')
xlabel('\phi / 2\pi'), ylabel('M_{final}'), axis([-Nperiods/2 Nperiods/2, -1 1])
legend('M_X', 'M_Y')

% sSTE

figure(3)
plot(tn,sqrt(1-Q.^2),'k--', tn, Q,'k-')
xlabel('\phi / 2\pi'), ylabel('M_{enc}'), axis([-Nperiods/2 Nperiods/2, -1 1])
legend('M_Y', 'M_Z')

figure(4)
plot(tn, cos(t).*Q,'k-', tn, Q.*sin(t),'k--')
xlabel('\phi / 2\pi'), ylabel('M_{final}'), axis([-Nperiods/2 Nperiods/2, -1 1])
legend('M_X', 'M_Y')


figure(5)
plot(tn, Q.^2,'k-')
xlabel('\phi / 2\pi'), ylabel('M_{final}'), axis([-Nperiods/2 Nperiods/2, -1 1])
legend('M_X')
