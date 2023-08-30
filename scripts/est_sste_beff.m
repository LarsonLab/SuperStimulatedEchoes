N = [18 2];
do_plots = 0;

GAMMA = 1071; % Hz/G
G = .4; % G/mm
TG = 2e-3;
a = 0;
TM = 0;
dk = GAMMA* 2*pi*G*TG;
Snorm = 1/(dk^2 * TG);

%% Trf test
TM =0;
a = [0:.1:1]*TG; % linear!
for atest = a

[S(:,atest==a) Dtest] = calc_sste_diff(N, dk, [TG atest], TM, [], [], do_plots);
pfit = polyfit(Dtest(1:4), log(S(1:4,atest==a)).', 1);
beff(atest==a) = -pfit(1);
Sfit = exp(-beff(atest==a)* Dtest);

figure(99), semilogy(Dtest, S(:,atest==a)/S(1,atest==a), Dtest, Sfit, '--') , ylim([1e-3 1])
drawnow
end

pfit_a = polyfit(a/TG, beff*Snorm, 1)
% SG Sa

plot(a/TG, beff*Snorm)

%btest = dk^2 * (19.4*TG + 17.3 * a)


%% TM test
a =0;
TM = [0:10]*TG;
for TMtest = TM

[S(:,TMtest==TM) Dtest] = calc_sste_diff(N, dk, [TG a], TMtest, [], [], do_plots);
pfit = polyfit(Dtest(1:4), log(S(1:4,TMtest==TM)).', 1);
beff(TMtest==TM) = -pfit(1);
Sfit = exp(-beff(TMtest==TM)* Dtest);

figure(99), semilogy(Dtest, S(:,TMtest==TM)/S(1,TMtest==TM), Dtest, Sfit, '--') , ylim([1e-3 1])
drawnow
end

pfit_a = polyfit(TM/TG, beff*Snorm, 1)

plot(TM/TG, beff*Snorm)
%btest = dk^2 * (19.4*TG + TM)

return
%% test values
