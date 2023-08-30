clear all
do_plots = 0;

GAMMA = 1071; % Hz/G
G = .4; % G/mm
TG = 10e-3;

%% EPG

dk = GAMMA* 2*pi*G*TG;
Snorm = 1/(dk^2 * TG);  % Normalization problems ?

%% super STE
Ntest = [2 6 10 18 30]; %?
TMtest = [0 10 20] * TG;  % problem
%Nfit = 4:30;
%% SLR pulses
for TM = TMtest
    for N = Ntest
        [S_SLR(:,N==Ntest,TM==TMtest) Dtest] = ...
            calc_sste_diff(N, dk, TG, TM, [], [], do_plots);
        pfit(1:2, N==Ntest,TM==TMtest) = polyfit(Dtest(1:5), log(S_SLR(1:5,N==Ntest, TM==TMtest)).', 1);
        beff_SLR(N==Ntest,TM==TMtest) = -pfit(1, N==Ntest, TM==TMtest);
        S_SLR_fit(:,N==Ntest, TM==TMtest) = exp(-beff_SLR(N==Ntest,TM==TMtest)* Dtest);
    end
end

if 0
       figure(100), semilogy(Dtest/Snorm, S_SLR(:,:,1)) , ylim([1e-3 1])
figure(200), semilogy(Dtest/Snorm, S_SLR(:,:,2)), ylim([1e-3 1])
figure(300), semilogy(Dtest/Snorm, S_SLR(:,:,3)), ylim([1e-3 1])
 
else
    figure(100), semilogy(Dtest/Snorm, S_SLR(:,:,1) ./ repmat(S_SLR(1,:,1), [length(Dtest), 1]), ...
        Dtest/Snorm, S_SLR_fit(:,:,1), '--'), ylim([1e-3 1])
figure(200), semilogy(Dtest/Snorm, S_SLR(:,:,2) ./ repmat(S_SLR(1,:,2), [length(Dtest), 1]), ...
        Dtest/Snorm, S_SLR_fit(:,:,2), '--'), ylim([1e-3 1])
figure(300), semilogy(Dtest/Snorm, S_SLR(:,:,3) ./ repmat(S_SLR(1,:,3), [length(Dtest), 1]), ...
        Dtest/Snorm, S_SLR_fit(:,:,3), '--'), ylim([1e-3 1])
end

%% sech pulses
params_pulse.method = 'sech';
mu = [3 1.5 2.5 4 8]; tbw = [3 1.1 2 4 7]; adiab_scale = [1 1 1.2 1.5 1.5];


for TM = TMtest
    for N = Ntest
        params_pulse.mu = mu(N==Ntest);     params_pulse.tbw = tbw(N==Ntest);
        params_pulse.adiab_scale = adiab_scale(N==Ntest);  
        [S_sech(:,N==Ntest,TM==TMtest) Dtest] = ...
            calc_sste_diff(N, dk, TG, TM, [], params_pulse, do_plots);
        pfit(1:2, N==Ntest,TM==TMtest) = polyfit(Dtest(1:5), log(S_sech(1:5,N==Ntest, TM==TMtest)).', 1);
        beff_sech(N==Ntest,TM==TMtest) = -pfit(1, N==Ntest, TM==TMtest);
        S_sech_fit(:,N==Ntest, TM==TMtest) = exp(-beff_sech(N==Ntest,TM==TMtest) * Dtest);
    end
%     pfit(1:2, TM==TMtest) = polyfit(Ntest, bsech(:, TM==TMtest).', 1);
%     bsech_fit(:, TM==TMtest) = polyval(pfit(1:2, TM==TMtest), Nfit);
end

figure(101),  semilogy(Dtest/Snorm, S_sech(:,:,1)./ repmat(S_sech(1,:,1), [length(Dtest), 1]), ...
        Dtest/Snorm, S_sech_fit(:,:,1), '--'), ylim([1e-3 1])
figure(201),  semilogy(Dtest/Snorm, S_sech(:,:,2)./ repmat(S_sech(1,:,2), [length(Dtest), 1]), ...
        Dtest/Snorm, S_sech_fit(:,:,2), '--'), ylim([1e-3 1])
figure(301),  semilogy(Dtest/Snorm, S_sech(:,:,3)./ repmat(S_sech(1,:,3), [length(Dtest), 1]), ...
        Dtest/Snorm, S_sech_fit(:,:,3), '--'), ylim([1e-3 1])

return

%%
Nfit = beff_SLR(:, 1) * Snorm


%%
T1 = Inf; T2 = Inf;
D = 0; %1e2; 
params = [dk, TG, T1, T2, D];


%% SE
[F,k] = epg({'T', 'EDS', 'T', 'EDS'},{pi/2, params, pi,params})
bSE = (GAMMA * 2*pi *G * TG)^2 * (2/3)*TG;
exp(-bSE*D)

%% STE
params_mix = [0, TM+eps, T1, eps/1e100, D];
% [F,k] = epg({'T', 'EDS', 'T', 'EDS','T', 'EDS'},...
%     {pi/2, params, pi/2,params_mix,pi/2, params})
%STEP
[F,k] = epg({'T', 'EDS', 'T', 'EDS','T', 'EDS','T'},...
    {pi/2, params, pi/2,params_mix,pi/2, params,pi/2})
bSTE = (GAMMA * 2*pi *G * TG)^2 * (TM + (2/3)*TG);
.5*exp(-bSTE*D)



%% validate with STE
N=2; % S_nom = 0.5025; SG = 2/3; SM = 1; S0 = 0;

TG = [1e-3 2e-3 4e-3 6e-3 8e-3];TGfit = linspace(min(TG), max(TG), 40);

TM = 1e-1; %[0 1e-3 1e-2];
p = zeros(2,length(TM));

for TMtest = TM
    
    for TGtest = TG
        disp(['TG = ' num2str(TGtest) ' TM = ' num2str(TMtest)])
        bhat(TGtest==TG,TMtest==TM) = ...
            calc_b_val(N, TGtest, TMtest, do_plots);
    end
    
    p(:,TMtest==TM)  = polyfit(TG.',bhat(:,TMtest==TM),1);
    figure(find(TMtest==TM))
    plot(TG, bhat(:,TMtest==TM),'x', TGfit, polyval(p(:,TMtest==TM),TGfit), '--')
    
end




%%
D = [.7 .8 3.2] *1e-3; % WM GM CSF

b = linspace(.1, 1e4, 100).';

semilogy(b, exp(-b * D))
ylim([1e-3 1])  % good limit of 1e-3 or 1e-2 (bmax ~ 6000)

%%
N = 6; S_nom = 0.703;  SG = 5.5; SM = 1; S0 = 0;  % good
%      N = 12; S_nom = 0.767; Sb = 13.;
%      N = 24; S_nom = 0.839;  Sb = 30;

% b_ste = (GAMMA * 2*pi *G * TG)^2 * (tm + (2/3)*TG);
% A_ste = 0.5 * exp(-b_ste * D / (2*pi*GAMMA *G*TG)^2)

% this maybe decent model for 180-180 (or prob 180-90) super STE
% underestimation, esp for higher N
% b_sste = (GAMMA * 2*pi *G * TG)^2 * (tm + (N -4/3)*TG);
% A_sste = S_nom * exp(-b_sste * D / (2*pi*GAMMA *G*TG)^2)

b_sste = (GAMMA * 2*pi *G * TG)^2 * (SM*tm + SG*TG + S0);
A_sste = S_nom * exp(-b_sste * D / (2*pi*GAMMA *G*TG)^2);

%%
%disp(['sim = ' num2str(Az) '; predict = ' num2str(A_sste)])