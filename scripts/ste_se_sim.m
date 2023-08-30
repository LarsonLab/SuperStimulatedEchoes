% Ste vs Se parameters:
%addpath ~/matlab/steam/
Gmax = 10; %G/mm
a = 1e-3; % deadtime around pulses
writeflag = 0;

% 13C
GAMMA = 1071; atom = '13C'; % Hz/G
T1 = 30; T2 = 100e-3;  % could plot for multiple T2s

if 1
    % 1H
    GAMMA = 4258; atom = '1H';
    T1 = 1.1; T2 = 70e-3; % WM at 3T
   % T1 = 1.4; T2 = 40e-3; % skeletal muscle at 3T
   % T1 = 1.8; T2 = 100e-3; % GM at 3T
   % T1 = .8; T2 = 40e-3; % Liver at 3T
    T1 = 1.2; T2 = 30e-3; % Cartilage at 3T
end

b = logspace(0, 4,20);
%dp = (b*3/2 / (2*pi*GAMMA*Gmax)^2).^(1/3);

% spin-echo
for bb = b
    bfcn = @(d) (2*pi*GAMMA*Gmax*d)^2 * (2*d/3 + a) - bb;
    dp(bb==b) = fzero(bfcn, (bb*3/2 / (2*pi*GAMMA*Gmax)^2).^(1/3) );
end
Sp = exp(-2*(dp+a)/T2);
Npulse = [2 6 10 18 30];
SG = [2/3 5.2 9.9 19.4 31.9];
Sa = [1 5.2 9.2 17.3 27.5];
% N=12, SG = 6.1; Sa = 6.2;

%%
N=2; E=0.5;
for bb=b
    %optimize over ds, Tm including T1.
    
   % Sfn = @(d) -E *exp(-2*abs(d)/T2) * exp(- (2/3 * dp(bb==b)^3./abs(d).^2 - (N-4/3)*abs(d) ) /T1);
    Sfn = @(d) -E *exp(-2*(abs(d)+a)/T2) * exp(- (bb / (2*pi*GAMMA*Gmax*d)^2  - ((N-1)*(abs(d)+a) - abs(d)/3) ) /T1);
    [ds(N==Npulse, bb==b) Ss(N==Npulse,bb==b)] = fminsearch(Sfn, dp(bb==b));

end

%%
for N = Npulse(2:end)
    for bb = b
        [ds(N==Npulse,bb==b) Ss(N==Npulse,bb==b)] = fminsearch(...
            @(d) -sSTE_signal(d, a, bb/Gmax^2, N, SG(N==Npulse), Sa(N==Npulse),T1, T2,GAMMA), abs(ds(1,bb==b)));
    end
end
Ss = -Ss;


figure
semilogx(b, Sp, b, Ss,'--')
xlabel('b_{max} (s/mm^2)'), ylabel('S_{opt}')
legend(['SE'; int2str(Npulse.')])
ylim([0 1])

if writeflag
    print('-depsc', sprintf('~/papers/steprep/other figs/ste_vs_se_%s_T2=%d_T1=%d_a=%d.eps', ...
        atom, T2*1e3, T1, a*1e3))
end
return


% optimize, neglecting T1
ds = eps;
TM = 2/3 * dp.^3./ds.^2 - (N-1/3)*ds;
Ss = E * exp(-2*ds/T2);


N=1;
E(N) = 0.5;
%E = [.5:.1:1];
dp = logspace(-2,1,200).';

ds = zeros(length(dp), length(E));
for n = 1:length(E)
% normalized times by T2
ds(:,n) = dp + .5*log(E(n));
end
ds(ds<0) = 0;
loglog(dp, ds)

% TM normalized by T2
N=1;
TM = 2/3 * dp.^3./ds.^2 - (N-1/3)*ds;

loglog(dp, TM)

%TM/T2 >> 10 is bad for T2/T1 = 1/10 ->TM/T1 >> 1