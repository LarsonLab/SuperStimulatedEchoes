% compare loss due to increased TE with gain due to square-encoding

Npulse = [6:4:30];
% for lin phase:
% just keeps getting better with more N! for most all T2 values.  at shorter T2 values (10 ms), N>7 all converge)
% for T2 < 20 ms, N >10 begin to do worse (but not by much)
% N = 22 is near optimal for T2 = 50ms

% for max and min phase:
% much more T2 dependance
% generally better for longer T2s (>200 ms) (but this threshold decreases
% with N)


dT = 3e-3; %time between RF pulses
Nsubpulse = 10;

fsim = linspace(-1/(2*dT), 1/(2*dT), 201);
T2s = logspace(0, 2, 201) * dT;
T1 = 30;
dt = 5e-5;


d1 = .01/8; d2 = sqrt(.01/2);
d = dinf(d1,d2);

Mz_stored = zeros(length(Npulse), length(T2s));

Mz_0 = 1;
GAMMA = 1071;

for N = Npulse
    tbw = max(floor(N/2), 3);
    if 0 % max/min phase
        n2 = 2*N-1;
        d = 0.5*dinf(2*d1,0.5*d2.*d2);
        fw = d/tbw;
        finv = [0 .5*(1-fw) .5*(1+fw) 1];a = [1 1 0 0];
        binv2 = firpm(n2-1,finv,a,[1 2*d1/(0.5*d2*d2)]);
        binv = fmp(binv2);
        binv = fliplr(binv);
    else % lin phase
        d = dinf(d1,d2);
        fw = d/tbw;
        finv = [0 .5*(1-fw) .5*(1+fw) 1];a = [1 1 0 0];
        binv = firls(N-1, finv, a, [1 d1/d2]);
    end
    
    rfinv = b2rf(binv);
    
    %rfinv = real(rfinv);
    rfinvg = rfscaleg(rfinv, Nsubpulse*N*dt*1e3, GAMMA);
    
    rfsim = [];
    for n = 1:N-1
        rfsim = [rfsim, rfinvg(n)*ones(1,Nsubpulse) zeros(1,dT/dt - Nsubpulse)];
    end
    rfsim = [rfsim, rfinvg(N)*ones(1,Nsubpulse)];
    
    for T2 = T2s
        [mx, my, mz] = bloch(rfsim, zeros(size(rfsim)), dt, T1, T2, ...
            fsim, 0, 0, 0, 0, Mz_0, GAMMA);
        Mz_store(find(N==Npulse),find(T2==T2s)) = sum(abs(mz))/length(fsim);
        if 0 && (T2==T2s(end))
            plot(fsim,mz), pause
            
        end
        
    end
end

% conventional...
mz_max = sum(abs(cos(2*pi*fsim*dT)))/length(fsim);
mz_conv = mz_max * exp(-dT./T2s);

% spin-echo, minimum decay though
mz_se = exp(-dT./T2s);


%%
semilogx(dT./T2s, Mz_store, dT./T2s, mz_conv, dT./T2s, mz_se)
%plot(T2s, Mz_store, T2s, mz_conv)
legend(int2str(Npulse.'))
axis([.01 1 0 1])
xlabel('\Delta T / T_2'), ylabel('Signal')