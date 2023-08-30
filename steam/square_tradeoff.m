% compare loss due to increased TE with gain due to square-encoding

T2s = logspace(-3, 0, 201);
Npulse = 3:2:7; % Npulse=5 starts to look good
T1 = 1;
dt = 5e-5;

dT = 10e-3; %time between RF pulses
Nsubpulse = 10;
fsim = linspace(-1/(2*dT), 1/(2*dT), 201);

d1 = .01/8; d2 = sqrt(.01/2);
d = dinf(d1,d2);
% real filter:

Mz_stored = zeros(length(Npulse), length(T2s));

for N = Npulse
  tbw = max(floor(N/2), 3);
  fw = d/tbw;
  finv = [0 .5*(1-fw) .5*(1+fw) 1];a = [1 1 0 0];
    binv = firls(N-1, finv, a, [1 d1/d2]);
    rfinv = b2rf(binv);

    %rfinv = real(rfinv);
    rfinvg = rfscaleg(rfinv, Nsubpulse*N*dt*1e3);
    
    rfsim = [];
    for n = 1:N-1
        rfsim = [rfsim, rfinvg(n)*ones(1,Nsubpulse) zeros(1,dT/dt - Nsubpulse)];
    end
    rfsim = [rfsim, rfinvg(N)*ones(1,Nsubpulse)];
    
    for T2 = T2s
        [mx, my, mz] = bloch(rfsim, zeros(size(rfsim)), dt, T1, T2, ...
            fsim, 0, 0);
        Mz_store(find(N==Npulse),find(T2==T2s)) = sum(abs(mz));
        if 0%(T2==T2s(end))
	  plot(fsim,mz), pause
	
	end
    
    end
end

% conventional...
mz_max = sum(abs(cos(2*pi*fsim*dT)));
mz_conv = mz_max * exp(-dT./T2s);

semilogx(T2s, Mz_store, T2s, mz_conv)
%plot(T2s, Mz_store, T2s, mz_conv)
legend(int2str(Npulse.'))
axis([.010 .1 40 1.5*mz_max])
xlabel('T_2'), ylabel('Signal')