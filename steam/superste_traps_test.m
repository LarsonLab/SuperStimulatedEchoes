%% SPSS

N = 1000;
w = linspace(N,-pi, pi);
Mzopt = [ones(1,N/4), -ones(1,N/2), ones(1,N/4)];
Zkopt = fft(Mzopt);
cplot(Zkopt)

% Maximize Zk, k=1,3,5,...
% target values

%%

N = 18;
tbw = N/2;
d1 = .01/8; d2 = sqrt(.01/2);

d = dinf(d1,d2);
fw = d/tbw;
finv = [0 .5*(1-fw) .5*(1+fw) 1];a = [1 1 0 0];
binv = firls(N-1, finv, a, [1 d1/d2]);
%binv(2:2:end) = fliplr(binv(2:2:end));

rfinv = b2rf(binv); % small-tip?

rfinv = real(rfinv);



train = cell(2*N-1,1); trainargs =train;

for w=1:N
    train{2*w-1} = 'T';
    trainargs{2*w-1} = rfinv(w);
    
    if w < N
        train{2*w} = 'S';
        trainargs{2*w} = round(exp(i*pi* w * 1));
    end
end
figure
[F,k] = epg(train, trainargs, [], [], 1);

