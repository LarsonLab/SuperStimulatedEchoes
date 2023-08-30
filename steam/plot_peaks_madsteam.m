addpath ~/matlab/scratch/steam/

expdir = '/data/vig1/C13/animal/Normal_Mice/M036_101213/';
fname = 'madsteam1';
writeflag = 1;

addpath(expdir)
eval([fname '_parameters']);

[data, spec, phi0] = process_madsteam([expdir fname], Iref,phi1,mets, Imets,TR, 1);

%extract mag, phase of peaks
Nwidth = 5;
Nt = size(spec,2);
t = (0:Nt-1) * 1;

Amets = zeros(length(Imets), Nt);
for n =1:length(Imets)
    Amets(n,:) = sum(spec(Imets(n) + (-Nwidth:Nwidth), :), 1);
    figure(n)
    plot(t,real(Amets(n,:)),'x--',t,imag(Amets(n,:)),'o--', ...
        t, abs(Amets(n,:)), ':');
    xlabel('time (s)'), legend('real', 'imag')
    title(mets{n})
    if writeflag
        print('-depsc', [expdir fname '_' mets{n}]);
    end
end

