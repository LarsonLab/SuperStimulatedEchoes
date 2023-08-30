% beff blurring in STEAM
GAMMA = 2*pi * 1071;
TR = .125;
Dtest = .2e-3;  % free water (2.3??), 0.7 Brain matter

G = [ sqrt(3)*0.11 sqrt(3)*0.35 .52 .69 .19 .19];
TM = [ 1 1 2 3 6 1];
d = [ 2.7 2.7 0.72 0.72 2.5 8.5]*1e-3;
D = TM + .01; % + TE/2

N = 8;
x = N/2; y=N/2;
count = 0;
for n=1:N-1
    dn = sign(mod(n,2)-.5);
    for in = 1:n
        acq_order(x,y) = count;
        count = count+1;
        x = x+dn;
    end
    for in = 1:n
        acq_order(x,y) = count;
        count = count+1;
        y = y+dn;
    end
end

dn = sign(mod(N,2)-.5);
for in = 1:N
    acq_order(x,y) = count;
    count = count+1;
    x = x+dn;
end

close all
for n = 1:length(G)
    
    
    b = (GAMMA * G(n) * d(n))^2 * (D(n) - d(n)/3 + acq_order * TR); % s/mm2
    
    PSF = exp(-b*Dtest);
    psf = ifft2c(PSF);
    %max(abs(psf(:)))
    
    b0(n) = b(4,4);
    bfit(n) = -log(abs(psf(5,5))/8)/Dtest;

    psffit = ifft2c(ones(size(b)) * exp(-bfit(n)*Dtest));
    
    if 0
        figure(98)
        subplot(2,2,2*n-1)
        imagesc(b)
        colorbar, axis square
        xlabel('k_x'),ylabel('k_y')
        title(['b(k_x,k_y) (s/mm^2)'])
%         
%         figure
%         imagesc(PSF)
%         colorbar, axis square
        subplot(2,2,2*n)
        imagesc(abs(psf)/8,[0 1])
        xlabel('x'),ylabel('y')
        colorbar, axis square
        title('PSF(x,y,D = 2.3e-3 mm^2/s)')
 
 
    end
    
    
    if 0
        figure(99)
        plot(1:size(psf,1), abs(psf(5,:)), 5, exp(-b(4,4)*Dtest)*8, 'x',5, exp(-bfit*Dtest)*8, 'o')
        hold on
    end
end