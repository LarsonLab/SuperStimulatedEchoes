function [rf_G, tau] = adiabatic_sech(N, B10, Om, mu);
% [rf_G, tau] = adiabatic_sech(N, B10, Om, mu);
%    Creates a hyperbolic secant adiabatic RF pulse.
%    rf = B10*sech(beta*tau)*exp(i*mu*log(sech(beta*tau))),
%    where beta = Om/mu
%
%    N - number of samples
%    B10 - maximum B1 (Gauss) 
%    Om - omega controls BW of the pulse
%       Set Om = bw*pw*2*pi/1.87 for an approximate FWHM bandwidth of bw 
%    mu - adiabaticity, controls sharpness of BW transition
%
%    rf_G - adiabatic pulse, scaled to Gauss
%    tau - normalized time
%
% PEZL 7/8/04, modified format 9/14/04

% Hyperbolic Secant shaped Adiabatic inversion pulse
GAMMA = 4258; % Hz/G
tau = [-1:2/(N-1):1];

beta = Om/mu;

B1 = B10*sech(beta*tau);
pB1 = mu*log(sech(beta*tau));

rf_G = B1.*exp(i*pB1);
