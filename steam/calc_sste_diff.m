function [S Dtest] = calc_sste_diff(N, dk, TG, TM, nmr_params, params_pulse, do_plots)
% super STE diffusion simulations
% specify N = [Nenc Nex] for different number of encoding and excitation
% pulses
% specify TG = [TG Trf] for including delays for RF pulses
% nmr_params include T1, T2, D, M0 (equilibrium mag, relative to 1, initial condition)
% if not specified, Dtest goes up to 1/(dk^2 * TG)

if length(N)==1
    N = [N N];
end

if length(TG)==2
    Trf = TG(2);
    TG = TG(1);
else
    Trf = 0;
end


if nargin < 5 || isempty(nmr_params)
    T1 = Inf; T2 = Inf;  M0 = 1;
    Dtest = [linspace(0, 1, 40)]/(dk^2 * TG);
else
    T1 = nmr_params.T1; T2 = nmr_params.T2;
    if isfield(nmr_params) == 'D'
        Dtest = [0 nmr_params.D];
    else
        Dtest = [linspace(0, 1, 40)]/(dk^2 * TG);
    end
    if isfield(nmr_params) == 'M0'
        M0 = nmr_params.M0;
    else
        M0 = 1;
    end
end

if nargin < 6 || isempty(params_pulse)
    params_pulse = [];
end
if nargin < 7 || isempty(do_plots)
    do_plots = 0;
end


%% parameters

%enc
if N(1)==2  %STE
    flips_enc = [90 90] * pi/180;
else
    [rftemp Nw Nwait flips_enc] = design_sSTE_inv(N(1), TG, params_pulse);
%    flips_enc = flips_enc*pi/sum(flips_enc);
end


if N(2)==2  %STE
    flips_ex = [90 90] * pi/180;
else
    [rftemp Nw Nwait flips_ex] = design_sSTE_inv(N(2), TG, params_pulse);
%    flips_enc = flips_enc*pi/sum(flips_enc);
end
if 0
    [rftemp Nw Nwait flips_ex] = design_sSTE_ex(floor(N/2)*2-1, TG);
end


for D = Dtest
    params = [dk, TG, T1, T2, D, M0];
    params_wait = [0, Trf/2, T1, T2, D, M0];
    params_mix = [0, TM+eps, T1, eps/1e100, D, M0];
    % initialize
    train = cell(1,4*(length(flips_enc) + length(flips_ex) - 2) + 3);
    trainargs = train;
    
    for n = 1:length(flips_enc)
        
        train{4*n-3} = 'T';
        trainargs{4*n-3} = flips_enc(n);
        if n ~= length(flips_enc)
           train{4*n-2} = 'EDS'; % decay, diffusion, etc
           trainargs{4*n-2} = params_wait;
           train{4*n-1} = 'EDS'; % decay, diffusion, etc
           trainargs{4*n-1} = params;
           train{4*n} = 'EDS'; % decay, diffusion, etc
           trainargs{4*n} = params_wait;
        end
    end
    train{4*(length(flips_enc)-1) + 2} = 'EDS';
    trainargs{4*(length(flips_enc)-1) + 2} = params_mix;

    for n = 1:length(flips_ex)
        In = 4*(n+length(flips_enc)-2) + 3;
        train{In} = 'T';
        trainargs{In} = flips_ex(n);
        if n ~= length(flips_ex)
            train{In+1} = 'EDS'; % decay, diffusion, etc
            trainargs{In+1} = params_wait;
            train{In+2} = 'EDS'; % decay, diffusion, etc
            trainargs{In+2} = params;
            train{In+3} = 'EDS'; % decay, diffusion, etc
            trainargs{In+3} = params_wait;
        end
    end
    
    
    [F,k] = epg(train, trainargs, [], [], do_plots);
    
    S(D==Dtest) = abs(F(2));
    
end

