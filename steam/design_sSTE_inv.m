function [rf Nw Nwait rfinv] = design_sSTE_inv(N, dT, params)

default_params = struct('method', 'SLR', 'tbw', max(3,N/2), 'd1', .01, 'd2', .01,...
    'mu', 3, 'adiab_scale', 1.5, 'peakB1', 0.15, 'GAMMA', 4257,...
    'dt', 40e-6, 'write', 0);

if nargin < 3 || isempty(params)
    params = default_params;
else
    names = fieldnames(default_params);
    for k = 1:length(names)
        if ~isfield(params, names(k))
            params.(names{k}) = default_params.(names{k});
        end
    end
end

switch params.method
    case 'SLR'
        
        d1n = params.d1/8; d2n = sqrt(params.d2/2);
        d = dinf(d1n,d2n);
        fw = d/params.tbw;
        finv = [0 .5*(1-fw) .5*(1+fw) 1];
        a = [1 1 0 0];
        binv = firls(N-1, finv, a, [1 d1n/d2n]);
        
        rfinv = b2rf(binv);
        
        rfinv = real(rfinv);
        
    case 'sech'
        
        rfinv = adiabatic_sech(N, 1, params.tbw * 2*pi/1.87, params.mu);
        rfinv = params.adiab_scale*pi*rfinv/sum(rfinv);
        
        
        
    otherwise
        error(['method: ' params.method ' not supported'])
end


TE = dT * N;
temp_B1 = max(abs(rfinv)) / (2*pi*params.GAMMA*params.dt);
Nw = ceil(temp_B1/params.peakB1);

Nwait = round((TE/ params.dt - (N*Nw+1)) / N);

rf = zeros(1, Nwait*(N-1) + N*Nw);

for w = 1:N
    rf((1:Nw) + (w-1)*(Nwait+Nw)) = rfinv(w)/Nw;
end

if params.write
    nompw = length(rf)*params.dt;
    rfwrite(rfscaleg(rf, nompw*1e3,params.GAMMA), nompw, pi, params.GAMMA, nompw/2, 0);
    disp('delta_rf1_lobes = ')
    disp((Nw+Nwait)*params.dt)
    disp('pw_rf1_lobes = ')
    disp(Nw*params.dt)
end



