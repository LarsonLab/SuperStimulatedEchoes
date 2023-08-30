function [rf Nw Nwait rfex] = design_sSTE_ex(N, dT, params)


default_params = struct('tbw', 3, 'd1', .01, 'd2', .01,...
    'peakB1', 0.15, 'GAMMA', 4257, 'dt', 40e-6, 'flip', pi/2);

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

d1n = params.d1/8; d2n = sqrt(params.d2/2);
d = dinf(d1n,d2n);
fw = d/params.tbw;
fex = [0 .5*(1-fw) .5*(1+fw) 1];
a = [1 1 -1 -1];


bex = firls(N-1, fex, a, [1 d1n/d2n]);
rfex = b2rf(bex*sin(params.flip/2));


TE = dT * N;
temp_B1 = max(abs(rfex)) / (2*pi*params.GAMMA*params.dt);
Nw = ceil(temp_B1/params.peakB1);

Nwait = round((TE/ params.dt - (N*Nw+1)) / N);

rf = zeros(1, Nwait*(N-1) + N*Nw);

for w = 1:N
    rf([1:Nw] + (w-1)*(Nwait+Nw)) = rfex(w)/Nw;
end



