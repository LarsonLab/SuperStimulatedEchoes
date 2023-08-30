function [F,k] = epg(train, trainargs, F0, k0, do_plots)
% [F,k] = epg(train, trainargs, F0, k0, do_plots)
%
% train | trainargs - cell arrays describing spin operations:
%   'T' | flip angle (radians) - RF
%   'E' | [t T1 T2] - relaxation decay
%   'D' | [dk t d] - diffusion
%   'S' | dk - dephasing
%   'EDS' | [dk t T1 T2 d] - combined relaxation, diffusion & dephasing
% F0, k0 - Initial state
% do_plots (optional)
%
% F = [F0 Z0 Fk1 F-k1 Zk1 Fk2 F-k2 Zk2 ...]
% k = [0 k1 k2 ...]

if nargin < 3 || isempty(F0)
    F = [0 1].'; k = 0;
else
    F = F0(:); k = k0(:);
end

if nargin < 4 || isempty(do_plots)
    do_plots = 0;
end



w = linspace(-pi, pi);

for n = 1:length(train)
    [F,k] = feval(train{n}, F,k, trainargs{n});
    if do_plots
        kmin = min(k(2:end));
        Mxy = F(1)*ones(size(w)); Mz = F(2)*ones(size(w));
        for Ik = 2:length(k)
            Mxy = Mxy + F(3*(Ik-1)) * exp(1i*k(Ik)/kmin * w);
            Mxy = Mxy + conj(F(3*(Ik-1)+1)) * exp(-1i*k(Ik)/kmin * w);
            Mz = Mz + F(3*Ik-1)* cos(k(Ik)/kmin * w)*2;
        end
        subplot(211)
        cplot(w, Mxy)
        subplot(212)
        plot(w, Mz)
        drawnow
    end
end

return


function [F,k] = EDS(F,k,args)
dk = args(1); t = args(2);
T1 = args(3); T2 = args(4); d = args(5);
if length(args) < 6
    M0 = 1;
else
    M0 = args(6);
end

if (dk == 0) && (t == 0)
    return;
end

F = E(F,k, [t,T1,T2, M0]);
F = D(F,k, [dk, t, d]);
[F,k] = S(F,k,dk);

return


function [Fnew,knew] = S(F,k, dk)

if dk == 0
    Fnew = F; knew = k;
    return;
end

knew = k;
Fnew = zeros(size(F));
Fnew(2:3:end) = F(2:3:end); % Z states preserved

% k=0
In = find(abs(dk)==k);
if isempty(In)
    knew(end+1) = abs(dk);
    if dk > 0
        Fnew(end+[1:3]) = [F(1),0,0];
    else
        Fnew(end+[1:3]) = [0,conj(F(1)),0];
    end
else
    if dk > 0
        Fnew( (In-1)*3 ) =  Fnew( (In-1)*3 ) +F(1);
    else
        Fnew( (In-1)*3+1 ) =  Fnew( (In-1)*3+1 ) +conj(F(1));
    end
end


for n = 2:length(k) % k > 0
    % Fk -> Fk+dk
    if k(n) + dk == 0
        Fnew(1) = Fnew(1)+F(3*(n-1));
    else
        In = find(abs(k(n)+dk) == knew);
        if isempty(In)
            knew(end+1) = abs(k(n)+dk);
            if k(n) + dk > 0
                Fnew(end+[1:3]) = [F(3*(n-1)),0,0];
            else
                Fnew(end+[1:3]) = [0,conj(F(3*(n-1))),0];
            end
        else
            if k(n) + dk > 0
                Fnew( (In-1)*3 ) =  Fnew( (In-1)*3 ) +F(3*(n-1));
            else
                Fnew( (In-1)*3+1 ) =  Fnew( (In-1)*3+1 ) +conj(F(3*(n-1)));
            end
            
        end
    end
    
    
    %        knew(end+1) = k(n)+dk;        F(end+[1:3]) = [F(3*(n-1)),0,0];
    F(3*(n-1)) = 0;
    
    % F-k -> F-k+dk
    if -k(n) + dk == 0
        Fnew(1) = Fnew(1)+conj(F(3*(n-1)+1));
    else
        In = find(abs(-k(n)+dk) == knew);
        if isempty(In)
            knew(end+1) = abs(-k(n)+dk);
            if -k(n) + dk > 0
                Fnew(end+[1:3]) = [conj(F(3*(n-1)+1)),0,0];
            else
                Fnew(end+[1:3]) = [0,F(3*(n-1)+1),0];
            end
        else
            if -k(n) + dk > 0
                Fnew( (In-1)*3 ) =  Fnew( (In-1)*3 ) +conj(F(3*(n-1)+1));
            else
                Fnew( (In-1)*3+1 ) =  Fnew( (In-1)*3+1 ) +F(3*(n-1)+1);
            end
            
        end
    end
end

return

function [F,k] = T(F,k,rf)
a =  abs(rf); P = angle(rf);

 T0 = [cos(a/2)^2               exp(2*1i*P)*sin(a/2)^2      -1i*exp(1i*P)*sin(a); ...
      -1i/2*exp(-1i*P)*sin(a)   1i/2*exp(1i*P)*sin(a)       cos(a)];
%T0 = [cos(a/2)^2 sin(a/2)^2 -sin(a); ...
%     -sin(a)/2   sin(a)/2   cos(a)];
F(1:2) = T0 * [F(1);conj(F(1)); F(2)];

 T1 = [cos(a/2)^2               exp(1i*2*P)*sin(a/2)^2      -1i*exp(1i*P)*sin(a); ...
       exp(-2*1i*P)*sin(a/2)^2  cos(a/2)^2                  1i*exp(-1i*P)*sin(a); ...
       -1i/2*exp(-1i*P)*sin(a)  1i/2*exp(1i*P)*sin(a)       cos(a)];
%T1 = [cos(a/2)^2 sin(a/2)^2 -sin(a); ...
%     sin(a/2)^2  cos(a/2)^2 sin(a); ...
%     -sin(a)/2   sin(a)/2   cos(a)];
for n = 2:length(k)
    F([0:2] + 3*(n-1)) = T1*F([0:2] + 3*(n-1));
end

return

function [F,k] = E(F,k,args)
t = args(1); T1 = args(2); T2 = args(3);
if length(args) < 4
    M0 = 1;
else
    M0 = args(4);
end

% recover to M0
F(1:2) = [exp(-t/T2) 0; 0 exp(-t/T1)] * F(1:2) + M0*[0; 1-exp(-t/T1)];

for n = 2:length(k)
    F([0:2] + 3*(n-1)) = ...
        [exp(-t/T2) 0 0; 0 exp(-t/T2) 0; 0 0 exp(-t/T1)]*F([0:2] + 3*(n-1));
end

return

function [F, k] = D(F,k,args)
dk = args(1); t = args(2); d = args(3);
% check b calculation (dk incorporation)
edk1 = exp(- ( (k + dk/2).^2 * t + dk^2*t/12 )*d);
edk2 = exp(- ( (-k + dk/2).^2 * t + dk^2*t/12 )*d);
ek = exp( - k.^2 * t * d);

F(1:2) = diag([edk1(1) 1])* F(1:2);
for n = 2:length(k)
    F([0:2] + 3*(n-1)) = ...
        diag([edk1(n),edk2(n),ek(n)]) *F([0:2] + 3*(n-1));
end

%
% F = F .* [repmat(exp(-bdk*d), [2, 1]); ...
%     exp(-bk*d)];
return
