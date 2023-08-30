
% Consult with Albert on frequency spec for pulse
% 
lac = 158;
pyr_hyd = 36;
ala = -53;
pyr = -235;
bicarb = -562;

% Constants
%
gamma = 1070;

% Set up vector of metabolites
%
str_met = {'Pyr', 'Ala', 'Lac', 'Pyr->Ala', 'Pyr->Lac'};
f_met = ([pyr ala lac]' - lac);

f_diff = [ala-pyr, lac-pyr];

t_phased  = 1./f_diff;

% Tmix: (for full echo)
% 66.3 is great for getting additive lac/ala that is created
% 73.2 gives nice phase (pyr->lac is at +90, pyr->ala is at -90)
% 70.6 gives pyr->lac/ala at +90

% Set up STEAM parameters
%
Tmix = 70.6e-3; %t_phased(2)*27.75;                           % Time between 90's (TE/2)
gamp = 1.1;                               % Gradient lobe amplitude (G/cm)
Tg = 2e-3;                                 % Gradient duration (s)
T1 = 30;                                % T1 of Pyruvate

% Set up position sampling of signal
%
xmax = 1;                               % +/- xmax cm
x = [-200:200]/200 * xmax;

% What is instantaneous frequency dependent on x position and metabolite
%
f_x = f_met * ones(1, length(x)) + ...
      ones(length(f_met),1) * (gamp*Tg/Tmix) * gamma * x ;

% Complex signal after mixing time... 
%
mxy_1 = exp(-i*2*pi*f_x*Tmix);

% Longitudinal mag after 2nd RF pulse
%
mz_2 = real(mxy_1);

% Allow for pyr->ala and pyr->lac metabolism
%
met_factor = 0.1;
A_met = [diag([(1-met_factor) 1 1]); met_factor 0 0; met_factor 0 0];
mz_met = A_met * mz_2;

% Get new frequencies
%
f_x2 = [f_x; f_x(2,:); f_x(3,:)];

% Get complex signal after 3rd RF pulse
%
mxy_3 = mz_met .* exp(-i*2*pi*f_x2*Tmix);

mxy_3_sum = mxy_3(1:3,:) + [zeros(1,length(x)); mxy_3(4,:); mxy_3(5,:)];

%Tmixes = ; 
subplot(221)
plot(x, real(mxy_3));
subplot(222)
plot(x, imag(mxy_3));
legend(str_met)

subplot(223)
plot(x, real(mxy_3_sum));
subplot(224)
plot(x, imag(mxy_3_sum));
legend(str_met(1:3))
