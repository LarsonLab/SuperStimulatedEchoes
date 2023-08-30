%%
addpath ~/matlab/C13_pulses/
mets = setup_C13_mets;

% Constants
%
gamma = 1070;

% Set up vector of metabolites
%
str_met = {'Pyr', 'Pyr->Pyrh', 'Pyr->Ala', 'Pyr->Lac'};
f_met = ([mets.pyr.Hz mets.pyr_hyd.Hz mets.ala.Hz mets.lac.Hz]' - mets.pyr.Hz);
% 181 391  % normal mice 10/2010, 10 Hz resolution

% Set up STEAM parameters
%
if 1
    Tmix_test = linspace(0, 10e-3, 5000);
    plot(Tmix_test*1e3, mod(f_met * Tmix_test,0.5),'-')%, ...
%        Tmix_test*1e3, mod((f_met+10) * Tmix_test,0.5),'--', [0 max(Tmix_test*1e3)], [1/4 1/4],':')
    legend(str_met)
    ylabel('\Delta \phi (\times \pi)')
    xlabel('TE/2 (ms)')
    % hard to get both with nice phase
end


%%
str_met = {'Pyr', 'Pyr->Ala', 'Pyr->Lac'};
f_met = ([mets.pyr.Hz mets.ala.Hz mets.lac.Hz]' - mets.pyr.Hz);

Tmix = (3/4 + 2) / f_met(3);  % tuned for lac phase     % Time between 90's
gamp = 2;                               % Gradient lobe amplitude (G/cm)

% Set up position sampling of signal
%
xmax = 0.1;                               % +/- 2cm
x = [-200:200]/200 * xmax;

% What is instantaneous frequency dependent on x position and metabolite
%
f_x = f_met * ones(1, length(x)) + ...
      ones(length(f_met),1) * gamp * gamma * x ;

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
A_met = [diag([1 1 1]); 1 0 0; 1 0 0];
mz_met = A_met * mz_2;

% Get new frequencies
%
f_x2 = [f_x; f_x(2,:); f_x(3,:)];

% Get complex signal after 3rd RF pulse
%
mxy_3 = mz_met .* exp(-i*2*pi*f_x2*Tmix);


figure(1)
subplot(321)
plot(x, mz_met(1,:)), title('pyr')
subplot(323)
plot(x, mz_met(3,:)), title('lac'); ylabel('M_Z')
subplot(325)
plot(x, mz_met(5,:)), title('pyr->lac'); xlabel('x')

subplot(322)
cplot(x, mxy_3(1,:)); legend('M_X', 'M_Y')
subplot(324)
cplot(x, mxy_3(3,:)); 
subplot(326)
cplot(x, mxy_3(5,:)); xlabel('x')


figure(2)
subplot(321)
plot(x, mz_met(1,:)), title('pyr')
subplot(323)
plot(x, mz_met(2,:)), title('ala'); ylabel('M_Z')
subplot(325)
plot(x, mz_met(4,:)), title('pyr->ala'); xlabel('x')

subplot(322)
cplot(x, mxy_3(1,:)); legend('M_X', 'M_Y')
subplot(324)
cplot(x, mxy_3(2,:)); 
subplot(326)
cplot(x, mxy_3(4,:)); xlabel('x')

