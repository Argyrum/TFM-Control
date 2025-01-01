%% Input parameters
fg = 50; % Hz - Output (grid) frequency
fs = 20e3; % Hz - Switching frequency
L = 300e-3; % H - LC LP^filter Inductance

%% Compute LC LP filter component values
fc = sqrt(fg*fs); % Hz - Cut-off frequency

% Compute C (fc = 1/(2*pi*sqrt(L*C)))
C = (2*pi*fc)^-2/L; % F - Capacitance C
