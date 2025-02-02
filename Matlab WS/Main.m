%% Input parameters
% Grid
E = 400; % V - DC supply voltage

% Filter
fg = 50; % Hz - Output (grid) frequency
fs = 20e3; % Hz - Switching frequency
L = 330e-6; % H - LC LP^filter Inductance
Rl = 75e-3; % Ohm - Inductor resistance


%% Compute LC LP filter component values
fc = sqrt(fg*fs); % Hz - Cut-off frequency

% Compute C (fc = 1/(2*pi*sqrt(L*C)))
C = (2*pi*fc)^-2/L; % F - Capacitance C


%% Compute plant
G = tf(E, [L*C, C*Rl 1]); % TF - Contribution to Vo from u
Zo = -tf([L, Rl], [L*C, C*Rl 1]); % TF - Contribution to Vo from Io