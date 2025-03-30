%% Clean environment
close all;
clear variables;
clc;

%% Set up workspace path
if ~contains(path,"LPRS")
    addpath("Custom scripts\", "LPRS scripts\", "Simulink models\");
end

%% Input parameters
% Supply (DC infeed)
E = 400; % V - DC supply voltage

% Reference (Output to grid)
fg = 50; % Hz - Reference frequency
vr = 230; % V - Reference RMS voltage

% Filter
fs = 20e3; % Hz - Switching frequency
Lf = 330e-6; % H - Filter inductance
Rf = 75e-3; % Ohm - Filter inductor resistance

% Relay
b = 0.01; % Relay input hysteresis amplitude (symmetrical)
c = 1; % Relay output amplitude (symmetrical, fixed by H-bridge topology)

%% Compute Low Pass LC filter component values
% Compute the cut-off as the geometric mean of grid and switching freq.
fc = sqrt(fg*fs); % Hz - Cut-off frequency

% Compute C (fc = 1/(2*pi*sqrt(L*C)))
% Alt: https://electronics.stackexchange.com/questions/695465/digikey-and-others-are-incorrect-about-cutoff-frequency-for-an-lc-filter
Cf = (2*pi*fc)^-2/Lf; % F - Filter capacitance


%% Compute plant
G = tf(E, [Lf*Cf, Cf*Rf 1]) % TF - Contribution to Vo from u
Gss = ss(G); %% SS - State-Space description of G
Zo = -tf([Lf, Rf], [Lf*Cf, Cf*Rf 1]) % TF - Contribution to Vo from Io

% Analyse plant
%isproper(G)
pole(G)

%% Plot plant
figure;
subplot(2, 2, 1), step(0.5*G); % Response for a 0.2 amplitude step action
subplot(2, 2, 2), bode(G);
subplot(2, 2, 3), nyquist(G);
subplot(2, 2, 4), rlocusplot(G);

%% Compute LPRS' J(w)
f = 0:1:500000; % Hz - Evaluate LPRS for frequencies 0 to 500 kHz
w = f*2*pi; % Rad/s - Angular frequency equivalent

% Find intersections with |pi*b/(4*c)|
J = arrayfun(@(w) lprsmatr(Gss.A, Gss.B, Gss.C, w), w); % Complex - LPRS of G(w)
imJ = imag(J);

icabove = find(1 == conv(imJ < pi*b/(4*c), [1 -1], 'same')); % Indeces - Crossings above -> below (index of above)
cabove = size(icabove, 2); % Nº of crossings above

icbelow = find(1 == conv(imJ > -pi*b/(4*c), [1 -1], 'same')); % index - Crossings below -> above (index of below)
cbelow = size(icbelow, 2); % Nº of crossings below

% Find the index(s) of limit frequency(s) and gain(s) needed for stability
JS = arrayfun(@(w) TestOrbitalStability(Gss.A, Gss.B, Gss.C, w), w); % Bool - Is J(w) stable
iwS = find(1 == conv(JS > 0, [1 -1], 'same')); % Find transition 0 -> 1 (index of 0)

%% Plot LPRS information

disp('Rad/s of crossings at pi*b/(4*c):')
for wi = w(icabove)
    fprintf('    w = %f (kn = %f, f = %.2f Hz)\n', wi, ComputeKn(Gss.A, Gss.B, Gss.C, wi), wi/(2*pi))
end
disp('Rad/s of crossings at -pi*b/(4*c):')
for wi = w(icbelow)
    fprintf('    w = %f (kn = %f, f = %.2f Hz)\n', wi, ComputeKn(Gss.A, Gss.B, Gss.C, wi), wi/(2*pi))
end
disp('Edge frequency(s) for stability:')
for wi = w(iwS)
    fprintf('    w = %f (kn = %f, f = %.2f Hz)\n', wi, ComputeKn(Gss.A, Gss.B, Gss.C, wi), wi/(2*pi))
end

figure;

%Plot J
subplot(2, 1, 1), plot(real(J), imJ);
title("LPRS plot"), xlabel("Re J(w)"), ylabel("Im J(w)");
yline(pi*b/(4*c).*[1 -1], 'r');
legend(["J(w)", "pi*b/(4*c)"])

%Plot Im J vs w
subplot(2, 1, 2), plot(w(JS), imJ(JS), w(~JS), imJ(~JS));
title(sprintf("Im J(w) crossing with |pi*b/(4*c)|: %d above and %d below", sum(cabove), sum(cbelow))), xlabel("w"), ylabel("Im J(w)");
yline(pi*b/(4*c).*[1 -1], 'r');
legend(["Im J(w) Stable", "Im J(w) Unstable", "pi*b/(4*c)"])
axis([0, 1e5, -1.5, pi*b/(4*c)*10]);

%% Compute PFC (Check ASPRness)