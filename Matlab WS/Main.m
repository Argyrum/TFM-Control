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
%disp('Poles:')
%pole(G)
disp('Zeros:');
disp(zero(G));
disp('Damping characteristics:');
damp(G);

%% Plot plant
figure;
subplot(2, 2, 1), step(G);
subplot(2, 2, 2), bode(G);
subplot(2, 2, 3), nyquist(G);
subplot(2, 2, 4), rlocusplot(G);
sgtitle('G response characteristics')

%% Compute LPRS' J(w)
f = 0:2:500000; % Hz - Evaluate LPRS for frequencies 0 to 500 kHz
w = f*2*pi; % Rad/s - Angular frequency equivalent

% Find intersections with |pi*b/(4*c)|
J = arrayfun(@(w) lprsmatr(Gss.A, Gss.B, Gss.C, w), w); % Complex - LPRS of G(w)
imJ = imag(J);

%icabove = find(1 == conv(imJ < pi*b/(4*c), [1 -1], 'same')); % Index - Crossings above -> below (index of above)
%cabove = size(icabove, 2); % Nº of crossings above

icbelow = find(1 == conv(imJ > -pi*b/(4*c), [1 -1], 'same')); % Index - Crossings below -> above (index of below)
cbelow = size(icbelow, 2); % Nº of crossings below

% Find the index(s) of the limit frequency(s)
%JS = arrayfun(@(w) TestOrbitalStability(Gss.A, Gss.B, Gss.C, w), w); % Bool - Is J(w) stable
%iwS = find(1 == conv(JS, [1 -1])) % Find transition 0 -> 1 (index of 0)

% Find the index of the minimum limit stable frequency (first value from 0)
[wS, iwS] = FindOrbitalStabilityLimit(Gss.A, Gss.B, Gss.C, w); % First w to be stable

%% Plot LPRS information

%disp('Rad/s of crossings at pi*b/(4*c):') <- Apparently, only crosses below
%for wi = w(icabove)
    %fprintf('    w = %f (kn = %f, f = %.2f Hz)\n', wi, ComputeKn(Gss.A, Gss.B, Gss.C, wi), wi/(2*pi))
%end
%disp('Rad/s of crossings at -pi*b/(4*c):') <- Superseded by FindWnPrecise
%for wi = w(icbelow)
    %fprintf('    w = %f (kn = %f, f = %.2f Hz)\n', wi, ComputeKn(Gss.A, Gss.B, Gss.C, wi), wi/(2*pi))
%end
fprintf('\nLimit frequency(s) for stability:\n')
for wi = w(iwS)
    fprintf('    w = %f (kn = %f, f = %.2f Hz)\n', wi, ComputeKn(Gss.A, Gss.B, Gss.C, wi), wi/(2*pi))
end
fprintf('\nFirst system equilibrium frequency:\n')
[~, ~, Kn] = FindWnPrecise(Gss.A, Gss.B, Gss.C, b, c, wS);

figure;

%Plot J
subplot(2, 1, 1), plot(real(J), imJ);
title("LPRS plot (Base plant)"), xlabel("Re J(w)"), ylabel("Im J(w)");
yline(-pi*b/(4*c), 'r');
legend(["J(w)", "-pi*b/(4*c)"]);

%Plot Im J vs w
JS = [false([1, iwS-1]) true([1, length(w)-(iwS-1)])]; % Stability mask for w (1 stable, 0 not)

subplot(2, 1, 2), plot(w(JS), imJ(JS), w(~JS), imJ(~JS));
%title(sprintf("Im J(w) crossing with |pi*b/(4*c)|: %d above and %d below", sum(cabove), sum(cbelow))), xlabel("w"), ylabel("Im J(w)");
title(sprintf("Im J(w) crossing with -pi*b/(4*c): %d", cbelow)), xlabel("w"), ylabel("Im J(w)");
%yline(pi*b/(4*c).*[1 -1], 'r');
yline(-pi*b/(4*c), 'r');
legend(["Im J(w) Stable", "Im J(w) Unstable", "-pi*b/(4*c)"]);
axis([0, 1e6, -1.5, pi*b/(4*c)*10]);

%% Compute PFC (Check ASPRness)

% Parallel (Ga = G + Gc)
% Hard requirement: Ga rd = 1, minimum phase
% since Ga = (nG*dGc + nGc*dG) / (dG*dGc)
% rd Ga = ord(dG*dGc) - max{ord(nG*dGc), ord(nGc*dG)} | ord(dG) = 2, ord(nG) = 0
% 1 = 2 + ord(dGc) - max{ord(dGc), 2 + ord(nGc)}
% 1 + ord(dGc) = max{ord(dGc), 2 + ord(nGc)}
% 1 = max{0, 2 + ord(nGc) - ord(dGc)}
% 1 = max{0, 2 - rd(Gc)}
% 1 = 2 - rd(Gc)
% rd(Gc) = 1
% Hard requirement: Gc rd = 1
% Hard requirement: Ga minimum phase
% Optimal: We want Ga -> G for f <= 50Hz
%   Option 1: Gc(50Hz)=0 should have a zero at f = 50 Hz
%   Option 2: Gc(<=50Hz)=0 should be high pass for f >= 50 Hz

% Optimization, fmincon, fminsearch
% Option 1: Gc(w)=0 -> Gc(w) = (s^2 + w^2)
% Since Gc must be rd 1: Gc = k * (s^2 + w^2) / (a1*s^3 + a2*s^2 + a3*s + 1)

%% Compensator by optimization

% G has 2 dof in denominator and 1 in numerator
% For full placement capabilities we need 3 dof compensator
% 3 dof + rd 1 => 2n order (2 poles + 1 zero)
% G = tf(b0, [a2 a1 1])
% Gc = tf([d1 1], [c2 c1 1])

% Case I: Restrict compensator dof to real numbers

GcD = [1 1];
GcN = 1;

fprintf('-- Compensator denominator --\n');

%options = optimset('PlotFcns',@optimplotfval, 'MaxFunEvals', 10000, 'MaxIter', 10000);
%options = optimset('Display', 'iter', 'MaxFunEvals', 10000, 'MaxIter', 10000);
options = optimset('MaxFunEvals', 10000, 'MaxIter', 10000);

function out = GcOptimDen(coefs, G)
    Gc = tf(1, [coefs 1]);
    Ga = G + Gc;
    out = max(real(pole(Ga)));
end

[GcD, ~, exitflag] = fminsearch(@(x)GcOptimDen(x, G), GcD, options);

if (exitflag ~= 1)
    fprintf('ERROR - Could not find solution (Den)\n');
else
    fprintf('OK - Solution found within Tol (Den)\n');
end

fprintf('-- Compensator numerator --\n');

%options = optimset('PlotFcns',@optimplotfval, 'MaxFunEvals', 10000, 'MaxIter', 10000);
%options = optimset('Display', 'iter', 'MaxFunEvals', 10000, 'MaxIter', 10000);
options = optimset('MaxFunEvals', 10000, 'MaxIter', 10000);

function out = GcOptimNum(coefs, G, GcD)
    Gc = tf([coefs 1], [GcD 1]);
    Ga = G + Gc;
    out = max(real(zero(Ga)));
end

[GcN, ~, exitflag] = fminsearch(@(x)GcOptimNum(x, G, GcD), GcN, options);

if (exitflag ~= 1)
    fprintf('ERROR - Could not find solution (Num)\n');
else
    fprintf('OK - Solution found within Tol (Num)\n');
end

fprintf('-- Compensator --\n');

Gc = tf([GcN 1], [GcD 1])

Ga = G + Gc;

disp('Zeros (Ga):');
disp(zero(Ga));
disp('Damping characteristics (Ga):');
damp(Ga);

%% Plot augmented plant
figure;
subplot(2, 2, 1), step(Ga);
subplot(2, 2, 2), bode(Ga);
subplot(2, 2, 3), nyquist(Ga);
subplot(2, 2, 4), rlocusplot(Ga);
sgtitle('Ga response characteristics')

%% Plot LPRS information (Ga)

% Find intersections with |pi*b/(4*c)|
Gass = ss(Ga);
Ja = arrayfun(@(w) lprsmatr(Gass.A, Gass.B, Gass.C, w), w); % Complex - LPRS of G(w)
imJa = imag(Ja);

icbelowa = find(1 == conv(imJa > -pi*b/(4*c), [1 -1], 'same')); % Index - Crossings below -> above (index of below)
cbelowa = size(icbelowa, 2); % Nº of crossings below

% Find the index of the minimum limit stable frequency (first value from 0)
[wSa, iwSa] = FindOrbitalStabilityLimit(Gass.A, Gass.B, Gass.C, w); % First w to be stable

fprintf('\nLimit frequency(s) for stability:\n')
for wi = w(iwSa)
    fprintf('    w = %f (kn = %f, f = %.2f Hz)\n', wi, ComputeKn(Gass.A, Gass.B, Gass.C, wi), wi/(2*pi))
end
fprintf('\nFirst system equilibrium frequency:\n')
FindWnPrecise(Gass.A, Gass.B, Gass.C, b, c, wSa);

figure;

%Plot J
subplot(2, 1, 1), plot(real(Ja), imJa);
title("LPRS plot (Amplified plant)"), xlabel("Re J(w)"), ylabel("Im J(w)");
yline(-pi*b/(4*c), 'r');
legend(["J(w)", "-pi*b/(4*c)"]);

%Plot Im J vs w
JSa = [false([1, iwSa-1]) true([1, length(w)-(iwSa-1)])]; % Stability mask for w (1 stable, 0 not)

subplot(2, 1, 2), plot(w(JSa), imJa(JSa), w(~JSa), imJa(~JSa));
%title(sprintf("Im J(w) crossing with |pi*b/(4*c)|: %d above and %d below", sum(cabove), sum(cbelow))), xlabel("w"), ylabel("Im J(w)");
title(sprintf("Im J(w) crossing with -pi*b/(4*c): %d", cbelowa)), xlabel("w"), ylabel("Im J(w)");
%yline(pi*b/(4*c).*[1 -1], 'r');
yline(-pi*b/(4*c), 'r');
legend(["Im J(w) Stable", "Im J(w) Unstable", "-pi*b/(4*c)"]);
axis([0, 1e6, -1.5, pi*b/(4*c)*10]);

%% Add 50 Hz notch filter (-+50 Hz zero)
wg = 2*pi*fg; % Grid (reference) angular frequency
rho = 1;

%GcF = Gc * tf([1 0 wg^2], [1 2*rho rho^2+wg^2]); % Second order notch filter
GcF = Gc * tf([1/wg 0], [1/wg 1]); % First order high pass filter
GaF = G + GcF;

disp('Zeros (GaF):');
disp(zero(GaF));
disp('Damping characteristics (GaF):');
damp(GaF);

%% Plot augmented plant
figure;
subplot(2, 2, 1), step(GaF);
subplot(2, 2, 2), bode(GaF);
subplot(2, 2, 3), nyquist(GaF);
subplot(2, 2, 4), rlocusplot(GaF);
sgtitle('GaF response characteristics')

%% Plot LPRS information (GaF)

% Find intersections with |pi*b/(4*c)|
GaFss = ss(GaF);
JaF = arrayfun(@(w) lprsmatr(GaFss.A, GaFss.B, GaFss.C, w), w); % Complex - LPRS of G(w)
imJaF = imag(JaF);

icbelowaF = find(1 == conv(imJaF > -pi*b/(4*c), [1 -1], 'same')); % Index - Crossings below -> above (index of below)
cbelowaF = size(icbelowaF, 2); % Nº of crossings below

% Find the index of the minimum limit stable frequency (first value from 0)
[wSaF, iwSaF] = FindOrbitalStabilityLimit(GaFss.A, GaFss.B, GaFss.C, w); % First w to be stable

fprintf('\nLimit frequency(s) for stability:\n')
for wi = w(iwSaF)
    fprintf('    w = %f (kn = %f, f = %.2f Hz)\n', wi, ComputeKn(GaFss.A, GaFss.B, GaFss.C, wi), wi/(2*pi))
end
fprintf('\nFirst system equilibrium frequency:\n')
FindWnPrecise(GaFss.A, GaFss.B, GaFss.C, b, c, wSaF);

figure;

%Plot J
subplot(2, 1, 1), plot(real(JaF), imJaF);
title("LPRS plot (Amplified filtered plant)"), xlabel("Re J(w)"), ylabel("Im J(w)");
yline(-pi*b/(4*c), 'r');
legend(["J(w)", "-pi*b/(4*c)"]);

%Plot Im J vs w
JSaF = [false([1, iwSaF-1]) true([1, length(w)-(iwSaF-1)])]; % Stability mask for w (1 stable, 0 not)

subplot(2, 1, 2), plot(w(JSaF), imJaF(JSaF), w(~JSaF), imJaF(~JSaF));
%title(sprintf("Im J(w) crossing with |pi*b/(4*c)|: %d above and %d below", sum(cabove), sum(cbelow))), xlabel("w"), ylabel("Im J(w)");
title(sprintf("Im J(w) crossing with -pi*b/(4*c): %d", cbelowaF)), xlabel("w"), ylabel("Im J(w)");
%yline(pi*b/(4*c).*[1 -1], 'r');
yline(-pi*b/(4*c), 'r');
legend(["Im J(w) Stable", "Im J(w) Unstable", "-pi*b/(4*c)"]);
axis([0, 1e6, -1.5, pi*b/(4*c)*10]);

%b = -imag(lprsmatr(GaFss.A, GaFss.B, GaFss.C, 238875))*4*c/pi;
return

%%
figure, bode(tf([1 0], [wg 1])), xline(wg);
figure, bode(tf([wg^2 0 0], [wg^2 2*1*wg 1])), xline(wg);
%%
wf = wg * 100;
Rtf = tf(wg, [1 0 wg^2]) + tf(1, [1 0]);
figure, impulse(Rtf), title("Sinusoid");
figure, impulse(Rtf*tf([wf 0], [wf 1])), title("1st Order");
figure, impulse(Rtf*tf([wf^2 0 0], [wf^2 2*1*wf 1])), title("2nd Order")