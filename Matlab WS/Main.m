%% Clean environment
close all;
clear variables;
clc;
format shorte;

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
fs = 50e3; % Hz - Switching frequency
Lf = 22e-6; % H - Filter inductance
Cf = 470e-6; % F - Filter capacitance
Rf = 37e-3; % Ohm - Filter inductor resistance

% TFs
tfverbose = true;
tfplots = true;

% LPRS
fm = 500e3; % Hz - Maximum evaluated switching frequency
b = 5e-3; % V - Relay hysteresis amplitude, 0 to disable analysis
c = 1; % Relay output amplitude (symmetrical, fixed by H-bridge topology)
lprsverbose = true; % Whether to print to console LPRS check verbose results
lprsplots = true; % Whether to plot LPRS check results

%% Compute Low Pass LC filter component values
% Compute the cut-off as the geometric mean of grid and switching freq.
fc = sqrt(fg*fs); % Hz - Cut-off frequency

% Compute optimal C (fc = 1/(2*pi*sqrt(L*C)))
Cfo = 1/(Lf*(2*pi*fc)^2); % F - Optimal filter capacitance

%% Compute plant
G = tf(E, [Lf*Cf, Cf*Rf 1]) % TF - Contribution to Vo from u
Zo = -tf([Lf, Rf], [Lf*Cf, Cf*Rf 1]) % TF - Contribution to Vo from Io

CheckTF(G, 'G', tfverbose, tfplots);

%% Set up LPRS
farr = 0:2:fm; % Hz - Evaluate LPRS for frequencies 0 to fm Hz
warr = farr*2*pi; % Rad/s - Angular frequency equivalent

%% Compute LPRS for base plant G
[wS, iwS, KnS, bS] = CheckLPRS(G, 'G', warr, c, b, lprsverbose, lprsplots);

%% Compute PFC (Check ASPRness)

wg = 2*pi*fg; % Grid (reference) angular frequency

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
% Since Gc must be rd 1 or 0: Gc = k * (s^2 + w^2) / (a1*s^3 + a2*s^2 + a3*s + 1)

%% Compensator 1: passive Twin T notch
% Twin T notch is rd = 0, add LP single pole for rd = 1

fprintf('\n-- Compensator 1: Passive Twin T notch + LP pole --\n');

% Implements: Gc1 = tf([1 0 wg^2], [1 4*wg wg^2])*tf(1, [tau 1]);
function out = GcOptimComp1(coefs, G) % coefs = [tau]
    global wg
    Gc = tf([1 0 wg^2], [1 4*wg wg^2])*tf(1, [coefs(1) 1]);
    Ga = G + Gc;
    out = real([pole(Ga); zero(Ga)]);
end

coefs1 = 1;
coefs1 = FindCompensator(@(x)GcOptimComp1(x, G), coefs1, 0)

Gc1 = tf([1 0 wg^2], [1 4*wg wg^2])*tf(1, [coefs1(1) 1]);
Ga1 = G + Gc1;

CheckTF(Gc1, 'Gc1', true, true);
CheckTF(Ga1, 'Ga1', true, true);

CheckLPRS(Ga1, 'Ga1', warr, c, b, lprsverbose, lprsplots);

%% Compensator 2: active Twin T notch
% Twin T notch is rd = 0, add LP single pole for rd = 1

fprintf('\n-- Compensator 2: Active Twin T notch + LP pole --\n');

% Implements: Gc2 = tf([1 0 wg^2], [1 wg/Q wg^2])*tf(1, [tau 1]);
function out = GcOptimComp2(coefs, G) % coefs = [Q, tau]
    global wg
    Gc = tf([1 0 wg^2], [1 wg/coefs(1) wg^2])*tf(1, [coefs(2) 1]);
    Ga = G + Gc;
    out = real([pole(Ga); zero(Ga)]);
end

coefs2 = [10, 1];
coefs2 = FindCompensator(@(x)max(GcOptimComp2(x, G)), coefs2, [-inf eps])

Gc2 = tf([1 0 wg^2], [1 wg/coefs2(1) wg^2])*tf(1, [coefs2(2) 1]);
Ga2 = G + Gc2;

CheckTF(Gc2, 'Gc2', true, true);
CheckTF(Ga2, 'Ga2', true, true);

% Solution not fit for LPRS calculations
%0

%% Compensator 3: active Bainter notch
% Bainter notch is rd = 0, add LP single pole for rd = 1

fprintf('\n-- Compensator 3: Active Bainter notch + LP pole --\n');

% Gc3 = tf([1 0 wg^2], [1 wo/Q wo^2])*tf(1, [tau 1]);
function out = GcOptimComp3(coefs, G) % coefs = [wo, Q, tau]
    global wg
    Gc = tf([1 0 wg^2], [1 coefs(1)/coefs(2) coefs(1)^2])*tf(1, [coefs(3) 1]);
    Ga = G + Gc;
    out = real([pole(Ga); zero(Ga)]);
end

coefs3 = [100, 0.707, 0.5];
coefs3 = FindCompensator(@(x)max(GcOptimComp3(x, G)), coefs3, [eps eps eps])

Gc3 = tf([1 0 wg^2], [1 coefs3(1)/coefs3(2) coefs3(1)^2])*tf(1, [coefs3(3) 1]);
Ga3 = G + Gc3;

CheckTF(Gc3, 'Gc3', true, true);
CheckTF(Ga3, 'Ga3', true, true);

CheckLPRS(Ga3, 'Ga3', warr, c, b, lprsverbose, lprsplots);

%% Compensator 3.2: active Bainter notch
% Bainter notch is rd = 0, add LP single pole for rd = 1

fprintf('\n-- Compensator 3.2: Active Bainter notch + LP pole --\n');

% Gc3.2 = H*tf([1 0 wg^2], [1 wo/Q wo^2])*tf(1, [tau 1]);
function out = GcOptimComp32(coefs, G) % coefs = [wo, Q, tau, H]
    global wg
    Gc = coefs(4)*tf([1 0 wg^2], [1 coefs(1)/coefs(2) coefs(1)^2])*tf(1, [coefs(3) 1]);
    Ga = G + Gc;
    out = real([pole(Ga); zero(Ga)]);
end

coefs32 = [100, 0.5, 0.5, 1];
coefs32 = FindCompensator(@(x)max(GcOptimComp32(x, G)), coefs32, [eps eps eps eps])

Gc32 = coefs32(4)*tf([1 0 wg^2], [1 coefs32(1)/coefs32(2) coefs32(1)^2])*tf(1, [coefs32(3) 1]);
Ga32 = G + Gc32;

CheckTF(Gc32, 'Gc32', true, true);
CheckTF(Ga32, 'Ga32', true, true);

CheckLPRS(Ga32, 'Ga32', warr, c, b, lprsverbose, lprsplots);

%% Compensator 4: active Boctor notch
% Boctor notch is rd = 0, add LP single pole for rd = 1

fprintf('\n-- Compensator 4: Active Boctor notch + LP pole --\n');

% Implements the same TF as compensator 3 (but more restricted)
% Implements: Gc4 = tf([1 0 wg^2], [1 wo/Q wo^2])*tf(1, [tau 1]);
function out = GcOptimComp4(coefs, G) % coefs = [wo, Q, tau]
    global wg
    Gc = tf([1 0 wg^2], [1 coefs(1)/coefs(2) coefs(1)^2])*tf(1, [coefs(3) 1]);
    Ga = G + Gc;
    out = real([pole(Ga); zero(Ga)]);
end

coefs4 = [1, 0.707, 0.1];
coefs4 = FindCompensator(@(x)max(GcOptimComp4(x, G)), coefs4, [eps eps eps], [wg inf inf])

Gc4 = tf([1 0 wg^2], [1 coefs4(1)/coefs4(2) coefs4(1)^2])*tf(1, [coefs4(3) 1]);
Ga4 = G + Gc4;

CheckTF(Gc4, 'Gc4', true, true);
CheckTF(Ga4, 'Ga4', true, true);

CheckLPRS(Ga4, 'Ga4', warr, c, b, lprsverbose, lprsplots);