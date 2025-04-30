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

% Filter - Calculator http://sim.okawa-denshi.jp/en/RLCtool.php
fs = 150e3; % Hz - Switching frequency
Lf = 22e-6; % H - Filter inductance
Cf = 150e-6; % F - Filter capacitance
Rf = 37e-3; % Ohm - Filter inductor resistance

% TFs
tfverbose = true;
tfplots = false;

% LPRS
fm = 1e6; % Hz - Maximum evaluated switching frequency
b = 0.1; % V - Relay hysteresis amplitude, 0 to disable analysis
c = 1; % Relay output amplitude (symmetrical, fixed by H-bridge topology)
lprsverbose = true; % Whether to print to console LPRS check verbose results
lprsplots = true; % Whether to plot LPRS check results

%% Compute Low Pass LC filter component values
% wc^2/(s^2 + 2chi wc + wc^2), LC = 1/(wc)^2, RC = 1/(2chi wc) = Q/wc

% Compute the cut-off as the geometric mean of grid and switching freq.
fc = sqrt(fg*fs); % Hz - Cut-off frequency

% Compute optimal C (fc = 1/(2*pi*sqrt(L*C)))
Cfo = 1/(Lf*(2*pi*fc)^2); % F - Optimal filter capacitance

% Compute R to avoid resonance (Rfo = Q/(C*wc) for Q = 0.707)
Rfo = 0.707*sqrt(Lf*Cf)/Cf;

%% Compute plant
G = tf(E, [Lf*Cf, Cf*Rf 1]) % TF - Contribution to Vo from u
Zo = -tf([Lf, Rf], [Lf*Cf, Cf*Rf 1]) % TF - Contribution to Vo from Io

CheckTF(G, 'G', tfverbose, tfplots);

%% Set up LPRS
farr = 0:2:fm; % Hz - Evaluate LPRS for frequencies 0 to fm Hz
warr = farr*2*pi; % Rad/s - Angular frequency equivalent

%% Compute LPRS for base plant G
[wS, iwS, KnS, bS] = CheckLPRS(G, 'G', warr, c, fs, lprsverbose, lprsplots);

%% Compute PFC (Check ASPRness)

wg = 2*pi*fg; % Grid (reference) angular frequency

% Parallel (Ga = G + Gc)
% Hard requirement: Gc rd = 1 (so Ga rd = 1)
% Hard requirement: Ga minimum phase
% Optimal: We want Ga -> G for f = fg
%   Option 1: Gc(50Hz)=0 should have a zero at f = 50 Hz
%   Option 2: Gc(<=50Hz)=0 should be high pass for f >= 50 Hz

% Option 1: Gc(w)=0 -> Gc(w) = (s^2 + w^2)
% Since Gc must be rd 1 or 0: Gc = k * (s^2 + w^2) / (s^3 + a2*s^2 + a1*s + a0)
% Gc implementable as notch filter + low-pass pole

% Twin T notch calculators
% http://sim.okawa-denshi.jp/en/TwinTCRkeisan.htm

% Single pole low-pass RC calculator
% http://sim.okawa-denshi.jp/en/CRtool.php

% Parameters
tau = 7.5e-6;
Q = 1.5;
H = 5;

%% Compensator 1: Twin T notch
% Twin T notch is rd = 0, add LP single pole for rd = 1

fprintf('\n-- Compensator 1: Twin T notch + LP pole --\n');

% Implements: Gc1 = tf([1 0 wg^2], [1 4*wg wg^2])*tf(1, [tau 1]);
coefs1 = tau; % [tau]

Gc1 = tf([1 0 wg^2], [1 4*wg wg^2])*tf(1, [coefs1(1) 1]);
Ga1 = G + Gc1;

CheckTF(Gc1, 'Gc1', tfverbose, tfplots);
CheckTF(Ga1, 'Ga1', tfverbose, tfplots);

CheckLPRS(Ga1, 'Ga1', warr, c, fs, lprsverbose, lprsplots);

%% Compensator 2: Twin T notch w/ feedback
% Twin T notch is rd = 0, add LP single pole for rd = 1

fprintf('\n-- Compensator 2: Twin T notch w/ feedback + LP pole --\n');

% Implements: Gc2 = tf([1 0 wg^2], [1 wg/Q wg^2])*tf(1, [tau 1]);
coefs2 = [tau Q]; %  [tau, Q]

Gc2 = tf([1 0 wg^2], [1 wg/coefs2(2) wg^2])*tf(1, [coefs2(1) 1]);
Ga2 = G + Gc2;

CheckTF(Gc2, 'Gc2', tfverbose, tfplots);
CheckTF(Ga2, 'Ga2', tfverbose, tfplots);

CheckLPRS(Ga2, 'Ga2', warr, c, fs, lprsverbose, lprsplots);

%% Compensator 3: Bainter notch
% Bainter notch is rd = 0, add LP single pole for rd = 1

fprintf('\n-- Compensator 3: Bainter notch + LP pole --\n');

% Gc3 = H*tf([1 0 wg^2], [1 wo/Q wo^2])*tf(1, [tau 1]);
coefs3 = [tau, Q, 5]; % [tau, Q, H]

Gc3 = coefs3(3)*tf([1 0 wg^2], [1 wg/coefs3(2) wg^2])*tf(1, [coefs3(1) 1]);
Ga3 = G + Gc3;

CheckTF(Gc3, 'Gc3', tfverbose, tfplots);
CheckTF(Ga3, 'Ga3', tfverbose, tfplots);

CheckLPRS(Ga3, 'Ga3', warr, c, fs, lprsverbose, lprsplots);

%% Compensator 4: Boctor notch
% Boctor notch is rd = 0, add LP single pole for rd = 1

fprintf('\n-- Compensator 4: Boctor notch + LP pole --\n');

% Implements: Gc4 = tf([1 0 wg^2], [1 wo/Q wo^2])*tf(1, [tau 1]);
coefs4 = [tau, Q, wg/2]; % [tau, Q, wo]

Gc4 = tf([1 0 wg^2], [1 coefs4(3)/coefs4(2) coefs4(3)^2])*tf(1, [coefs4(1) 1]);
Ga4 = G + Gc4;

CheckTF(Gc4, 'Gc4', tfverbose, tfplots);
CheckTF(Ga4, 'Ga4', tfverbose, tfplots);

CheckLPRS(Ga4, 'Ga4', warr, c, fs, lprsverbose, lprsplots);