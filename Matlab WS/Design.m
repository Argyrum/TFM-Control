% 1. Create model of the plant
% 2. Design a PFC
% 3. Implement PFC
% 4. Implement SSR
% 5. Implement error amplifier

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
Vg = 230*sqrt(2); % V - Reference maximum voltage

% Filter - Calculator http://sim.okawa-denshi.jp/en/RLCtool.php
fs = 150e3; % Hz - Switching frequency
Lf = 37e-6; % H - Filter inductance
Cf = 90e-6; % F - Filter capacitance
Rf = 3.2e-3; % Ohm - Filter inductor resistance

% TFs
tfverbose = true;
tfplots = false;

% LPRS
fm = 1e6; % Hz - Maximum evaluated switching frequency
fmstep = 2; % Hz - Evaluated switching frequency step
c = 1; % Relay output amplitude (symmetrical, fixed by H-bridge topology)
lprsverbose = true; % Whether to print to console LPRS check verbose results
lprsplots = false; % Whether to plot LPRS check results

% Implementation
Vc = 12; % V - OpAmp maximum output voltage

%% Compute Low Pass LC filter component values
% wc^2/(s^2 + 2chi wc + wc^2), LC = 1/(wc)^2, RC = 1/(2chi wc) = Q/wc

% Compute the cut-off as the geometric mean of grid and switching freq.
fc = sqrt(fg*fs); % Hz - Cut-off frequency

% Compute optimal C (fc = 1/(2*pi*sqrt(L*C)))
Cfo = 1/(Lf*(2*pi*fc)^2); % F - Optimal filter capacitance
fprintf('\nCf offset: %f F\n', Cf-Cfo);

% Compute R to avoid resonance (Rfo = Q/(C*wc) for Q = 0.707)
Rfo = 0.707*sqrt(Lf*Cf);
fprintf('\nRf offset: %f Ohm\n', Rf-Rfo);

clear fc

%% Compute plant
G = tf(E, [Lf*Cf, Cf*Rf 1]) % TF - Contribution to Vo from u
Zo = -tf([Lf, Rf], [Lf*Cf, Cf*Rf 1]) % TF - Contribution to Vo from Io

CheckTF(G, 'G', tfverbose, tfplots);

%% Set up LPRS
farr = 0:fmstep:fm; % Hz - Evaluate LPRS for frequencies 0 to fm Hz
warr = farr*2*pi; % Rad/s - Angular frequency equivalent

clear fm fmstep farr

%% Compute LPRS for base plant G
[~, ~, KnS, bS] = CheckLPRS(G, 'G', warr, c, fs, lprsverbose, lprsplots);

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

% Parameters
tau = 7.5e-6;
Q = 1;
H = 5;

%% Compensator 1: Twin T notch
% Twin T notch is rd = 0, add LP single pole for rd = 1
if true
    fprintf('\n-- Compensator 1: Twin T notch + LP pole --\n');
    
    % Implements: Gc1 = tf([1 0 wg^2], [1 4*wg wg^2])*tf(1, [tau 1]);
    coefs1 = tau; % [tau]
    
    Gc1 = tf([1 0 wg^2], [1 4*wg wg^2])*tf(1, [coefs1(1) 1]);
    Ga1 = G + Gc1;
    
    CheckTF(Gc1, 'Gc1', tfverbose, tfplots);
    CheckTF(Ga1, 'Ga1', tfverbose, tfplots);
    
    [~, ~, ~, ~, bE, KnE] = CheckLPRS(Ga1, 'Ga1', warr, c, fs, lprsverbose, lprsplots);
    
    clear coefs1
end

%% Compensator 2: Twin T notch w/ feedback
% Twin T notch is rd = 0, add LP single pole for rd = 1
if false
    fprintf('\n-- Compensator 2: Twin T notch w/ feedback + LP pole --\n');
    
    % Implements: Gc2 = tf([1 0 wg^2], [1 wg/Q wg^2])*tf(1, [tau 1]);
    coefs2 = [tau Q]; %  [tau, Q]
    
    Gc2 = tf([1 0 wg^2], [1 wg/coefs2(2) wg^2])*tf(1, [coefs2(1) 1]);
    Ga2 = G + Gc2;
    
    CheckTF(Gc2, 'Gc2', tfverbose, tfplots);
    CheckTF(Ga2, 'Ga2', tfverbose, tfplots);
    
    CheckLPRS(Ga2, 'Ga2', warr, c, fs, lprsverbose, lprsplots);
    
    clear coefs2
end

%% Compensator 3: Bainter notch
% Bainter notch is rd = 0, add LP single pole for rd = 1
if true
    fprintf('\n-- Compensator 3: Bainter notch + LP pole --\n');
    
    % Gc3 = H*tf([1 0 wg^2], [1 wo/Q wo^2])*tf(1, [tau 1]);
    coefs3 = [tau, Q, H]; % [tau, Q, H]
    
    Gc3 = coefs3(3)*tf([1 0 wg^2], [1 wg/coefs3(2) wg^2])*tf(1, [coefs3(1) 1]);
    Ga3 = G + Gc3;
    
    CheckTF(Gc3, 'Gc3', tfverbose, tfplots);
    CheckTF(Ga3, 'Ga3', tfverbose, tfplots);
    
    CheckLPRS(Ga3, 'Ga3', warr, c, fs, lprsverbose, lprsplots);
    
    clear coefs3
end

%% Compensator 4: Boctor notch
% Boctor notch is rd = 0, add LP single pole for rd = 1
if false
    fprintf('\n-- Compensator 4: Boctor notch + LP pole --\n');
    
    % Implements: Gc4 = tf([1 0 wg^2], [1 wo/Q wo^2])*tf(1, [tau 1]);
    coefs4 = [tau, Q, wg/2]; % [tau, Q, wo]
    
    Gc4 = tf([1 0 wg^2], [1 coefs4(3)/coefs4(2) coefs4(3)^2])*tf(1, [coefs4(1) 1]);
    Ga4 = G + Gc4;
    
    CheckTF(Gc4, 'Gc4', tfverbose, tfplots);
    CheckTF(Ga4, 'Ga4', tfverbose, tfplots);
    
    CheckLPRS(Ga4, 'Ga4', warr, c, fs, lprsverbose, lprsplots);
    
    clear coefs4
end

%% Compensator implementation
% Compensator 1 chosen
% Twin T notch - http://sim.okawa-denshi.jp/en/TwinTCRkeisan.htm
% RC low-pass - http://sim.okawa-denshi.jp/en/CRtool.php

Rc1 = 47e3;
Rc2 = 47e3;
Rc3 = 12e3;

Cc1 = 100e-9;
Cc2 = 100e-9;
Cc3 = 100e-9;

Rclp = 750;

Cclp = 10e-9;

Rc = mean([Rc1, Rc2, 2*Rc3]);
Cc = mean([Cc1, Cc2, Cc3/2]);

Gc = tf([1 0 (1/(Rc*Cc))^2], [1 4/(Rc*Cc) (1/(Rc*Cc))^2]) * tf(1, [Rclp*Cclp 1]);

figure, hold on;
bode(Gc1, Gc), xline(wg);
title("Compensator realization"), legend(["Designed", "Realization", "50 Hz"]);

clear Rc Cc

%% Cleanup compensator design stage
clear tau Q H

%% Relay realization
% Schmitt trigger

% bE/Vc = Rst1/Rst2
Rst1 = 1e3; % Ohm
Rst2 = 56e3; % Ohm

fprintf('\nRelay hysteresis offset : %.2f mV (actual %.2f)\n', 1e3*(Vc*Rst1/Rst2-bE), bE*1e3);
fprintf('\nRst2 needed offset: %.f Ohm (%.f)\n', Rst2-Vc*Rst1/bE, Rst2);

%% Error calculation
% Inverting summing amplifier

% Vout = -sum(Vin*Rout/Rin)
Rin1 = 27e3; % Ohm - Output feedback
Rin2 = 27e3; % Ohm - Reference
Rin3 = 12e3; % Ohm - Compensator feedback

Rout = 1e3; % Ohm

fprintf('\nRin1 needed offset: %.f Ohm (actual %.f)\n', Rin1-Rout*Vg/Vc, Rin1);
fprintf('\nRin2 needed offset: %.f Ohm (actual %.f)\n', Rin2-Rout*Vg/Vc, Rin2);
fprintf('\nRin3 needed offset: %.f Ohm (actual %.f)\n', Rin3-Rout*Vc, Rin3);

%% Prepare for analysis

clear bE bS c G Ga1 Ga2 Ga3 Ga4 Gc Gc1 Gc2 Gc3 Gc4 KnE KnS lprsplots lprsverbose tfplots tfverbose warr wg Zo