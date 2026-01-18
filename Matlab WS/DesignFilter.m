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
Rf = 30e-3; % Ohm - Filter inductor resistance

% LPRS
c = 1; % Relay output amplitude (symmetrical, fixed by H-bridge topology)

%% Compute Low Pass LC filter component values

% Generate possible filter ratios (r/100 = L/C)
r = 1:100;

% Compute the cut-off as the geometric mean of grid and switching freq.
fc = sqrt(fg*fs); % Hz - Cut-off frequency

% Compute possible L and C pairs (1/(2*pi*fc)^2 = L*C)
LC = 1/(2*pi*fc)^2;
Carr = 10*sqrt(LC./r);
Larr = Carr.*r./100;

clear LC fc

% Plot results
figure;
tiledlayout('vertical');

nexttile;
hold on;
semilogy(r, Larr, 'DisplayName', 'L');
semilogy(r, Carr, 'DisplayName', 'C');
hold off;
title('LP component values');
xlabel('Ratio L/C (%)');
legend;
yscale log;
grid on;

%% Check b and kn values across pairs with LPRS (with G)
barr = zeros(size(r));
karr = zeros(size(r));

for i = r    
    % Compute plant
    Lf = Larr(i);
    Cf = Carr(i);
    G = tf(E, [Lf*Cf, Cf*Rf 1]); % TF - Contribution to Vo from u
    
    % Check LPRS for stability at fs
    wn = fs*2*pi; % Rad/s - Angular frequency equivalent
    Gss = ss(G);

    [barr(i), ~, karr(i), ~] = Compute_b(Gss.A, Gss.B, Gss.C, c, wn, false);
end

% Plot results
nexttile;
hold on;
yyaxis left;
semilogy(r, barr, 'DisplayName', 'b');
yyaxis right;
semilogy(r, karr, 'DisplayName', 'kn');
hold off;
title('LPRS results for G(s)');
xlabel('Ratio L/C (%)');
legend;
yscale log;
grid on;

return

%% Check b and kn values across pairs with LPRS (with G + PFC)
b2arr = zeros(size(r));
k2arr = zeros(size(r));

for i = r    
    % Compute plant
    Lf = Larr(i);
    Cf = Carr(i);
    G = tf(E, [Lf*Cf, Cf*Rf 1]); % TF - Contribution to Vo from u

    % Compute PFC (Check ASPRness)
    wg = 2*pi*fg; % Grid (reference) angular frequency
    tau = 7.5e-6;
    Q = 1.5;

    % Compensator 1: Twin T notch
    % coefs = tau;    
    % Gc = tf([1 0 wg^2], [1 4*wg wg^2])*tf(1, [coefs(1) 1]);
    % Ga = G + Gc;

    % Compensator 2: Twin T notch w/ feedback
    % coefs = [tau Q]; %  [tau, Q]
    % Gc = tf([1 0 wg^2], [1 wg/coefs(2) wg^2])*tf(1, [coefs(1) 1]);
    % Ga = G + Gc;

    % Compensator 3: Bainter notch
    coefs = [tau, Q, 5]; % [tau, Q, H]
    Gc = coefs(3)*tf([1 0 wg^2], [1 wg/coefs(2) wg^2])*tf(1, [coefs(1) 1]);
    Ga = G + Gc;
    
    % Compensator 4: Boctor notch
    % coefs = [tau, Q, wg/2]; % [tau, Q, wo]
    % Gc = tf([1 0 wg^2], [1 coefs(3)/coefs(2) coefs(3)^2])*tf(1, [coefs(1) 1]);
    % Ga = G + Gc;
    
    % Check LPRS for stability at fs
    wn = fs*2*pi; % Rad/s - Angular frequency equivalent
    Gss = ss(Ga);

    [b2arr(i), ~, k2arr(i), ~] = Compute_b(Gss.A, Gss.B, Gss.C, c, wn, false);
end

clear fm fmstep farr coefs

% Plot results
nexttile;
hold on;
yyaxis left;
semilogy(r, b2arr, 'DisplayName', 'b');
yyaxis right;
semilogy(r, k2arr, 'DisplayName', 'kn');
hold off;
title('LPRS results for Ga(s)');
xlabel('Ratio L/C (%)');
legend;
yscale log;
grid on;