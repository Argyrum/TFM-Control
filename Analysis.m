%% Clean environment
close all;
clear variables;
clc;
format shorte;

%% Set up workspace path
if ~contains(path,"Data\")
    addpath("Simulink models\", "Data\", "Custom scripts\");
end

%% Input parameters and component values
load Data\SimulationWorkspace.mat

%% Configure simulation and analysis
% Simulation
ncycles_end = 8; % Number of electrical cycles to simulate
outDec = 1; % Decimation factor for simulation out values

% Analysis
ncycles_skip = ncycles_end-3; % Number of electrical cycles to skip at start
%ncycles_skip = 2; % Long simulation case (used in sweeps)
%ncycles_skip = 0; % Long simulation case (used in sweeps)
ncycles_max = ncycles_end; % Final electrical cycle for statistics and plots
fResampling = 5e6; % Uniform time resample frequency (used in THD)

%% Calculate timings
assert(ncycles_max > ncycles_skip);
tend = ncycles_end/fg; % Simulation end time
tskip = ncycles_skip/fg; % Analysis start time
tmax = ncycles_max/fg; % Analisys end time

%% Open model
%open_system("CLRFSPFCRealizedSimscape");

%% Scenarios (Uncomment to select)
% Static cases
%label = 'Noload';      % 0. Response to No load
%label = 'Rload';       % 1. Response to R 22 Ohm load
%label = 'Lload';       % 2. Response to L 1 Ohm 400 mH load
%label = 'CPLload';     % 3. Response to CPL Po 2500 W Qo 500 W
%label = 'LRectRCload'; % 4. Response to L 22 uH + rectified R 22 Ohm C 100 uF load
% Dynamic cases
%label = 'Rstepsload';  % 5. Response to R 10k -> 22 -> 11 load step change (0%->50%->100% of rated Po max)
%label = 'CPLPoSweep';  % 6. Response to CPL Po 0->5000 W load ramp change
%label = 'CPLQoSweep';  % 7. Response to CPL Qo 0->2500 W load ramp change
%label = 'rAmpl';       % 8. Response to reference shift (amplitude) step change
%label = 'rAmplSweep';  % 9. Response to reference (amplitude) 0->230rms ramp change
%label = 'rPhase';      % 10. Response to reference shift (phase) step change
%label = 'BusStep';       % 11. Response to supply shift (constant) step change
label = 'BusOsc';        % 12. Response to supply shift (oscilating) step change

%% Save response
%save(['Data/SimOut', label, '.mat'], "out")

%% Load saved response
load(['Data\SimOut', label, '.mat'])

%% Select data ------------------------------------------------------------
idx = (tskip <= out.tout) & (out.tout <= tmax);

out.tout = out.tout(idx);
out.e = out.e(idx);
out.Gc = out.Gc(idx);
out.Ibus = out.Ibus(idx);
out.Io = out.Io(idx);
out.Po = out.Po(idx);
out.Qo = out.Qo(idx);
out.r_inv = out.r_inv(idx);
out.sigma = out.sigma(idx);
out.u = out.u(idx);
out.Vbus = out.Vbus(idx);
out.Vo = out.Vo(idx);
out.VoRect = out.VoRect(idx);
out.VoRMS = out.VoRMS(idx);

clear idx;

%% Calculate statistics ---------------------------------------------------
% Compute THD and Spectrum
toutResampled = out.tout(1):1/fResampling:out.tout(end);
VoResampled = interp1(out.tout, out.Vo, toutResampled);
[f1, harmpow, harmfreq] = thd(VoResampled, fResampling);
fprintf('THD = %.2f dB, Fundamental and harmonics (Hz):\n', f1);
disp(harmfreq);

% Compute ripple
Vripple = out.e;
[Vrup, Vrlo] = envelope(Vripple, 5000, 'peak');
VrEnv = Vrup-Vrlo;
Vripplerms = rms(Vripple);
Vripplepkpk = max(VrEnv);
fprintf([label, ' - Ripple RMS: %.1f V (%.2f %% of Vg_rms), worst pkpk: %.1f V (%.2f %% of Vg_pkpk)\n'], Vripplerms, Vripplerms*sqrt(2)/Vg*100, Vripplepkpk,  Vripplepkpk/(2*Vg)*100)

% Power (bus must be weighted since steps are not constant)
P_in_inst = out.Vbus.*out.Ibus;
t_weight = out.tout(2:end) - out.tout(1:end-1);
P_in = sum(P_in_inst(2:end).*t_weight) / sum(t_weight);
P_out = out.Po(end);

% Efficiency
eff = P_out/P_in*100;
fprintf([label, ' - Power IN: %.1f W, OUT: %.1f W, LOST: %.1f (%.1f %% Eff.)\n'], P_in, P_out, P_in-P_out, eff);

% Instantaneous switching frequency
[ufinst, uftinst] = instfreq(out.u, out.tout); % <--- not exact
fsmax = max(ufinst);
fsmin = min(ufinst);

%% Plot response ----------------------------------------------------------
figure, tiledlayout("vertical");

nexttile, hold on;
plot(out.tout, out.Vo, DisplayName='V_o');
%plot(out.tout, out.VoRect, DisplayName='VoRect');
plot(out.tout, out.Vbus, DisplayName='Vbus');
plot(out.tout, out.r_inv, '--', DisplayName='-r');
%plot(out.tout, out.VoRMS, DisplayName='Vo rms');
ylabel('Voltage (V)');
legend(Location='eastoutside');
grid on;

nexttile, hold on;
%plot(out.tout, out.Ibus, DisplayName='Ibus');
plot(out.tout, out.Io, DisplayName='I_o');
%ylabel('Current (A)');
ylabel('I_o (A)');
%legend();
grid on;

%nexttile, hold on;
%plot(out.tout, out.Vbus.*out.Ibus, DisplayName='Pbus (inst.)');
%plot(out.tout, out.Po, DisplayName='Po (cum.)');
%plot(out.tout, out.Qo, DisplayName='Qo (cum.)');
%ylabel('Power (W)');
%legend(Location='eastoutside');
%grid on;

%nexttile, hold on;
%plot(out.tout, out.u/Vc, DisplayName='u');
%plot(out.tout, out.Gc/Vc, DisplayName='Gc');
%plot(out.tout, out.e/Vc, DisplayName='e');
%plot(out.tout, out.sigma/Vc, DisplayName='σ');
%ylabel('(scaled)');
%legend();

nexttile, hold on;
plot(out.tout, out.u, DisplayName='u');
plot(out.tout, out.Gc, DisplayName='Gc');
plot(out.tout, out.e, DisplayName='e');
plot(out.tout, out.sigma, DisplayName='σ');
ylabel('Voltage (V)');
legend(Location='eastoutside');

nexttile, hold on;
yyaxis left;
plot(out.tout, Vripple, DisplayName='V_e');
%plot(out.tout, Vrup, 'g--', DisplayName='Upper bound');
%plot(out.tout, Vrlo, 'm-.', DisplayName='Lower bound');
xlabel('Time (s)');
%ylabel('Voltage (V)');
ylabel('V_e (V)');
yyaxis right;
plot(uftinst, ufinst/1e3, DisplayName='f_s');
ylabel('f_s (approx, kHz)');
%legend();
grid on;

%sgtitle('Simulated response');

%% Export Plot
ExportPlot(['SimPlot', label], 1, 1.2);

%% Plot Spectrum
figure, hold on;
thd(VoResampled, fResampling, 5);
%xline(fsmin/1e6, 'g--', DisplayName='fs min');
%xline(31600/1e6, 'g--', DisplayName='31.6 kHz');
%xline(fsmax/1e6, 'm-.', DisplayName='fs max');
%xline(91400/1e6, 'm-.', DisplayName='91.4 kHz');
%legend(Location='eastoutside');
legend('off');
axis([0 1 -inf inf]);
xscale log;

%% Export plot
%ExportPlot(['SimSpectrumPlotExtra', label], 1, 0.5);
ExportPlot(['SimSpectrumPlot', label], 1, 0.5);
return

%% Cleanup
clear Vripple P_in P_out eff fResampling tstart tend toutResampled VoResampled f1 harmfreq out tend outDec sys;