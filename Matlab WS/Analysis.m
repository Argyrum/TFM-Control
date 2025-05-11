% 1. Response to R load
% 2. Response to RL load
% 3. Response to L + rectified C load
% 4. Response to reference shift (amplitude)
% 5. Response to reference shift (phase)
% 5. Response to supply shift (constant)
% 5. Response to supply shift (oscilating)

%% Clean environment
close all;
clear variables;
clc;
format shorte;

%% Set up workspace path
if ~contains(path,"Data\")
    addpath("Simulink models\", "Data\");
end

%% Input parameters and component values
load Data\SimulationWorkspace.mat

%% Configure simulation
tstep = 1e-7; % Simulation time step
tmax = (1/fg)*4; % Simulation end time

outDec = 10; % Decimation factor for simulation out values

%% Load simulation results (skip execution)

load Data\OutDemoSim.mat

%sys = "CLRFSPFCRealizedSimscape";
%load_system(sys);
%close_system;

%% Response analysis
figure, tiledlayout("vertical");

nexttile, hold on;
plot(out.tout, out.Vo, 'DisplayName', 'Vo');
plot(out.tout, out.VoRect, 'DisplayName', 'VoRect');
plot(out.tout, out.E, 'DisplayName', 'E');
plot(out.tout, out.r_inv, 'DisplayName', '-r');
legend();

nexttile, hold on;
plot(out.tout, out.Io, 'DisplayName', 'Io');
legend();

nexttile, hold on
plot(out.tout, out.Po, 'DisplayName', 'Po');
plot(out.tout, out.Qo, 'DisplayName', 'Qo');
legend();

nexttile, hold on;
plot(out.tout, out.u/Vc, 'DisplayName', 'u (scaled)');
plot(out.tout, out.Gc/Vc, 'DisplayName', 'Gc (scaled)');
plot(out.tout, out.e, 'DisplayName', 'e');
legend();

nexttile, hold on;
plot(out.tout, out.THD, 'DisplayName', 'THD');
yscale log;
legend();

sgtitle('Simulated response');