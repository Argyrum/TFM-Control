%% OL Ripple Test ---------------------------------------------------------
mdl = "OLRippleTest";
load_system(mdl);
simIn = Simulink.SimulationInput(mdl);
load("Data\OLRippleTestOP.mat");
simIn = setInitialState(simIn, OLRippleTestOP);
simIn = setModelParameter(simIn,"StopTime","1.0001");
out = sim(simIn);
close_system(mdl, 0);

%% Plot results
figure;
tiledlayout('vertical');
idx = out.u.Time > 1;

nexttile;
hold on;
plot(out.u.Time(idx)-1, out.u.Data(idx, 1));
plot(out.u.Time(idx)-1, out.u.Data(idx, 2));
ylabel('u');
grid on
axis padded
hold off;

nexttile;
hold on;
plot(out.V.Time(idx)-1, out.V.Data(idx, 1));
plot(out.V.Time(idx)-1, out.V.Data(idx, 2));
ylabel('Vo (V)');
grid on
axis padded
hold off;

nexttile;
hold on;
plot(out.I.Time(idx)-1, out.I.Data(idx, 1));
plot(out.I.Time(idx)-1, out.I.Data(idx, 2));
ylabel('Io (A)');
xlabel('Time (s), from t=1')
grid on
axis padded
hold off;

%% Export plot
ExportPlot('OLRippleTest');

%% Cleanup
clearvars mdl out simIn OLRippleTestOP idx;

%% LP tau Design ----------------------------------------------------------
% Compensator 1: Twin T notch
tauarr = logspace(-2, -8, 200);
karr = zeros(size(tauarr));
barr = zeros(size(tauarr));
ws = fs*2*pi; % Rad/s - Angular switching frequency equivalent

for i = 1:size(tauarr, 2)
    coefs1 = tauarr(i); % [tau]
    Gc1 = tf([1 0 wg^2], [1 4*wg wg^2])*tf(1, [coefs1(1) 1]);
    Ga1 = G + Gc1;    
    [isStable, ~, ~, ~] = CheckTF(Gc1, 'Gc1', false, false);

    if ~isStable
        barr(i) = 0;
        karr(i) = 0;
        continue
    end

    [isStable, ~, ~, isASPR] = CheckTF(Ga1, 'Ga1', false, false);

    if ~isStable | ~isASPR
        barr(i) = 0;
        karr(i) = 0;
        continue
    end

    Gass = ss(Ga1);
    [barr(i), isStable, karr(i), ~] = Compute_b(Gass.A, Gass.B, Gass.C, c, ws, false);
    if ~isStable
        barr(i) = 0;
        karr(i) = 0;
        continue
    end
end

%% Plot results
figure;
yyaxis left;
plot(tauarr, barr, 'DisplayName', 'b (V)');
ylabel('b (V)');
yyaxis right;
plot(tauarr, karr, 'DisplayName', 'kn');
ylabel('kn')
%title('Relay values vs. LP add-on tau');
xlabel('τ (s)');
xscale log;
grid on;

%% Export plot
ExportPlot('LPtauDesign', 1, 0.5);

%% Cleanup
clear coefs1 ws tauarr barr karr

%% Gc Rejection Comparison ------------------------------------------------
t = 0:0.00001:0.25;
u = Vg*sin(wg*t);

y1 = lsim(Gc1, u, t);
y2 = lsim(Gc2, u, t);
y3 = lsim(Gc3, u, t);
y4 = lsim(Gc4, u, t);

%% Plot results
figure, hold on;
plot(t, u, 'DisplayName', 'Mains');
plot(t, y1, 'DisplayName', 'Gc1');
plot(t, y2, 'DisplayName', 'Gc2');
plot(t, y3, 'DisplayName', 'Gc3');
plot(t, y4, 'DisplayName', 'Gc4');
hold off;
legend;
ylabel('Amplitude (V)');
xlabel('Time (s)');
grid on;

%% Export plot
ExportPlot('GcRejectionComparison');

%% Cleanup
clear t u y1 y2 y3 y4;

%% Gc offset Bode plot ----------------------------------------------------
% Compensator 1: Twin T notch
figure, hold on;
bp = bodeplot(Gc1, Gc, {wg/fg, 2*wg});
hold off;
bp.FrequencyScale = "linear";
bp.FrequencyUnit = "Hz";
bp.LegendLocation = "southeast";
bp.Title.String = "";
legend(["Designed", "Realization"]);
bp.LegendVisible = false;

%% Export
ExportPlot('GcOffsetBodePlot');

%% Cleanup
clear bp;

%% G LPRS characteristics -------------------------------------------------
CheckLPRS(G, 'G', warr, c, fs, true, true);
ExportPlot('LPRS_G');
CheckLPRS(Ga, 'Ga', warr, c, fs, true, true);
ExportPlot('LPRS_Ga');

%% SSR Test ---------------------------------------------------------------
mdl = "SSRTest";
load_system(mdl);
simIn = Simulink.SimulationInput(mdl);
simIn = setModelParameter(simIn,"StopTime","6.67e-6");
out = sim(simIn);
close_system(mdl, 0);

%% Plot results
figure;
tiledlayout('horizontal');

nexttile;
hold on;
plot(out.Vin.Time, out.Vin.Data, DisplayName='Vin');
plot(out.Vout.Time, out.Vout.Data, DisplayName='Vout');
yline(bS, DisplayName='Design +b');
yline(-bS, DisplayName='Design -b');
yline(Vc*Rst1/Rst2, '--', DisplayName='Realization +b');
yline(-Vc*Rst1/Rst2, '--', DisplayName='Realization -b');
hold off;
lgd = legend(Location='north', NumColumns=3);
lgd.Layout.Tile = 'north';
xlabel('Time (s)');
ylabel('(V)');
grid on;
axis([0 6.67e-6 -inf inf]);

nexttile;
hold on;
plot(out.Vin.Time, out.Vin.Data, DisplayName='Vin');
plot(out.Vout.Time, out.Vout.Data, DisplayName='Vout');
yline(bS, DisplayName='Design +b');
yline(-bS, DisplayName='Design -b');
yline(Vc*Rst1/Rst2, '--', DisplayName='Realization +b');
yline(-Vc*Rst1/Rst2, '--', DisplayName='Realization -b');
axis([1.5e-6 1.86e-6 -1 1])
hold off;
grid on;

nexttile;
hold on;
plot(out.Vin.Time, out.Vin.Data, DisplayName='Vin');
plot(out.Vout.Time, out.Vout.Data, DisplayName='Vout');
yline(bS, DisplayName='Design +b');
yline(-bS, DisplayName='Design -b');
yline(Vc*Rst1/Rst2, '--', DisplayName='Realization +b');
yline(-Vc*Rst1/Rst2, '--', DisplayName='Realization -b');
axis([4.86e-6 5.18e-6 -1 1])
hold off;
grid on;

%% Export plot
ExportPlot('SSRTest', 1.5, 1);

%% Cleanup
clearvars mdl out;