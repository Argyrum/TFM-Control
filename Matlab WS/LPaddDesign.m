%% Compensator 1: Twin T notch
% Twin T notch is rd = 0, add LP single pole for rd = 1
tauarr = logspace(-2, -8, 200);
karr = zeros(size(tauarr));
barr = zeros(size(tauarr));
ws = fs*2*pi; % Rad/s - Angular switching frequency equivalent

for i = 1:size(tauarr, 2)
    % Implements: Gc1 = tf([1 0 wg^2], [1 4*wg wg^2])*tf(1, [tau 1]);
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

clear coefs1
%% 

% Plot results
figure;
yyaxis left;
plot(tauarr, barr, 'DisplayName', 'b (V)');
ylabel('b (V)');
yyaxis right;
plot(tauarr, karr, 'DisplayName', 'kn');
ylabel('kn')
%title('Relay values vs. LP add-on tau');
xlabel('tau (s)');
%legend;
xscale log;
grid on;

ExportPlot('LPaddDesign');