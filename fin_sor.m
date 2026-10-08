%% Steady 1-D fin: SOR with temperature-dependent conductivity
% Source: MA6801 Heat Conduction lecture notes, PDF pages 35-37.
% Run this script. x, T, iterations and residualK remain in the workspace.
% Each iteration is ONE left-to-right sweep with frozen face conductivity.
% Insulated tip: zero end-face flux, but lateral convection over dx/2.
clear; clc;

%% Inputs (SI units; temperatures in kelvin)
Tb = 300; Ta = 20; hc = 100;
D = 0.005; L = 0.10;
N = 51;                         % Includes base and tip; editable
omega = 1.5;                  % 1 = GS; 1 < omega < 2 = SOR
tol = 1e-3;                    % K: maximum temperature change per sweep
maxIter = 200000;
assert(N >= 3 && N == floor(N), 'N must be an integer >= 3.');
assert(omega > 0 && omega < 2, 'Require 0 < omega < 2.');
A = pi*D^2/4; P = pi*D;
dx = L/(N-1); x = linspace(0,L,N)';
kfun = @(temp) 10.3 + (112-10.3)/(300-20)*(temp-20);
C = hc*P*dx^2/A;                % Units W/(m K), same as conductivity
T = ((Tb+Ta)/2)*ones(N,1);       % Initial guess: 160 K at unknown nodes
T(1) = Tb;                     % Prescribed boundary overrides initial guess
converged = false;

%% Nonlinear iteration: update conductivity BEFORE each complete sweep
for iterations = 1:maxIter
    Told = T;
    kface = kfun((Told(1:end-1)+Told(2:end))/2);
    if any(~isfinite(kface)) || any(kface <= 0)
        error('Invalid conductivity. Reduce omega or check input data.');
    end
    for i = 2:N-1
        kw = kface(i-1); ke = kface(i);
        Tgs = (kw*T(i-1) + ke*Told(i+1) + C*Ta)/(kw+ke+C);
        T(i) = Told(i) + omega*(Tgs-Told(i));
    end
    % Tip half-control-volume: kw*(T(N-1)-T(N)) + C/2*(Ta-T(N)) = 0
    kw = kface(end);
    Tgs = (kw*T(N-1) + (C/2)*Ta)/(kw+C/2);
    T(N) = Told(N) + omega*(Tgs-Told(N));
    T(1) = Tb;
    if any(~isfinite(T))
        error('Iteration diverged. Reduce omega.');
    end
    maxChange = max(abs(T-Told));

    % Independent balance check using UPDATED conductivity.
    % Divide equation imbalance by diagonal coefficient to express in K.
    knew = kfun((T(1:end-1)+T(2:end))/2);
    kw = knew(1:end-1); ke = knew(2:end);
    r = kw.*(T(1:end-2)-T(2:end-1)) ...
        + ke.*(T(3:end)-T(2:end-1)) + C*(Ta-T(2:end-1));
    rtip = knew(end)*(T(N-1)-T(N)) + C/2*(Ta-T(N));
    residualK = max([abs(r)./(kw+ke+C); abs(rtip)/(knew(end)+C/2)]);
    if maxChange <= tol && residualK <= tol
        converged = true;
        break
    end
end
if ~converged
    error('No convergence after %d sweeps (change %.3g K, residual %.3g K).', ...
        maxIter,maxChange,residualK);
end

%% Report and plot every node
fprintf('SOR: N = %d, omega = %.2f\n',N,omega);
fprintf('Iterations (complete sweeps): %d\n',iterations);
fprintf('Maximum final change: %.6g K\n',maxChange);
fprintf('Maximum diagonal-normalized residual: %.6g K\n',residualK);
fprintf('Tip temperature: %.6f K\n',T(end));
figure('Color','w','Name','SOR fin temperature');
plot(x,T,'-o','LineWidth',1.4,'MarkerSize',4);
grid on; xlabel('Position x (m)'); ylabel('Temperature T (K)');
title(sprintf('SOR: %d nodes, %d iterations',N,iterations));
xlim([0 L]);
% Optional: writetable(table(x,T),'fin_temperatures.csv');
% The 1e-3 K iteration tolerance is not a bound on discretization/global error.
% Refine N and tighten tol to check numerical accuracy.
