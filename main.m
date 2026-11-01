clear; clc; close all;

%% 1. Problem definition

t0 = 0;
tf = 1;
y0 = 1;
N = 10;  % Number of interior Legendre points

% Reference solution
yExact = @(t) 1 ./ (1 + 2*t);

%% 2. Solve the ODE using ode23
tPlot = linspace(t0, tf, 201);

optsODE = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);

[tODE, yODE] = ode23( ...
    @(t,y) -2*y.^2, tPlot, y0, optsODE);

%% 3. Build Legendre-based collocation nodes
z = legendreRoots(N);
gaussNodes = sort((z + 1)/2);

% Include the initial and final endpoints
x = [0; gaussNodes; 1];
M = numel(x);

%% 4. Construct the differentiation matrix
% Barycentric weights
w = ones(M,1);

for i = 1:M
    for j = 1:M
        if i ~= j
            w(i) = w(i) / (x(i) - x(j));
        end
    end
end

D = zeros(M,M);

for i = 1:M
    for j = 1:M
        if i ~= j
            D(i,j) = w(j) / (w(i)*(x(i)-x(j)));
        end
    end

    D(i,i) = -sum(D(i,[1:i-1, i+1:M]));
end

%% 5. Solve the nonlinear algebraic system with fsolve
% The initial condition y(0) = 1 is already known.
% Unknowns are y at all remaining collocation nodes.
residual = @(zUnknown) ...
    D(2:end,:) * [y0; zUnknown] ...
    + 2 * zUnknown.^2;

z0 = yExact(x(2:end));

optsFS = optimoptions('fsolve', ...
    'Display','off', ...
    'FunctionTolerance',1e-12, ...
    'StepTolerance',1e-12);

[zSol, fval, exitflag] = fsolve( ...
    residual, z0, optsFS);

yNodes = [y0; zSol];

%% 6. Evaluate the collocation polynomial
yColl = barycentricEval(x, yNodes, w, tPlot);

yRef = yExact(tPlot);

%% 7. Compare errors
% Evaluate both methods on the same time grid
yODEcommon = interp1(tODE, yODE, tPlot, 'pchip');
yExactCommon = yExact(tPlot);

errODE = abs(yODEcommon - yExactCommon);
errColl = abs(yColl - yExactCommon);

figure;

semilogy(tPlot, max(errODE, eps), ...
    'b-', 'LineWidth', 1.5);
hold on;

semilogy(tPlot, max(errColl, eps), ...
    'r--', 'LineWidth', 1.5);

xlabel('t');
ylabel('Absolute error');
title('Numerical Error Comparison');
legend('ode23', 'Collocation', 'Location', 'best');
grid on;

%% 8. Plot solutions
figure;

plot(tPlot, yRef, 'k-', 'LineWidth', 2);
hold on;
plot(tODE, yODE, 'b--', 'LineWidth', 1.5);
plot(tPlot, yColl, 'r:', 'LineWidth', 2);
plot(x, yNodes, 'ro', 'MarkerFaceColor','r');

xlabel('t');
ylabel('y(t)');
title('Orthogonal Collocation vs ode23');
legend('Exact solution', 'ode23', ...
    'Collocation polynomial', 'Collocation nodes', ...
    'Location','best');
grid on;

fprintf('ode23 error = %.3e\n', max(errODE));
fprintf('Collocation error = %.3e\n', max(errColl));
fprintf('Nonlinear residual = %.3e\n', norm(fval, inf));
fprintf('fsolve exitflag = %d\n', exitflag);