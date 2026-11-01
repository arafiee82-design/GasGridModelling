%% Steady-state hydrogen flow in a gas grid
% Network: 200 km main pipe (4 x 50 km) + 4 connections (20 km) to 4 GTs
clear; clc;

%% ---------- Gas and environment properties ----------
p.M    = 2.016e-3;     % Molar mass of H2 [kg/mol]
p.R    = 8.314;        % Universal gas constant [J/(mol K)]
p.T    = 288.15;       % Gas temperature [K] (assumed)
p.Z    = 1.0;          % Compressibility factor (ideal gas assumption)
p.mu   = 9.0e-6;       % Dynamic viscosity of H2 [Pa s]
p.eps  = 0.045e-3;     % Pipe roughness (steel) [m]
LHV    = 120e6;        % Lower heating value of H2 [J/kg]

%% ---------- Network geometry and operating limits ----------
p.P_in = 70e5;         % Inlet pressure [Pa] (70 bar)
p.P_min = 20e5;        % Minimum allowable GT inlet pressure [Pa] (assumed)
p.L1   = 50e3;         % Main segment length [m]
p.L2   = 20e3;         % Connection length [m]
p.v_max = 40;          % Maximum allowed gas velocity [m/s]
nGT    = 4;

%% ---------- Single case: given diameters ----------
GT_MW = 50;                          % GT size [MW]
D1    = 0.600;                       % Main pipe diameter [m]
D2    = 0.150;                       % Connection diameter [m]

mGT   = (GT_MW*1e6/LHV) * ones(1,nGT);   % Mass flow per GT [kg/s]
out   = simulateGrid(mGT, D1, D2, p);

fprintf('--- Single case: GT = %d MW, D1 = %d mm, D2 = %d mm ---\n', GT_MW, D1*1e3, D2*1e3);
fprintf('Main segment pressures [bar]: '); fprintf('%.3f  ', out.P_main/1e5); fprintf('\n');
fprintf('GT inlet pressures    [bar]: '); fprintf('%.3f  ', out.P_GT/1e5);   fprintf('\n');
fprintf('Max velocity          [m/s]: %.2f\n\n', out.Vmax);

%% ---------- Sweep: smallest diameters that satisfy constraints ----------
GT_list = 50:50:500;
D1_list = [600 750 900]*1e-3;
D2_list = [150 200 250 300 400 450 500 600 650 750 900]*1e-3;

fprintf('%-8s %-10s %-10s %-12s %-12s\n', 'GT[MW]', 'D1[mm]', 'D2[mm]', 'Vmax[m/s]', 'P_GT,min[bar]');
for g = GT_list
    mGT = (g*1e6/LHV) * ones(1,nGT);
    found = false;
    for d1 = D1_list
        for d2 = D2_list
            o = simulateGrid(mGT, d1, d2, p);
            if o.Vmax <= p.v_max && min(o.P_GT) >= p.P_min
                fprintf('%-8d %-10d %-10d %-12.2f %-12.2f\n', g, d1*1e3, d2*1e3, o.Vmax, min(o.P_GT)/1e5);
                found = true;
                break;
            end
        end
        if found, break; end
    end
    if ~found
        fprintf('%-8d  no feasible combination in the given diameter list\n', g);
    end
end

%% ========================================================================
function out = simulateGrid(mGT, D1, D2, p)
% mGT : 1 x nGT vector of GT mass flows [kg/s]
% D1  : main pipe diameter [m]
% D2  : connection diameter [m]
nGT = numel(mGT);

% Main segment k (from node k-1 to node k) carries demand of GT k..nGT
mainFlow = fliplr(cumsum(fliplr(mGT)));

P_main = zeros(1, nGT+1);
P_main(1) = p.P_in;
P_GT = zeros(1, nGT);
Vmax = 0;

for k = 1:nGT
    % Main transmission segment k (50 km)
    [P_main(k+1), vM] = segmentOutlet(P_main(k), mainFlow(k), D1, p.L1, p);

    % Connection from node k to GT k (20 km)
    [P_GT(k), vC] = segmentOutlet(P_main(k+1), mGT(k), D2, p.L2, p);

    Vmax = max([Vmax, vM, vC]);
end

out.P_main = P_main;
out.P_GT   = P_GT;
out.Vmax   = Vmax;
end

%% ========================================================================
function [Pout, vout] = segmentOutlet(Pin, mdot, D, L, p)
% Isothermal pressure drop for a single pipe segment (steady state)
%   Pin^2 - Pout^2 = 16 f mdot^2 Z R T L / (pi^2 D^5 M)

A  = pi*D^2/4;                       % Cross-section area [m^2]
Re = 4*mdot/(pi*D*p.mu);             % Reynolds number (independent of pressure)

% Swamee-Jain explicit friction factor (approximates Colebrook)
f  = 0.25 / (log10(p.eps/(3.7*D) + 5.74/Re^0.9))^2;

Pout = sqrt(Pin^2 - 16*f*mdot^2*p.Z*p.R*p.T*L/(pi^2*D^5*p.M));

% Velocity at outlet (lowest pressure, highest velocity in the segment)
vout = mdot*p.Z*p.R*p.T/(Pout*p.M*A);
end