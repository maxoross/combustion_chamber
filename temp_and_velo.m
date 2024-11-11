%% Combustion Chamber of an Afterburner (With Velocity Profile)
clear; clc;
% Constants
R = 287.05; % Specific gas constant for air (J/kg·K)
T_air_inlet = 851; % Inlet air temperature (K)
P_air_inlet = 101325; % Inlet air pressure (Pa)
m_dot_air = 0.5; % Mass flow rate of air (kg/s)
m_dot_fuel = 0.01; % Mass flow rate of fuel (kg/s) (Kerosene)

% Chemical reaction: C12H26 + 18.5 O2 -> 12 CO2 + 13 H2O
fuel_formula = 'C12H26';  % Kerosene

% Fuel properties
% We assume Jet-A (C12H26) molecular weight is approximately 170 g/mol
M_fuel = 170 / 1000;  % Molar mass of kerosene (kg/mol)

% Oxygen-to-fuel mass ratio (stoichiometric)
O2_to_fuel = 18.5 * (32 / 170);  % Mass of O2 per mass of fuel (kg O2/kg fuel)

% Air-fuel equivalence ratio (stoichiometric ratio: 1, can vary for excess air)
equivalence_ratio = 1; 

% Calculate fuel-to-air ratio
m_dot_O2 = O2_to_fuel * m_dot_fuel;  % Mass flow rate of O2 (kg/s)
m_dot_N2 = 0.79 * m_dot_air;         % Mass flow rate of N2 (kg/s) (air is 21% O2, 79% N2)

% Reaction products (simplified)
% For stoichiometric combustion:
m_dot_CO2 = 12 * m_dot_fuel;  % Mass flow rate of CO2 (kg/s)
m_dot_H2O = 13 * m_dot_fuel;  % Mass flow rate of H2O (kg/s)

% Exit temperature (simplified energy balance)
% Assuming complete combustion and no heat losses:
Cp_air = 1005; % Specific heat capacity of air (J/kg·K)
Cp_CO2 = 846;  % Specific heat capacity of CO2 (J/kg·K)
Cp_H2O = 1850; % Specific heat capacity of H2O (J/kg·K)

% Energy balance:
LHV_fuel = 43e6; % Lower Heating Value of kerosene (J/kg)

Q_combustion = m_dot_fuel * LHV_fuel; % Total energy released (J/s)

% Assume the temperature rise in the combustion chamber
T_exit = T_air_inlet + (Q_combustion / (m_dot_air * Cp_air));  % Approximate exit temperature

% Velocity Profile Calculation

% Inlet velocity calculation
% Using ideal gas law: p = rho * R * T, so rho = p / (R * T)
rho_air_inlet = P_air_inlet / (R * T_air_inlet); % Air density (kg/m^3)
A_chamber = 0.0033; % Assume a constant cross-sectional area of combustion chamber (m^2)

% Velocity at the inlet (m/s), using continuity equation: v = m_dot / (rho * A)
v_inlet = 180;

% Exit velocity calculation
% Assuming the exit gas behaves ideally, we can estimate the exit density using the exit temperature and pressure
rho_exit = P_air_inlet / (R * T_exit);  % Exit density (kg/m^3)
v_exit = (m_dot_air+m_dot_fuel) / (rho_exit * A_chamber); % Exit velocity (m/s)

% Plotting the velocity profile across the combustion chamber
% Assume the velocity changes linearly along the length of the chamber
% For simplicity, we'll assume a linear temperature gradient, leading to a linear density gradient
length_chamber = 1; % Length of the combustion chamber (m)
n_points = 100; % Number of points to calculate velocity profile

x = linspace(0, length_chamber, n_points); % Distance along the chamber
T_profile = T_air_inlet + (Q_combustion / (m_dot_air * Cp_air)) * (x / length_chamber);  % Temperature profile
rho_profile = P_air_inlet ./ (R * T_profile); % Density profile along the chamber
v_profile = m_dot_air ./ (rho_profile * A_chamber); % Velocity profile along the chamber
% Adjust the velocity profile to include the initial fixed inlet velocity
v_profile = v_inlet + (v_profile - v_inlet) .* (x / length_chamber); % Transition to the new velocity profile

% Plot velocity profile
figure;
subplot(2, 1, 1);
plot(x, v_profile, 'LineWidth', 2);
xlabel('Distance along the combustion chamber (m)');
ylabel('Velocity (m/s)');
title('Velocity Profile in the Afterburner Combustion Chamber');
grid on;

% Display Results
fprintf('Combustion Results:\n');
fprintf('Mass flow rate of air: %.2f kg/s\n', m_dot_air);
fprintf('Mass flow rate of fuel (Kerosene): %.2f kg/s\n', m_dot_fuel);
fprintf('Mass flow rate of O2: %.2f kg/s\n', m_dot_O2);
fprintf('Mass flow rate of N2: %.2f kg/s\n', m_dot_N2);
fprintf('Mass flow rate of CO2: %.2f kg/s\n', m_dot_CO2);
fprintf('Mass flow rate of H2O: %.2f kg/s\n', m_dot_H2O);
fprintf('Estimated exit temperature (T_exit): %.2f K\n', T_exit);
fprintf('Inlet velocity: %.2f m/s\n', v_inlet);
fprintf('Exit velocity: %.2f m/s\n', v_exit);

subplot(2, 1, 2);
plot(x, T_profile, 'LineWidth', 2);
xlabel('Distance along the combustion chamber (m)');
ylabel('Temperature (K)');
title('Temperature Profile in the Afterburner Combustion Chamber');
grid on;
