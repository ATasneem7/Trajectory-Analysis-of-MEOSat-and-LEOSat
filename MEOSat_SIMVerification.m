%% Verification of the Simulation Results of the Satellite Orbital Motion

% Earth Parameters
R_E = 6378;                                % in km
go  = 9.81e-3;                             % km/s^2

% Initial Conditions
s01 = 3858.213+R_E;                        % in km
s02 = -5798.143-R_E;                       % in km
s03 = -0.863;                              % in km/s                
s04 = -0.542;                              % in km/s
initial_states = [s01;s02;s03;s04];

%% 
% Time Interval
tspan = 0:0.01:30000;               % in secs% 

% Error Tolerance
tolerance = 1e-9;
options = odeset("RelTol",tolerance, "AbsTol", tolerance);

% Implementation of ODE45 Numerical Solver
[t, S] = ode45(@MEOSat, tspan, initial_states, options, go, R_E);

%% 
% Extracting the Position and Velocity Data from the State Vector
X_Positions = S(:, 1);
Y_Positions = S(:, 2);
V_Xcomp = S(:, 3);
V_Ycomp = S(:,4);


% Computing Total Specific Energy at each time step
Energy = zeros(length(t), 1);
for i = 1:length(t)
    r = norm(S(i, 1:2));
    v = norm(S(i, 3:4));
    Energy(i) = 0.5 * v^2 - go*R_E^2 / r;       % in km^2/s^2    
end

%% 
% Calculate the initial specific energy
Energy_initial = Energy(1);

% Calculate the relative error in specific energy
Energy_error = abs((Energy - Energy_initial) / Energy_initial);
% Visualization of the Orbital Specific Energy and its Relative Error

%% 
% Plot energy over time
fig1 = figure();
plot(t, Energy, 'b', 'LineWidth', 0.00001);
grid on;
xlim([2500 30000])

xlabel('Time (s)');
ylabel('Specific Energy (km^2/s^2)');
title('Specific Orbital Energy Over Time')

% Plot Relative Error in Specific Energy
fig2 = figure();
plot(t, Energy_error, 'r', 'LineWidth', 0.001);
grid on;
xlim([200 30000]);
ylim([0 7e-9])
xlabel('Time (s)');
ylabel('Relative Energy Error');
title('Relative Error in Specific Orbital Energy');