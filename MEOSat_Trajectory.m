%% Numerical Simulation of a Satellite Orbital Motion in a 2D Plane

% Earth Parameters
R_E = 6378;                       % in km
go  = 9.81e-3;                    % km/s^2

% Initial States of MEOSat Motion
s01 = 3858.213 + R_E;                              % in km
s02 = -5798.143 - R_E;                             % in km
s03 = -0.863;                                      % in km/s                
s04 = -0.542;                                      % in km/s
initial_states = [s01;s02;s03;s04];

%% 
% Time Interval
tspan = 0:0.01:30000;               % in secs

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


% Computation of the Final Position and Velocity of the Satellite
Sat_Positions  = sqrt(X_Positions.^2 + Y_Positions.^2);
Sat_Final_Position = Sat_Positions(length(t));
disp('Position (km) at t = 30,000 seconds:');
disp(Sat_Final_Position);

Sat_Velocities = sqrt(V_Xcomp.^2 + V_Ycomp.^2);
Sat_Final_Velocity = Sat_Velocities(length(t));
disp('Velocity (km/s) at t = 30,000 seconds:');
disp(Sat_Final_Velocity);


%% 
% Plotting the Orbital Plane in a 2D Space
fig1 = figure();

plot(X_Positions,Y_Positions,'b-','LineWidth',2);
xlabel('x (km)');
ylabel('y (km)');
title('2D Orbit Propagation Using ode45');
grid on

% Plotting the Time Evolution of the Satellite Position
fig2 = figure();
plot(t, Sat_Positions, "LineWidth",2, "Color","k");
ylim([0 20000]);
xlabel('time (s)');
ylabel('Satellite Position (km)');
title(' Time Evolution of Satellite Position');
grid on