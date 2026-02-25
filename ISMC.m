
clear all; close all; clc;

% System Parameters 
M = 0.23;       % Cart mass (kg)
m = 2.4;        % Pendulum mass (kg)
l = 0.45;       % Length of pendulum (m)
I = 0.38;       % Inertia (kg*m^2)
g = 9.8;        % Gravity (m/s^2)
Bx = 0.005;     % Friction of cart (N)

%Control Parameters 
lambda = 10;    % Tuning parameter Lambda
k = 3;          % Gain k 
zeta = 1;       % Gain zeta 
epsilon = 1;  % Boundary layer thickness for chattering reduction

%Simulation Setup
tspan = 0:0.001:50;             
x0 = [0; 0; 0.1; 0];            % Initial states 
xd = 1;                         % Desired cart position 

%% Run Simulation
[t, states] = ode45(@(t, y) pendulum_dynamics_ismc(t, y, M, m, l, I, g, Bx, xd, ...
    lambda, k, zeta, epsilon), tspan, x0);

%% Plotting Results
figure('Color', 'w');
subplot(2,1,1);
plot(t, states(:,1), 'r', 'LineWidth', 2); hold on;
yline(xd, '--k', 'Desired');
title('Step Response of Cart Position (ISMC)');
ylabel('Position (m)'); grid on;

subplot(2,1,2);
plot(t, states(:,3), 'b', 'LineWidth', 2);
title('Step Response of Pendulum Angle (ISMC)');
ylabel('Angle (rad)'); xlabel('Time (sec)'); grid on;

%% Dynamics & ISMC Control Function
function dy = pendulum_dynamics_ismc(t, y, M, m, l, I, g, Bx, xd, ...
    lambda, k, zeta, eps)

    % States
    x = y(1);    % Cart position
    v = y(2);    % Cart velocity
    theta = y(3); % Pendulum angle
    omega = y(4); % Angular velocity

    % Nonlinear Dynamics Matrices (Newton approach)
    % Mass matrix M(q) and Vector V(q, q_dot)
    detM = (M + m)*(I + m*l^2) - (m*l*cos(theta))^2;
    
    % Tracking Error
    e = x - xd;
    
    
    
    s = omega + lambda^3*e + 3*lambda^2*v + 3*lambda*theta;

    
    
    if abs(s) > eps
        sat_s = sign(s);
    else
        sat_s = s/eps;
    end
    u_dis = -k * sat_s - zeta * s;

    % Control Force F
    F = - (lambda^3*v + 3*lambda^2*omega) + u_dis; 

    % Solving for accelerations
    dy = zeros(4,1);
    dy(1) = v;
    dy(3) = omega;
    
    % M * [v_dot; omega_dot] = [F + m*l*omega^2*sin - Bx*v; -m*g*l*sin]
    A = [ (M+m), m*l*cos(theta); m*l*cos(theta), (I + m*l^2) ];
    B = [ F + m*l*omega^2*sin(theta) - Bx*v; -m*g*l*sin(theta) ];
    accels = A \ B;
    
    dy(2) = accels(1); % x_ddot
    dy(4) = accels(2); % theta_ddot
end