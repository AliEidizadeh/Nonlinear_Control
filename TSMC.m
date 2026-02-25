
clear all; close all; clc;

% System Parameters 
M = 0.23;       % Cart mass (kg)
m = 2.4;        % Pendulum mass (kg)
l = 0.45;       % Length of pendulum (m)
I = 0.38;       % Inertia (kg*m^2)
g = 9.8;        % Gravity (m/s^2)
Bx = 0.005;     % Friction of cart (N)

% Control Parameters 
lambda = 20;    
k = 10;         
zeta = 20;      
alpha = 30;     
beta = 40;      
p_tsmc = 5;     
q_tsmc = 3;     

% Simulation Setup
tspan = 0:0.01:50;              
x0 = [0; 0; 0.1; 0];           
xd = 1;                         

%% Run Simulation
[t, states] = ode45(@(t, y) pendulum_dynamics(t, y, M, m, l, I, g, Bx, xd, ...
    lambda, k, zeta, alpha, beta, p_tsmc, q_tsmc), tspan, x0);

%% Plotting Results
figure;
subplot(2,1,1);
plot(t, states(:,1), 'm', 'LineWidth', 1.5); hold on;
yline(xd, '--k');
title('Step Response of Cart Position (TSMC)');
ylabel('Position (m)'); grid on; 

subplot(2,1,2);
plot(t, states(:,3), 'm', 'LineWidth', 1.5);
title('Step Response of Pendulum Angle (TSMC)');
ylabel('Angle (rad)'); xlabel('Time (sec)'); grid on;

%% Dynamics & Control Function
function dy = pendulum_dynamics(t, y, M, m, l, I, g, Bx, xd, ...
    lambda, k, zeta, alpha, beta, p, q)

    % States
    x1 = y(1); % Cart position
    x2 = y(2); % Cart velocity
    x3 = y(3); % Pendulum angle
    x4 = y(4); % Angular velocity

    % Nonlinear Functions
    den = M + m * sin(x3)^2;
    o2 = ((M+m)*g*sin(x3) - m*l*x4^2*cos(x3)) / (l * den);
    i2 = -cos(x3) / (l * den);

    % Tracking Errors
    e = x1 - xd;
    
    % TSMC Sliding Surface S
    s = x4 + lambda^3*e + 3*lambda^2*x2 + beta*sign(e)*abs(e)^(q/p); 

    % Equivalent Control Part (ueq)
    u_eq = -(o2 + lambda^3*x2 + beta*(q/p)*abs(e)^(q/p-1)*x2) / i2;

    % Discontinuous Control Part (udis)
    u_dis = -k * sign(s) - zeta * s;

    % Total Control Input
    F = u_eq + u_dis;

    % Nonlinear State Equations
    dy = zeros(4,1);
    dy(1) = x2;
    dy(3) = x4;
    
    mat_A = [ (M+m), m*l*cos(x3); m*l*cos(x3), (I + m*l^2) ];
    vec_B = [ F + m*l*x4^2*sin(x3) - Bx*x2; -m*g*l*sin(x3) ];
    accels = mat_A \ vec_B;
    
    dy(2) = accels(1); % x_ddot
    dy(4) = accels(2); % theta_ddot
end