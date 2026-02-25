clc; clear; close all;

%% System Parameters
m = 2.4;       % Pendulum mass [kg]
M = 0.23;      % Cart mass [kg]
I = 0.38;      % Pendulum inertia [kg*m^2]
g = 9.8;       % Gravity [m/s^2]
Bx = 0.005;    % Cart friction [N*s/m]
l = 0.45;      % Pendulum length [m]


x0 = [0.5; 0; 0.1; 0];  
tspan = 0:0.01:10;  

r = ones(size(tspan));  % Unit step reference

%% SMC

% Sliding Mode Control Gains
k3 = 70;   % derivative gain
k2 = 37;
k = 20;    % sliding surface gain

% SMC Control Function
u_func = @(t,x) smc_control(x, r(round(t/0.01)+1), k, k2, k3, m, M, I, g, Bx, l);

%% ODE45
odefun = @(t,x) nonlinear_IP(x, u_func(t,x), m, M, I, g, Bx, l);

[t_sim, x_sim] = ode45(odefun, tspan, x0);

% Compute control input history
u_sim = zeros(length(t_sim),1);
for i = 1:length(t_sim)
    u_sim(i) = u_func(t_sim(i), x_sim(i,:)');
end

%% Plot results
figure;
plot(t_sim, x_sim(:,1), 'LineWidth',1.5); hold on;
plot(tspan, r, '--', 'LineWidth',1.2);
xlabel('Time [s]'); ylabel('Cart Position [m]');
title('SMC: Cart Position'); grid on;
legend('x(t)','Reference');

figure;
plot(t_sim, x_sim(:,3), 'LineWidth',1.5);
xlabel('Time [s]'); ylabel('Pendulum Angle [rad]');
title('SMC: Pendulum Angle'); grid on;

figure;
plot(t_sim, u_sim, 'LineWidth',1.5);
xlabel('Time [s]'); ylabel('Control Input [N]');
title('SMC: Control Input'); grid on;


%% Nonlinear inverted pendulum dynamics
function dx = nonlinear_IP(x,u,m,M,I,g,Bx,l)
    s = sin(x(3)); c = cos(x(3));
    D = I*(M+m) + M*m*l^2 + m^2*l^2*(1-c^2);  % Denominator

    dx = zeros(4,1);
    dx(1) = x(2);
    dx(2) = (I + m*l^2)*(u - Bx*x(2) + m*l*s*x(4)^2) + m^2*l^2*g*s*c;
    dx(2) = dx(2)/D;
    dx(3) = x(4);
    dx(4) = -(m*l*c)*(u - Bx*x(2) + m*l*s*x(4)^2) + (M+m)*m*g*l*s;
    dx(4) = dx(4)/D;
end

%% Sliding Mode Control function
function u = smc_control(x, xd, k, k2, k3, m, M, I, g, Bx, l)
    % Tracking errors (cart position)
    e = x(1) - xd;  
    ed = x(2);       % cart velocity error
    
    % Sliding surface
    s = ed + k3*e;   % simple first-order sliding for single input
    
    % Dynamics terms (for equivalent control)
    s_theta = sin(x(3)); c_theta = cos(x(3));
    D = I*(M+m) + M*m*l^2 + m^2*l^2*(1-c_theta^2);
    
    % f(x) term for cart acceleration
    f = (I + m*l^2)*( -Bx*x(2) + m*l*s_theta*x(4)^2 ) + m^2*l^2*g*s_theta*c_theta;
    f = f / D;
    
    g_term = (I + m*l^2)/D;  % input gain
    
    % Equivalent control
    u_eq = (-f - k3*ed)/g_term;
    
    % Discontinuous control
    u_dis = -k*sign(s);
    
    % Total control
    u = u_eq + u_dis;
end