clc; clear; close all;

% Define system parameters 
m = 2.4;       % Pendulum mass [kg]
M = 0.23;      % Cart mass [kg]
I = 0.38;      % Pendulum inertia [kg*m^2]
g = 9.8;       % Gravity [m/s^2]
Bx = 0.005;    % Cart friction [N*s/m]
l = 0.45;      % Pendulum length [m]
L = 0.91;      % Cart length [m]

% Denominator factor for transfer functions 
q = (M + m)*(I + M*l^2) - (m*l^2);  

%% Linearized state-space matrices

A = [0       1       0       0;
     0   -0.0194   0.2188    0;
     0       0       0       1;
     0   -0.0128   6.5848    0];

B = [0;
     0.3887;
     0;
     0.2552];

C = [1 0 0 0;   
     0 0 1 0];  

D = [0;
     0];

sys = ss(A,B,C,D);

open_loop_poles = pole(sys);
disp('Open-loop poles (Eq. 8):');
disp(open_loop_poles);

%% LQR design

Q = diag([50 0 100 0]);  
R = 0.5;                  

% Compute LQR gain K 
K = lqr(A,B,Q,R);

Acl = A - B*K;
sys_cl = ss(Acl,B,C,D);

% Closed-loop poles
cl_poles = pole(sys_cl);
disp('Closed-loop poles with LQR:');
disp(cl_poles);

%% Simulation 

t = 0:0.01:10;          
x0 = [0.5; 0; 0.1; 0];  

r = [ones(size(t)); zeros(size(t))];  

u_func = @(t,x) -K*(x - [r(1,round(t/0.01)+1);0;0;0]);


[t_sim, x_sim] = ode45(@(t,x) A*x + B*u_func(t,x), t, x0);

%% Plot results

figure;
plot(t_sim, x_sim(:,1), 'LineWidth', 1.5); hold on;
plot(t, r(1,:), '--', 'LineWidth', 1.2);
xlabel('Time [s]'); ylabel('Cart Position [m]');
title('Cart Position - LQR Control'); grid on;
legend('x(t)','Reference');

figure;
plot(t_sim, x_sim(:,3), 'LineWidth', 1.5); hold on;
plot(t, r(2,:), '--', 'LineWidth', 1.2);
xlabel('Time [s]'); ylabel('Pendulum Angle [rad]');
title('Pendulum Angle - LQR Control'); grid on;
legend('\theta(t)','Reference');

u_sim = zeros(length(t_sim),1);
for i = 1:length(t_sim)
    u_sim(i) = -K*x_sim(i,:)';
end

figure;
plot(t_sim, u_sim, 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('Control Input [N]');
title('LQR Control Input'); grid on;