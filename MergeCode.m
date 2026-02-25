clc; clear; close all;

%  System Parameters

m = 2.4;
M = 0.23;
I = 0.38;
g = 9.8;
Bx = 0.005;
l = 0.45;


%  Simulation Setup
dt    = 0.01;
tspan = 0:dt:10;
x0    = [0.5; 0; 0.1; 0];
xd    = 1;                     % reference cart position
r     = xd * ones(size(tspan));


%%  LQR Controller

A = [0 1 0 0;
     0 -0.0194 0.2188 0;
     0 0 0 1;
     0 -0.0128 6.5848 0];

B = [0; 0.3887; 0; 0.2552];

Q = diag([50 0 100 0]);
R = 0.5;
K = lqr(A,B,Q,R);

u_lqr = @(t,x) -K*(x - [xd;0;0;0]);
[t_lqr,x_lqr] = ode45(@(t,x) A*x + B*u_lqr(t,x), tspan, x0);

uLQR = zeros(length(t_lqr),1);
for i=1:length(t_lqr)
    uLQR(i) = u_lqr(t_lqr(i),x_lqr(i,:)');
end

%%  SMC Controller

k = 20; k2 = 37; k3 = 70;

u_smc = @(t,x) smc_control(x, xd, k, k2, k3, m, M, I, g, Bx, l);
[t_smc,x_smc] = ode45(@(t,x) nonlinear_IP(x,u_smc(t,x),m,M,I,g,Bx,l), tspan, x0);

uSMC = zeros(length(t_smc),1);
for i=1:length(t_smc)
    uSMC(i) = u_smc(t_smc(i),x_smc(i,:)');
end

%%  ISMC Controller

lambda = 10; k_i = 3; zeta = 1; eps = 1;

[t_ismc,x_ismc] = ode45(@(t,y) pendulum_dynamics_ismc( ...
    t,y,M,m,l,I,g,Bx,xd,lambda,k_i,zeta,eps), tspan, x0);


%%  TSMC Controller
lambda = 20; k_t = 10; zeta = 20;
alpha = 30; beta = 40; p = 5; q = 3;

[t_tsmc,x_tsmc] = ode45(@(t,y) pendulum_dynamics_tsmc( ...
    t,y,M,m,l,I,g,Bx,xd,lambda,k_t,zeta,alpha,beta,p,q), tspan, x0);


%%                     PLOTS


% Figure 1: Cart Position
figure;
plot(t_lqr,x_lqr(:,1),'k','LineWidth',1.5); hold on;
plot(t_smc,x_smc(:,1),'r','LineWidth',1.5);
plot(t_ismc,x_ismc(:,1),'b','LineWidth',1.5);
plot(t_tsmc,x_tsmc(:,1),'m','LineWidth',1.5);
plot(tspan,r,'--','LineWidth',1.5);
xlabel('Time (s)'); ylabel('Cart Position (m)');
legend('LQR','SMC','ISMC','TSMC','Reference');
title('Cart Position Comparison'); grid on;

% Figure 2: Pendulum Angle
figure;
plot(t_lqr,x_lqr(:,3),'k','LineWidth',1.5); hold on;
plot(t_smc,x_smc(:,3),'r','LineWidth',1.5);
plot(t_ismc,x_ismc(:,3),'b','LineWidth',1.5);
plot(t_tsmc,x_tsmc(:,3),'m','LineWidth',1.5);
xlabel('Time (s)'); ylabel('\theta (rad)');
legend('LQR','SMC','ISMC','TSMC');
title('Pendulum Angle Comparison'); grid on;

% Figure 3: Control Input
figure;
plot(t_lqr,uLQR,'k','LineWidth',1.5); hold on;
plot(t_smc,uSMC,'r','LineWidth',1.5);
xlabel('Time (s)'); ylabel('Control Input (N)');
legend('LQR','SMC');
title('Control Effort'); grid on;

% Figure 4: Cart Position Error
figure;
plot(t_lqr,x_lqr(:,1)-xd,'k','LineWidth',1.5); hold on;
plot(t_smc,x_smc(:,1)-xd,'r','LineWidth',1.5);
plot(t_ismc,x_ismc(:,1)-xd,'b','LineWidth',1.5);
plot(t_tsmc,x_tsmc(:,1)-xd,'m','LineWidth',1.5);
xlabel('Time (s)'); ylabel('Error (m)');
legend('LQR','SMC','ISMC','TSMC');
title('Tracking Error Comparison'); grid on;


%%                  LOCAL FUNCTIONS

% (SMC)
function dx = nonlinear_IP(x,u,m,M,I,g,Bx,l)
    s = sin(x(3)); c = cos(x(3));
    D = I*(M+m) + M*m*l^2 + m^2*l^2*(1-c^2);

    dx = zeros(4,1);
    dx(1) = x(2);
    dx(2) = ((I + m*l^2)*(u - Bx*x(2) + m*l*s*x(4)^2) + m^2*l^2*g*s*c)/D;
    dx(3) = x(4);
    dx(4) = (-(m*l*c)*(u - Bx*x(2) + m*l*s*x(4)^2) + (M+m)*m*g*l*s)/D;
end

% SMC Control Law
function u = smc_control(x, xd, k, k2, k3, m, M, I, g, Bx, l)
    e  = x(1) - xd;
    ed = x(2);

    s = ed + k3*e;

    s_theta = sin(x(3)); c_theta = cos(x(3));
    D = I*(M+m) + M*m*l^2 + m^2*l^2*(1-c_theta^2);

    f = ((I + m*l^2)*(-Bx*x(2) + m*l*s_theta*x(4)^2) + m^2*l^2*g*s_theta*c_theta)/D;
    g_term = (I + m*l^2)/D;

    u_eq  = (-f - k3*ed)/g_term;
    u_dis = -k*sign(s);

    u = u_eq + u_dis;
end

% ISMC Dynamics
function dy = pendulum_dynamics_ismc(t,y,M,m,l,I,g,Bx,xd,lambda,k,zeta,eps)
    x = y(1); v = y(2); th = y(3); w = y(4);

    e = x - xd;
    s = w + lambda^3*e + 3*lambda^2*v + 3*lambda*th;

    if abs(s) > eps
        sat_s = sign(s);
    else
        sat_s = s/eps;
    end

    u = -k*sat_s - zeta*s - (lambda^3*v + 3*lambda^2*w);

    dy = zeros(4,1);
    dy(1) = v;
    dy(3) = w;

    A = [(M+m), m*l*cos(th); m*l*cos(th), (I+m*l^2)];
    B = [u + m*l*w^2*sin(th) - Bx*v; -m*g*l*sin(th)];
    acc = A\B;

    dy(2) = acc(1);
    dy(4) = acc(2);
end

%TSMC Dynamics
function dy = pendulum_dynamics_tsmc(t,y,M,m,l,I,g,Bx,xd,lambda,k,zeta,alpha,beta,p,q)
    x1=y(1); x2=y(2); x3=y(3); x4=y(4);

    den = M + m*sin(x3)^2;
    o2 = ((M+m)*g*sin(x3) - m*l*x4^2*cos(x3))/(l*den);
    i2 = -cos(x3)/(l*den);

    e = x1 - xd;
    s = x4 + lambda^3*e + 3*lambda^2*x2 + beta*sign(e)*abs(e)^(q/p);

    u_eq = -(o2 + lambda^3*x2 + beta*(q/p)*abs(e)^(q/p-1)*x2)/i2;
    u_dis = -k*sign(s) - zeta*s;
    F = u_eq + u_dis;

    dy = zeros(4,1);
    dy(1)=x2; dy(3)=x4;

    A=[(M+m), m*l*cos(x3); m*l*cos(x3), (I+m*l^2)];
    B=[F + m*l*x4^2*sin(x3) - Bx*x2; -m*g*l*sin(x3)];
    acc=A\B;

    dy(2)=acc(1);
    dy(4)=acc(2);
end