clear; clc; close all;

tic


K = 1.0e+04 *[0.0393   -1.2409   -0.0013];
K = 1.0e+04 *[-9.1155   -4.3490   -0.0973];
K = 1.0e+04 *[6.7696   -9.0340   -0.1148];

Kp = K(2);
Kd = K(3);
alpha = K(1);
F_max = 300;inf;4.5;

q1d=pi;

% % % Parametros iniciales del problema
m1 = 1;
m2 = 1;
lc1 = 0.5;
lc2 = 0.5;
l1 = 1;
l2 = 1;
I1 = 1;
I2 = 1;
g = 9.8;
dt_controller = 0.002;
dt = dt_controller/4;

t_total = 40;

% % % Inicio de parametros % % %
q1_ini = rand*100;0.01;pi+rand*2.8;
q2_ini = rand*100;0;
dq1_ini = rand*100;0;
dq2_ini = rand*100;0;
ddq1_ini = 0;
ddq2_ini = 0;
F_ini = 0;

F = F_ini; q1 = q1_ini; q2 = q2_ini; dq1 = dq1_ini; dq2 = dq2_ini;
q1_t = [];
q1d_t = [];
q2_t = [];
F_t = [];

for t = 0:dt_controller:t_total
    
    q1_ant = q1; q2_ant = q2; dq1_ant = dq1; dq2_ant = dq2;
    
    %% controlador

    
    d11 = m1*lc1^2 + m2*(l1^2+lc2^2+2*l1*lc2*cos(q2)) + I1 + I2;
    d22 = m2*lc2^2 + I2;
    d12 = m2*(lc2^2 + l1*lc2*cos(q2)) + I2;
    h1 = -m2*l1*lc2*dq2^2*sin(q2) - 2*m2*l1*lc2*dq2*dq1*sin(q2);
    h2 = m2*l1*lc2*dq1^2*sin(q2);
    phi2 = m2*lc2*g*cos(q1+q2-pi/2);
    phi1 = (m1*lc1 + m2*l1)*g*cos(q1-pi/2) + phi2;
    
    
    d22_ = d22 - d12^2/d11;
    h2_ = h2 - d12/d11*h1;
    phi2_ = phi2 - d12/d11*phi1;
        
    error = mod(q2-(q1-q1d)-pi,2*pi)-pi;
    Act = Kp*(error) - Kd*(dq2 - dq1);

    F = d22_ * Act + h2_ +phi2_;
    
    F = max(-F_max, min(F_max, F));
    
    
    
    %% bucle de robot
    for times = 1:dt_controller/dt
        d11 = m1*lc1^2 + m2*(l1^2+lc2^2+2*l1*lc2*cos(q2)) + I1 + I2;
        d22 = m2*lc2^2 + I2;
        d12 = m2*(lc2^2 + l1*lc2*cos(q2)) + I2;
        h1 = -m2*l1*lc2*dq2^2*sin(q2) - 2*m2*l1*lc2*dq2*dq1*sin(q2);
        h2 = m2*l1*lc2*dq1^2*sin(q2);
        phi2 = m2*lc2*g*cos(q1+q2-pi/2);
        phi1 = (m1*lc1 + m2*l1)*g*cos(q1-pi/2) + phi2;
        
        
        ddq2 = (F + d12/d11 * (h1 + phi1) - h2 - phi2) / (d22 - d12^2/d11);
        ddq1 = -(d12* ddq2 + h1 + phi1) / d11;
        
        dq1 = max(-4*pi, min(4*pi, dq1 + ddq1*dt));
        dq2 = max(-9*pi, min(9*pi, dq2 + ddq2*dt));
        q1 = mod(q1 + dq1*dt, 2*pi);
        q2 = mod(q2 + dq2*dt, 2*pi);
        
%         q1d = mod(q1d+0.003, 2*pi);
        
        q1_t = [q1_t q1];
        q1d_t = [q1d_t q1d];
        q2_t = [q2_t q2];
        F_t = [F_t F];
    end
    
    
    
    
    
    
    
end



for i=1:length(q1_t)
    x1(i) = l1 * cos( q1_t(i) - pi/2 );
    y1(i) = l1 * sin( q1_t(i) - pi/2 );
    x2(i) = x1(i) + l2 * cos( q1_t(i) + q2_t(i) - pi/2 );
    y2(i) = y1(i) + l2 * sin( q1_t(i) + q2_t(i) - pi/2 );
    
end
% 
figure(1);
 plot(1:length(q1_t),q1_t,1:length(q1d_t),q1d_t);
% figure(2);
% plot(1:length(y1),y1,1:length(y2),y2);
% figure(4);
% plot(x1,y1,x2,y2);
% figure(5);
% plot(1:length(F_t),F_t);

for i=1:dt_controller*10/dt:length(q1_t)
    figure(3);
    plot([0;x1(i)],[0;y1(i)],[x1(i);x2(i)],[y1(i);y2(i)]);
    xlim([-2.5 2.5]);
    ylim([-2.5 2.5]);
    pause(0.03);
end




toc