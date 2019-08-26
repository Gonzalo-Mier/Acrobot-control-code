% clear; clc; close all;
% 
% tic
% 
% K = [1 44 2.1]; %Hecho a mano 
% K = [-3.1547   57.1824    7.9877]; %%% 656
% K = [1.3658   53.4288    0.2733];
% K = [-1.9800   47.8122    2.9995];

alpha = K(1);
Kp = K(2);
Kd = K(3);
% F_max = 1;4.5;

% % % Parametros iniciales del problema
% m1 = 1;
% m2 = 1;
% lc1 = 0.5;
% lc2 = 0.5;
% l1 = 1;
% l2 = 1;
% I1 = 1;
% I2 = 1;
% g = 9.8;
% dt_controller = 0.01;
% dt = dt_controller/4;

t_total = 40;


% % % Inicio de parametros % % %
q1_ini = 0+rand*.1;
q2_ini = 0;
dq1_ini = 0;
dq2_ini = 0;
ddq1_ini = 0;
ddq2_ini = 0;
F_ini = 0;

F = F_ini; q1 = q1_ini; q2 = q2_ini; dq1 = dq1_ini; dq2 = dq2_ini;
q1_t = [];
q2_t = [];
F_t = [];

x1 = [];
x2 = [];
y1 = [];
y2 = [];


q1d_t = [];
q2d_t = [];

for t = 0:dt_controller:t_total
    
    q1_ant = q1; q2_ant = q2; dq1_ant = dq1; dq2_ant = dq2;
    
    %% controlador
 
%     q1_m = mod(q1 +(2*rand-1)*2*pi*error_m - pi,2*pi)-pi;
%     q2_m = mod(q2 +(2*rand-1)*2*pi*error_m - pi,2*pi)-pi; 
 
    q1_m = q1 +(2*rand-1)*2*pi*error_m;
    q2_m = q2 +(2*rand-1)*2*pi*error_m;    
    
    dq1_m =max(-4*pi, min(4*pi, dq1 +(2*rand-1)*8*pi*error_m));
    dq2_m = max(-9*pi, min(9*pi,dq2 +(2*rand-1)*18*pi*error_m));
    
    
    d11_est = m1_est*lc1_est^2 + m2_est*(l1_est^2+lc2_est^2+2*l1_est*lc2_est*cos(q2_m)) + I1_est + I2_est;
    d22_est = m2_est*lc2_est^2 + I2_est;
    d12_est = m2_est*(lc2_est^2 + l1_est*lc2_est*cos(q2_m)) + I2_est;
    h1_est = -m2_est*l1_est*lc2_est*dq2_m^2*sin(q2_m) - 2*m2_est*l1_est*lc2_est*dq2_m*dq1_m*sin(q2_m);
    h2_est = m2_est*l1_est*lc2_est*dq1_m^2*sin(q2_m);
    phi2_est = m2_est*lc2_est*g*cos(q1_m+q2_m-pi/2);
    phi1_est = (m1_est*lc1_est + m2_est*l1_est)*g*cos(q1_m-pi/2) + phi2_est;
    
    
    d22__est = d22_est - d12_est^2/d11_est;
    h2__est = h2_est - d12_est/d11_est*h1_est;
    phi2__est = phi2_est - d12_est/d11_est*phi1_est;
    
    q2d = alpha*atan(dq1_m);
    q2d_t = [q2d_t q2d]; 
    
    
    %     F = d22__est * (Kp*(alpha*atan(dq1_m)- q2) - Kd*dq2_m) + h2__est +phi2__est;
    
    q2d1 = q2d;
    q2d = alpha*atan(dq1_m);
    dq2d = (q2d - q2d1)/dt_controller;
    F = d22__est * (Kp*(q2d - q2) + Kd*(dq2d - dq2_m)) + h2__est +phi2__est;

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
%         q1 = mod(q1 + dq1*dt, 2*pi);
%         q2 = mod(q2 + dq2*dt, 2*pi);
        q1 = q1 + dq1*dt;
        q2 = q2 + dq2*dt;       
        
        
        q1_t = [q1_t q1];
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
% figure(1);
% plot(1:length(q1_t),q1_t,1:length(q1_t),q2_t);
% figure(2);
% plot(1:length(y1),y1,1:length(y2),y2);
% figure(4);
% plot(x1,y1,x2,y2);
%     pause(0.03);
% figure(3);
% xlim([-2.5 2.5]);
% ylim([-2.5 2.5]);
% for i=1:dt_controller/dt:length(q1_t)
%     figure(3);
%     plot([0;x1(i)],[0;y1(i)],[x1(i);x2(i)],[y1(i);y2(i)]);
%     xlim([-2.5 2.5]);
%     ylim([-2.5 2.5]);
%     pause(0.03);
% 
%     
% end




toc