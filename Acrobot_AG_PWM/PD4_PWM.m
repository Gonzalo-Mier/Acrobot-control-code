clear; clc; close all;

tic


K = 1.0e+05 *[0.3913    0.1278    1.0200    0.1647   -0.3239   -0.1644   -0.7094    0.4668]


K1 = K(1);
K2 = K(2);
K3 = K(3);
K4 = K(4);
Ki1 = K(5);
Ki2 = K(6);
Ki3 = K(7);
Ki4 = K(8);
F_max = 300;inf;4.5;

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
dt_controller = 0.05;
dt = dt_controller/4;
filtro = 2/dt;
t_total = 40;

% % % Inicio de parametros % % %
q1_ini = 0.01;
q2_ini = 0;
dq1_ini = 0;
dq2_ini = 0;
ddq1_ini = 0;
ddq2_ini = 0;
F_ini = 0;

F = F_ini; q1 = q1_ini; q2 = q2_ini; dq1 = dq1_ini; dq2 = dq2_ini;
q1_t = [];
q1d_t = [];
q2_t = [];
q2d_t = [];
F_t = [];

Ei1 = 0;
Ei2 = 0;
Ei3 = 0;
Ei4 = 0;

q1d = pi/2;
q2d = 0;
dq1d = 0;
dq2d = 0;


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
        
%       E1 = mod((q1d-q1)-pi,2*pi)-pi;
%     E2 = mod((q2d-q2)-pi,2*pi)-pi;
    E1 = q1d-q1;
    E2 = q2d-q2;    
    E3 = dq1d-dq1;
    E4 = dq2d-dq2;
    Ei1 = Ei1 + E1;
    Ei2 = Ei2 + E2;
    Ei3 = Ei3 + E3;
    Ei4 = Ei4 + E4; 
    
    Act = K1*E1 + K2*E2 + K3*E3 + K4*E4 + Ki1*Ei1 + Ki2*Ei2 + Ki3*Ei3 + Ki4*Ei4 ;

    F = Act; d22_ * Act + h2_ +phi2_;
    
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
        
%         q1d = mod(q1d+0.5*dt, 2*pi);
%         q1d = q1d+0.2*dt;
        
% if t<t_total/4
%     q1d = 0;
% elseif t<t_total*3/4
%     q1d = pi * (t - t_total/4)/(t_total/2)
% else
%     q1d = pi;
% end
%     







        q1_t = [q1_t q1];
        q1d_t = [q1d_t q1d];        
        q2_t = [q2_t q2];
        q2d_t = [q2d_t q2d];             
        F_t = [F_t F];
    end
    
    
    
    
    
    
    
end



for i=1:length(q1_t)
    x1(i) = l1 * cos( q1_t(i) - pi/2 );
    y1(i) = l1 * sin( q1_t(i) - pi/2 );
    x2(i) = x1(i) + l2 * cos( q1_t(i) + q2_t(i) - pi/2 );
    y2(i) = y1(i) + l2 * sin( q1_t(i) + q2_t(i) - pi/2 );
    
end

for i = 1:length(F_t)-filtro
    F_m(i) = sum(F_t(i:i+filtro))/(filtro+1); 
end

% 
figure(1);
plot(1:length(q1_t),q1_t,1:length(q1d_t),q1d_t);
figure(7);
plot(1:length(F_t),F_t,'b',1:length(F_m),F_m,'r');

figure(2);
plot(1:length(y1),y1,1:length(y2),y2);
figure(4);
plot(x1,y1,x2,y2);

for i=1:dt_controller/dt:length(q1_t)
    figure(3);
    plot([0;x1(i)],[0;y1(i)],[x1(i);x2(i)],[y1(i);y2(i)]);
    xlim([-2.5 2.5]);
    ylim([-2.5 2.5]);
    pause(0.03);
end




toc