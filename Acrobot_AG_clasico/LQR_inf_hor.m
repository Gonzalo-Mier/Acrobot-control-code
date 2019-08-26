clear; clc; close all;

tic


% K=  -1.0e+03 * [-1.5170   -0.4435   -0.9808   -0.3690];

F_max = 10;inf;4.5;



% % % Parametros iniciales del problema
m1 = 1;
m2 = 1;
lc1 = 0.5;
lc2 = 0.5;
l1 = 1;
l2 = 1;
I1 = 1;
I2 = 1;
g = 9.81;
dt_controller = 0.01;
dt = dt_controller/4;
t_total = 10;


% % % Inicio de parametros % % %

q1_ini = pi+.04;
q2_ini = -.05;
dq1_ini = -.2;
dq2_ini = +0.04;
ddq1_ini = 0;
ddq2_ini = 0;
F_ini = 0;

F = F_ini; q1 = q1_ini; q2 = q2_ini; dq1 = dq1_ini; dq2 = dq2_ini;
F_t = [];
q1_t = [];
q2_t = [];

q1d = pi;
q2d = 0;
dq1d = 0;
dq2d= 0;


p_ini = [q1_ini, q2_ini, dq1_ini, dq2_ini, ddq1_ini, ddq2_ini, F_ini];
params = [m1, m2, lc1, lc2, l1, l2, I1, I2, g, dt_controller, dt, t_total, F_max];


                 
 
m11 = m1*lc1^2+I1+m2*lc2^2+I2+m2*l1^2+2*m2*l1*lc2;
m22 = m2*lc2^2+I2;
m12 = m2*lc2^2+I2+m2*l1*lc2;
m21 = m12;
g_1 = -(m1*lc1+m2*l1+ m2*lc2)*g;
g_2 = -m2*lc2*g;
Z = m11*m22 - m12*m21;                 


A = [                       0,               0, 1, 0; ...
                            0,               0, 0, 1; ...
         (-m22*g_1+m12*g_2)/Z, (m12-m22)*g_2/Z, 0, 0; ...
         (+m21*g_1-m11*g_2)/Z, (m21-m11)*g_2/Z, 0, 0];
                 
B = [0; 0; -m12/Z; m11/Z];

Q = [1000, -500, 0 0; -500, 1000, 0, 0; 0, 0, 1000, -500; 0 0 -500 1000];

R = [10000];


K = lqr(A,B,Q,R);









for t = 0:dt_controller:t_total
    
    q1_ant = q1; q2_ant = q2; dq1_ant = dq1; dq2_ant = dq2;
    
    %% controlador
    x= [mod(q1,2*pi)-pi; mod(q2-pi,2*pi)-pi; dq1; dq2];
    F = -K*x;
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
        q1 = mod(q1 + dq1*dt - pi, 2*pi)-pi;
        q2 = mod(q2 + dq2*dt - pi, 2*pi)-pi;
        
        
        F_t = [F_t F];
        q1_t = [q1_t q1];
        q2_t = [q2_t q2];
    end
    
    

    
    
    
end

for i=1:length(q1_t)
    x1(i) = l1 * cos( q1_t(i) - pi/2 );
    y1(i) = l1 * sin( q1_t(i) - pi/2 );
    x2(i) = x1(i) + l2 * cos( q1_t(i) + q2_t(i) - pi/2 );
    y2(i) = y1(i) + l2 * sin( q1_t(i) + q2_t(i) - pi/2 );
    
end
% 
% 
% 
% figure(4);
% plot(1:length(q1_t),mod(q1_t,2*pi)-pi,1:length(q2_t),q2_t);
% legend('Primera art.','Segunda art.')
% title('Error del controlador')
% xlabel('Tiempo (s)')
% ylabel('Error (rad)')
% ylim([-pi,pi]);
% 
% 
% fitness = sqrt((median(x2)-0)^2 + (median(y2)-2)^2)
% 
% figure(5);
% plot(x1,y1,x2,y2,0,0,'o');
% legend('Primera art.','Segunda art.','Base del robot')
% title('Recorrido del brazo')
% xlabel('Desplazamiento en x (m)')
% ylabel('Desplazamiento en y (m)')
% ylim([-.1,2.1]);
% xlim([-2.1,2.1]);
% 
% figure(6);
% plot(F_t);


% figure(7);
% xlim([-2.5 2.5]);
% ylim([-2.5 2.5]);
% % for i=1:dt_controller/dt:length(q1_t)
% for i=1:5:length(q1_t)
%     figure(7);
%     plot([0;x1(i)],[0;y1(i)],[x1(i);x2(i)],[y1(i);y2(i)]);
%     xlim([-2.5 2.5]);
%     ylim([-2.5 2.5]);
%     pause(0.03);
% end


Carpeta = 'LQR_hor_inf';
grafica_res



toc