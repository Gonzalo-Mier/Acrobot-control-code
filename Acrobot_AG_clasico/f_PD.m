function [ fitness ] = f_PD( K, p_ini, params, params_est )

m1 = params(1);
m2 = params(2);
lc1 = params(3);
lc2 = params(4);
l1 = params(5);
l2 = params(6);
I1 = params(7);
I2 = params(8);
g = params(9);
dt_controller = params(10);
dt = params(11);
t_total = params(12);
F_max = params(13);

m1_est = params_est(1);
m2_est = params_est(2);
lc1_est = params_est(3);
lc2_est = params_est(4);
l1_est = params_est(5);
l2_est = params_est(6);
I1_est = params_est(7);
I2_est = params_est(8);
error_m = params_est(9);



q1 = p_ini(1);
q2 = p_ini(2);
dq1 = p_ini(3);
dq2 = p_ini(4);
ddq1 = p_ini(5);
ddq2 = p_ini(6);
F = p_ini(7);

alpha = K(1);
Kp = K(2);
Kd = K(3);


F_t = [];
q1_t = [];
q2_t = [];
q2d = 0;

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
%         q1 = mod(q1 + dq1*dt - pi, 2*pi)-pi;
%         q2 = mod(q2 + dq2*dt - pi, 2*pi)-pi;
        q1 = q1 + dq1*dt;
        q2 = q2 + dq2*dt;       
        
        F_t = [F_t F];
        q1_t = [q1_t q1];
        q2_t = [q2_t q2];
    end
    
    y1 = l1 * sin( q1 - pi/2 );
    y2 = y1 + l2 * sin( q1 + q2 - pi/2 );
    
    
    if y2 >= 1.0
        break;
    end
    
    
end

fitness = length(q1_t);


end

