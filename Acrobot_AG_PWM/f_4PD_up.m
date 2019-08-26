function [ fitness ] = f_4PD_up( K, p_ini, params )

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

q1 = p_ini(1);
q2 = p_ini(2);
dq1 = p_ini(3);
dq2 = p_ini(4);
ddq1 = p_ini(5);
ddq2 = p_ini(6);
F = p_ini(7);

K1 = K(1);
K2 = K(2);
K3 = K(3);
K4 = K(4);
Ki1 = K(5);
Ki2 = K(6);
Ki3 = K(7);
Ki4 = K(8);


F_t = [];
q1_t = [];
q2_t = [];

Ei1 = 0;
Ei2 = 0;
Ei3 = 0;
Ei4 = 0;

q1d = pi;
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
    
    
%     E1 = mod((q1d-q1)-pi,2*pi)-pi;
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
        q1 = q1 + dq1*dt;
        q2 = q2 + dq2*dt;
        
%         q1 = mod(q1 + dq1*dt - pi, 2*pi)-pi;
%         q2 = mod(q2 + dq2*dt - pi, 2*pi)-pi;
        
        
        F_t = [F_t F];
        q1_t = [q1_t q1];
        q2_t = [q2_t q2];
    end
    y1 = l1 * sin( q1 - pi/2 );
    y2 = y1 + l2 * sin( q1 + q2 - pi/2 );
    if y2 <= 0.0
        break;
    end
    
    
end

if y2 <= 0.0
    fitness = 10000000-length(q1_t);
else
    
    
    
    for i=1:length(q1_t)
        x1(i) = l1 * cos( q1_t(i) - pi/2 );
        y1(i) = l1 * sin( q1_t(i) - pi/2 );
        x2(i) = x1(i) + l2 * cos( q1_t(i) + q2_t(i) - pi/2 );
        y2(i) = y1(i) + l2 * sin( q1_t(i) + q2_t(i) - pi/2 );
    end
    
    for i=1:length(q1_t)-1
        c_x2(i) = (x2(i+1) - x2(i)).^2;
        c_y2(i) = (y2(i+1) - y2(i)).^2;
        c_x1(i) = (x1(i+1) - x1(i)).^2;
        c_y1(i) = (y1(i+1) - y1(i)).^2;
        
        d_x2(i) = (x2(i) - 0).^2;
        d_y2(i) = (y2(i) - 2).^2;
        d_x1(i) = (x1(i) - 0).^2;
        d_y1(i) = (y1(i) - 1).^2;
    end
    w = [1 1 0.1 0.1 1 10 0.1 0.1 1/10000]/length(q1_t);
    w = [0 0 0 0 1 10 0.1 0.1 1/10000]/length(q1_t);
    
    
    fitness = w*[sum(c_x2) ; sum(c_y2) ; sum(c_x1) ; sum(c_y1) ; sum(d_x2) ; sum(d_y2) ;sum(d_x1) ; sum(d_y1) ; sum(abs(F_t))];
    
end


end