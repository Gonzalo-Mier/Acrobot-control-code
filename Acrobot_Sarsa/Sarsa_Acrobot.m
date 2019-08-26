% % Modelado del problema del Acrobot
% clear all,
% close all, 
% clc


% % % % Parametros iniciales del problema
% m1 = 1;
% m2 = 1;
% lc1 = 0.5;
% lc2 = 0.5;
% l1 = 1;
% l2 = 1;
% I1 = 1;
% I2 = 1;
% g = 9.8;
% dt = 0.05;
% t_total = 1000;
% 
% % % % Inicio de parametros % % %
% q1_ini = 0;
% q2_ini = 0;
% dq1_ini = 0;
% dq2_ini = 0;
% ddq1_ini = 0;
% ddq2_ini = 0;
% F_ini = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%                              %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%             SARSA            %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%                              %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Parametros Controlador
% 
% alpha = 0.2/48;
% lambda = 0.9;
% gamma = 1;
% eps = 0;
% MaxIt=600;


%%%% Variables para el controlador SARSA %%%%%%%%



actionList = [-1.0; 0; 1.0]*1;



for tablas = 1:1
    for tabla_4v =1:12
        Q_tabla{tabla_4v,1} = [0, 2*pi/5*rand:2*pi/5:2*pi];
        Q_tabla{tabla_4v,2} = [0, 2*pi/5*rand:2*pi/5:2*pi];
        Q_tabla{tabla_4v,3} = [0,8*pi/6*rand+0:8*pi/6:8*pi] - 4*pi;
        Q_tabla{tabla_4v,4} = [0,18*pi/6*rand+0:18*pi/6:18*pi] - 9*pi;
    end
    
    for tabla_3v =1:3
        Q_tabla{(tabla_3v-1)*4 + 13,1} = [0, 2*pi/5*rand:2*pi/5:2*pi];
        Q_tabla{(tabla_3v-1)*4 + 13,2} = [0, 2*pi/5*rand:2*pi/5:2*pi];
        Q_tabla{(tabla_3v-1)*4 + 13,3} = [0,8*pi/6*rand+0:8*pi/6:8*pi] - 4*pi;
        Q_tabla{(tabla_3v-1)*4 + 13,4} = -inf;
        
        Q_tabla{(tabla_3v-1)*4 + 14,1} = [0, 2*pi/5*rand:2*pi/5:2*pi];
        Q_tabla{(tabla_3v-1)*4 + 14,2} = [0, 2*pi/5*rand:2*pi/5:2*pi];
        Q_tabla{(tabla_3v-1)*4 + 14,3} = -inf;
        Q_tabla{(tabla_3v-1)*4 + 14,4} = [0,18*pi/6*rand+0:18*pi/6:18*pi] - 9*pi;
        
        Q_tabla{(tabla_3v-1)*4 + 15,1} = [0, 2*pi/5*rand:2*pi/5:2*pi];
        Q_tabla{(tabla_3v-1)*4 + 15,2} = -inf;
        Q_tabla{(tabla_3v-1)*4 + 15,3} = [0,8*pi/6*rand+0:8*pi/6:8*pi] - 4*pi;
        Q_tabla{(tabla_3v-1)*4 + 15,4} = [0,18*pi/6*rand+0:18*pi/6:18*pi] - 9*pi;
        
        Q_tabla{(tabla_3v-1)*4 + 16,1} = -inf;
        Q_tabla{(tabla_3v-1)*4 + 16,2} = [0, 2*pi/5*rand:2*pi/5:2*pi];
        Q_tabla{(tabla_3v-1)*4 + 16,3} = [0,8*pi/6*rand+0:8*pi/6:8*pi] - 4*pi;
        Q_tabla{(tabla_3v-1)*4 + 16,4} = [0,18*pi/6*rand+0:18*pi/6:18*pi] - 9*pi;
    end
    
    for tabla_2v =1:2
        Q_tabla{(tabla_2v-1)*6 + 25,1} = [0, 2*pi/5*rand:2*pi/5:2*pi];
        Q_tabla{(tabla_2v-1)*6 + 25,2} = [0, 2*pi/5*rand:2*pi/5:2*pi];
        Q_tabla{(tabla_2v-1)*6 + 25,3} = -inf;
        Q_tabla{(tabla_2v-1)*6 + 25,4} = -inf;
        
        Q_tabla{(tabla_2v-1)*6 + 26,1} = [0, 2*pi/5*rand:2*pi/5:2*pi];
        Q_tabla{(tabla_2v-1)*6 + 26,2} = -inf;
        Q_tabla{(tabla_2v-1)*6 + 26,3} = [0,8*pi/6*rand+0:8*pi/6:8*pi] - 4*pi;
        Q_tabla{(tabla_2v-1)*6 + 26,4} = -inf;
        
        Q_tabla{(tabla_2v-1)*6 + 27,1} = [0, 2*pi/5*rand:2*pi/5:2*pi];
        Q_tabla{(tabla_2v-1)*6 + 27,2} = -inf;
        Q_tabla{(tabla_2v-1)*6 + 27,3} = -inf;
        Q_tabla{(tabla_2v-1)*6 + 27,4} = [0,18*pi/6*rand+0:18*pi/6:18*pi] - 9*pi;
        
        Q_tabla{(tabla_2v-1)*6 + 28,1} = -inf;
        Q_tabla{(tabla_2v-1)*6 + 28,2} = [0, 2*pi/5*rand:2*pi/5:2*pi];
        Q_tabla{(tabla_2v-1)*6 + 28,3} = [0,8*pi/6*rand+0:8*pi/6:8*pi] - 4*pi;
        Q_tabla{(tabla_2v-1)*6 + 28,4} = -inf;
        
        Q_tabla{(tabla_2v-1)*6 + 29,1} = -inf;
        Q_tabla{(tabla_2v-1)*6 + 29,2} = [0, 2*pi/5*rand:2*pi/5:2*pi];
        Q_tabla{(tabla_2v-1)*6 + 29,3} = -inf;
        Q_tabla{(tabla_2v-1)*6 + 29,4} = [0,18*pi/6*rand+0:18*pi/6:18*pi] - 9*pi;
        
        Q_tabla{(tabla_2v-1)*6 + 30,1} = -inf;
        Q_tabla{(tabla_2v-1)*6 + 30,2} = -inf;
        Q_tabla{(tabla_2v-1)*6 + 30,3} = [0,8*pi/6*rand+0:8*pi/6:8*pi] - 4*pi;
        Q_tabla{(tabla_2v-1)*6 + 30,4} = [0,18*pi/6*rand+0:18*pi/6:18*pi] - 9*pi;
    end
    
    for tabla_1v =1:3
        Q_tabla{(tabla_1v-1)*4 + 37,1} = [0, 2*pi/5*rand:2*pi/5:2*pi];
        Q_tabla{(tabla_1v-1)*4 + 37,2} = -inf;
        Q_tabla{(tabla_1v-1)*4 + 37,3} = -inf;
        Q_tabla{(tabla_1v-1)*4 + 37,4} = -inf;
        
        Q_tabla{(tabla_1v-1)*4 + 38,1} = -inf;
        Q_tabla{(tabla_1v-1)*4 + 38,2} = [0, 2*pi/5*rand:2*pi/5:2*pi];
        Q_tabla{(tabla_1v-1)*4 + 38,3} = -inf;
        Q_tabla{(tabla_1v-1)*4 + 38,4} = -inf;
        
        Q_tabla{(tabla_1v-1)*4 + 39,1} = -inf;
        Q_tabla{(tabla_1v-1)*4 + 39,2} = -inf;
        Q_tabla{(tabla_1v-1)*4 + 39,3} = [0,8*pi/6*rand+0:8*pi/6:8*pi] - 4*pi;
        Q_tabla{(tabla_1v-1)*4 + 39,4} = -inf;
        
        Q_tabla{(tabla_1v-1)*4 + 40,1} = -inf;
        Q_tabla{(tabla_1v-1)*4 + 40,2} = -inf;
        Q_tabla{(tabla_1v-1)*4 + 40,3} = -inf;
        Q_tabla{(tabla_1v-1)*4 + 40,4} = [0,18*pi/6*rand+0:18*pi/6:18*pi] - 9*pi;
        
        
    end
    
end



for tabla = 1:48
    Q_matrix{tabla} = zeros(length(Q_tabla{tabla,1}),length(Q_tabla{tabla,2}),length(Q_tabla{tabla,3}),length(Q_tabla{tabla,4}),length(actionList));
end


% % % Ejecutar modelo
for it=1:MaxIt

    total_reward=0;
    flag = 0; first_move = 1;
    q1 = q1_ini; q2 = q2_ini; dq1 = dq1_ini; dq2 = dq2_ini;
    
    w = zeros(1,length(actionList));
    for tabla = 1:48
        I(tabla) = find( q1 >= Q_tabla{tabla,1}, 1, 'last' );
        J(tabla) = find( q2 >= Q_tabla{tabla,2}, 1, 'last' );
        K(tabla) = find( dq1 >= Q_tabla{tabla,3}, 1, 'last' );
        L(tabla) = find( dq2 >= Q_tabla{tabla,4}, 1, 'last' );
        w = w + reshape(Q_matrix{tabla}(I(tabla),J(tabla),K(tabla),L(tabla),:),1,length(actionList));
    end
    
    if (rand >= eps)
        [~, pos] = max(w);
        F = actionList(pos);
    else
        pos = randi(length(actionList));
        F = actionList(pos);
    end
    
    
    
    
    t=[]; q1_t = []; q2_t = []; F_t=[]; n_pasos = 0;
    moves = 1;
    
    for i = 0:dt*4:t_total
        
        I_a = I; J_a = J; K_a = K; L_a = L;  pos_a = pos; 
        
        for j=1:4
            q1_ant = q1; q2_ant = q2; dq1_ant = dq1; dq2_ant = dq2;
            
            
            phi2 = m2*lc2*g*cos(q1+q2-pi/2);
            phi1 = -m2*l1*lc2*dq2^2*sin(q2) - 2*m2*l1*lc2*dq2*dq1*sin(q2) + (m1*lc1 + m2*l1)*g*cos(q1-pi/2) + phi2;
            d2 = m2*(lc2^2 + l1*lc2*cos(q2)) + I2;
            d1 = m1*lc1^2 + m2*(l1^2+lc2^2+2*l1*lc2*cos(q2)) + I1 + I2;
            ddq2_n = (F + d2/d1 * phi1 - m2*l1*lc2*dq1^2*sin(q2) - phi2) / (m2*lc2^2 + I2 - d2^2/d1);
            ddq1_n = -(d2* ddq2_n + phi1) / d1;
            
            dq1 = max(-4*pi, min(4*pi, dq1 + ddq1_n*dt));
            dq2 = max(-9*pi, min(9*pi, dq2 + ddq2_n*dt));
            q1 = mod(q1 + dq1*dt, 2*pi);

            q2 = mod(max(min(mod(q2-pi,2*pi) + dq2*dt-pi,pi-1^-6),-pi+1^-6),2*pi);
%             q2 = mod(q2-pi,2*pi) + dq2*dt -pi;
%             if  q2 > pi-1^-6
%                 q2 = pi-1^-6;
%                 dq2 = 0;
%             elseif  q2 < -pi+1^-6
%                 q2 = -pi+1^-6;
%                 dq2 = 0;
%             end
%             q2 = mod(q2,2*pi);
        end
        
        
        F_ant = F;
        w_a = w;
        w = zeros(1,length(actionList));
        for tabla = 1:48
            I(tabla) = find( q1 >= Q_tabla{tabla,1}, 1, 'last' );
            J(tabla) = find( q2 >= Q_tabla{tabla,2}, 1, 'last' );
            K(tabla) = find( dq1 >= Q_tabla{tabla,3}, 1, 'last' );
            L(tabla) = find( dq2 >= Q_tabla{tabla,4}, 1, 'last' );
            w = w + reshape(Q_matrix{tabla}(I(tabla),J(tabla),K(tabla),L(tabla),:),1,length(actionList));
        end
        
        
        if (rand >= eps)
            [~, pos] = max(w);
            F = actionList(pos);
        else
            pos = randi(length(actionList));
            F = actionList(pos);
        end
        
        if first_move == 0
            if flag == 0
                y1 = sin( q1 - pi/2 );
                y2 = y1 + sin( q1 + q2 - pi/2 );
                if y2 >= 1.0
                    flag = 1;
                    reward = 0; %t_total*dt/4;
                    
                else
                    reward = -1;%/t_total*dt/4;
                end
                
                total_reward = total_reward + reward;
                
                for tabla = 1:48
                    Q_matrix{tabla}(I_a(tabla),J_a(tabla),K_a(tabla),L_a(tabla),pos_a) =  ...
                        Q_matrix{tabla}(I_a(tabla),J_a(tabla),K_a(tabla),L_a(tabla),pos_a) + ...
                        alpha * ( reward + gamma.*w(pos) - w_a(pos_a));
                end
                
            else
                break;
            end
            
        else
            first_move = 0;
        end
        
        
        
        if F_ant ~= F
            moves = moves+1;
        end
        
        t = [t i];
        q1_t = [q1_t q1];
        q2_t = [q2_t q2];
        F_t = [F_t F];
        n_pasos = n_pasos + 1;
        
    end
    reward_t(it) = total_reward;
    it_t(it) = n_pasos;
    eps = eps*0.99;
    
end




for i=1:length(t)
    x1(i) = l1 * cos( q1_t(i) - pi/2 );
    y1(i) = l1 * sin( q1_t(i) - pi/2 );
    x2(i) = x1(i) + l2 * cos( q1_t(i) + q2_t(i) - pi/2 );
    y2(i) = y1(i) + l2 * sin( q1_t(i) + q2_t(i) - pi/2 );
end
q1_t = mod(q1_t-pi,2*pi)-pi;
q2_t = mod(q2_t-pi,2*pi)-pi;

for i=1:MaxIt-9
    it_median(i) = sum(it_t(i:i+9))/10;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(3);
% xlim([-2.5 2.5]);
% ylim([-2.5 2.5]);
% for i=1:length(t)
%     plot([0;x1(i)],[0;y1(i)],[x1(i);x2(i)],[y1(i);y2(i)]);
%     xlim([-2.5 2.5]);
%     ylim([-2.5 2.5]);
%     pause(0.03);
% end
% % % % % % %
% 
% figure(1);
% plot(t,q1_t,t,q2_t);
% figure(6);
% plot(1:it_t(end),F_t);
% figure(4);
% plot(x1,y1,x2,y2);

% figure(2);
% plot(1:MaxIt, reward_t);
% figure(5);
% plot(1:MaxIt, it_t);
% 
% figure(7);
% plot(1:MaxIt-9, it_median,1:MaxIt, it_t);


