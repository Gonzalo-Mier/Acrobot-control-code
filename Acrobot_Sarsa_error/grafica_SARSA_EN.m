close all

figure(6);

semilogy((1:MaxIt),it_S_exp,1:length(it_S_exp_median),it_S_exp_median);
legend('Steps', 'Filtered steps')
title('Steps per Episode')
xlabel('Episode')
ylabel('Steps')
 ylim([10, 100000])


% 
% 
% q1d1 =  mod(pi-0.8-pi,2*pi)-pi;
% q1d2 =  mod(pi+0.8-pi,2*pi)-pi;
% 
% figure('units','normalized','outerposition',[0 0 1 1]);
% 
% subplot(3,4,[1,2]);
% plot(t, mod(q1_t-pi,2*pi)-pi,'b',[0,max(t)],[q1d1,q1d1],'r',[0,max(t)],[q1d2,q1d2],'r');
% title('Recorrido primera articulación');
% legend('q1','q1_d')
% ylabel('Q1(rad)')
% xlabel('Tiempo (s)')
% xlim([0,max(t)]);
% ylim([-pi-0.1,pi+0.1]);
% grid;
% 
% subplot(3,4,[3,4]);
% plot(t, mod(q2_t-pi,2*pi)-pi);
% title('Recorrido segunda articulación');
% legend('q2','q2_d')
% ylabel('Q2(rad)')
% xlabel('Tiempo (s)')
% xlim([0,max(t)]);
%  ylim([-pi-0.1,pi+0.1]);
% grid;
% 
% subplot(3,4,5);
% plot(t, x1);
% title('Posición 1ªart en x');
% ylabel('x_1 (m)')
% xlabel('Tiempo (s)')
% xlim([0,max(t)]);
% ylim([-2-0.1,2+0.1]);
% grid;
% 
% subplot(3,4,6);
% plot(t, y1, [0,max(t)], [1,1]);
% title('Posición 1ªart en y');
% legend('y1','y1_d')
% ylabel('y_1 (m)')
% xlabel('Tiempo (s)')
% xlim([0,max(t)]);
% ylim([-2-0.1,2+0.1]);
% grid;
% 
% subplot(3,4,7);
% plot(t, x2);
% title('Posición 2ªart en x');
% ylabel('x_2 (m)')
% xlabel('Tiempo (s)')
% xlim([0,max(t)]);
% ylim([-2-0.1,2+0.1]);
% grid;
% 
% subplot(3,4,8);
% plot(t, y2, [0,max(t)], [1,1]);
% title('Posición 2ªart en y');
% legend('y2','y2_d')
% ylabel('y_2 (m)')
% xlabel('Tiempo (s)')
% xlim([0,max(t)]);
% ylim([-2-0.1,2+0.1]);
% grid;
% 
% subplot(3,4,[9,10,11,12]);
% plot(t, F_t);
% title('Torque de entrada');
% ylabel('F (Nm)')
% xlabel('Tiempo (s)')
% xlim([0,max(t)]);
% 
% grid;
