
% close all;









q1d1 =  mod(pi-0.8-pi,2*pi)-pi;
q1d2 =  mod(pi+0.8-pi,2*pi)-pi;

figure('units','normalized','outerposition',[0 0 1 1]);

subplot(3,4,[1,2]);
plot([1:length(q1_t)]*dt, mod(q1_t-pi,2*pi)-pi,'b',[0,t_total],[q1d1,q1d1],'r',[0,t_total],[q1d2,q1d2],'r');
title('Recorrido primera articulación');
legend('q1','q1_d')
ylabel('Q1(rad)')
xlabel('Tiempo (s)')
xlim([0,t_total]);
ylim([-pi-0.1,pi+0.1]);
grid;

subplot(3,4,[3,4]);
plot([1:length(q2_t)]*dt, mod(q2_t-pi,2*pi)-pi,[1:length(q2d_t)]*dt_controller,q2d_t);
title('Recorrido segunda articulación');
legend('q2','q2_d')
ylabel('Q2(rad)')
xlabel('Tiempo (s)')
xlim([0,t_total]);
 ylim([-pi-0.1,pi+0.1]);
grid;

subplot(3,4,5);
plot([1:length(q1_t)]*dt, x1);
title('Posición 1ªart en x');
ylabel('x_1 (m)')
xlabel('Tiempo (s)')
xlim([0,t_total]);
ylim([-2-0.1,2+0.1]);
grid;

subplot(3,4,6);
plot([1:length(q1_t)]*dt, y1, [0,t_total], [1,1]);
title('Posición 1ªart en y');
legend('y1','y1_d')
ylabel('y_1 (m)')
xlabel('Tiempo (s)')
xlim([0,t_total]);
ylim([-2-0.1,2+0.1]);
grid;

subplot(3,4,7);
plot([1:length(q2_t)]*dt, x2);
title('Posición 2ªart en x');
ylabel('x_2 (m)')
xlabel('Tiempo (s)')
xlim([0,t_total]);
ylim([-2-0.1,2+0.1]);
grid;

subplot(3,4,8);
plot([1:length(q1_t)]*dt, y2, [0,t_total], [1,1]);
title('Posición 2ªart en y');
legend('y2','y2_d')
ylabel('y_2 (m)')
xlabel('Tiempo (s)')
xlim([0,t_total]);
ylim([-2-0.1,2+0.1]);
grid;

subplot(3,4,[9,10,11,12]);
plot([1:length(F_t)]*dt, F_t);
title('Torque de entrada');
ylabel('F (Nm)')
xlabel('Tiempo (s)')
xlim([0,t_total]);

grid;



% mkdir(Carpeta);
% cd(Carpeta);
% print('grafica_completa','-djpeg')
% cd ..




