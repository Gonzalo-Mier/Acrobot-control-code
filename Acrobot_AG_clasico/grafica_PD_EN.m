
% close all;

figure(1);
plot(fitness_ex/MAX_EXP);
title('Fitness per Episode')
xlabel('Episode')
ylabel('Fitness')
% ylim([0,LIM_MAX]);




q1d1 =  mod(pi-0.8-pi,2*pi)-pi;
q1d2 =  mod(pi+0.8-pi,2*pi)-pi;

figure('units','normalized','outerposition',[0 0 1 1]);

subplot(3,4,[1,2]);
plot([1:length(q1_t)]*dt, mod(q1_t-pi,2*pi)-pi,'b',[0,t_total],[q1d1,q1d1],'r',[0,t_total],[q1d2,q1d2],'r');
title('Angle of the first link');
legend('q1','q1_d')
ylabel('Q1(rad)')
xlabel('Time (s)')
xlim([0,t_total]);
ylim([-pi-0.1,pi+0.1]);
grid;

subplot(3,4,[3,4]);
plot([1:length(q2_t)]*dt, mod(q2_t-pi,2*pi)-pi,[1:length(q2d_t)]*dt_controller,q2d_t);
title('Angle of the second link');
legend('q2','q2_d')
ylabel('Q2(rad)')
xlabel('Time (s)')
xlim([0,t_total]);
 ylim([-pi-0.1,pi+0.1]);
grid;

subplot(3,4,5);
plot([1:length(q1_t)]*dt, x1);
title('Position on x of the first link');
ylabel('x_1 (m)')
xlabel('Time (s)')
xlim([0,t_total]);
ylim([-2-0.1,2+0.1]);
grid;

subplot(3,4,6);
plot([1:length(q1_t)]*dt, y1, [0,t_total], [1,1]);
title('Position on y of the first link');
legend('y1','y1_d')
ylabel('y_1 (m)')
xlabel('Time (s)')
xlim([0,t_total]);
ylim([-2-0.1,2+0.1]);
grid;

subplot(3,4,7);
plot([1:length(q2_t)]*dt, x2);
title('Position on x of the second link');
ylabel('x_2 (m)')
xlabel('Time (s)')
xlim([0,t_total]);
ylim([-2-0.1,2+0.1]);
grid;

subplot(3,4,8);
plot([1:length(q1_t)]*dt, y2, [0,t_total], [1,1]);
title('Position on y of the second link');
legend('y2','y2_d')
ylabel('y_2 (m)')
xlabel('Time (s)')
xlim([0,t_total]);
ylim([-2-0.1,2+0.1]);
grid;

subplot(3,4,[9,10,11,12]);
plot([1:length(F_t)]*dt, F_t);
title('Input torque');
ylabel('F (Nm)')
xlabel('Time (s)')
xlim([0,t_total]);

grid;



% mkdir(Carpeta);
% cd(Carpeta);
% print('grafica_completa','-djpeg')
% cd ..




