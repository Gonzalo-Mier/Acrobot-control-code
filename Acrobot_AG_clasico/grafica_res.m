
close all;

figure('units','normalized','outerposition',[0 0 1 1]);

subplot(3,2,1);
plot([1:length(q1_t)]*dt, mod(q1_t,2*pi),[0,t_total],[q1d,q1d]);
title('Recorrido primera articulación');
legend('q1','q1_d')
ylabel('Q1(rad)')
xlabel('Tiempo (s)')
xlim([0,t_total]);
ylim([0-0.1,2*pi+0.1]);
grid;

subplot(3,2,2);
plot([1:length(q2_t)]*dt, mod(q2_t-pi,2*pi)-pi,[0,t_total],[q2d,q2d]);
title('Recorrido segunda articulación');
legend('q2','q2_d')

ylabel('Q2(rad)')
xlabel('Tiempo (s)')
xlim([0,t_total]);
ylim([-pi-0.1,pi+0.1]);
grid;

subplot(3,2,3);
plot([1:length(q1_t)]*dt, mod(q1d - q1_t-pi,2*pi)-pi);
title('Error primera articulación');
ylabel('Q1_d - Q1 (rad)')
xlabel('Tiempo (s)')
xlim([0,t_total]);
ylim([-pi-0.1,pi+0.1]);
grid;

subplot(3,2,4);
plot([1:length(q2_t)]*dt, mod(q2d - q2_t-pi,2*pi)-pi);
title('Error segunda articulación');
ylabel('Q2_d - Q2 (rad)')
xlabel('Tiempo (s)')
ylim([-pi-0.1,pi+0.1]);

xlim([0,t_total]);

grid;

subplot(3,2,[5,6]);
plot([1:length(F_t)]*dt, F_t);
title('Torque de entrada');
ylabel('F (Nm)')
xlabel('Tiempo (s)')
xlim([0,t_total]);

grid;



mkdir(Carpeta);
cd(Carpeta);
print('grafica_completa','-djpeg')



 
figure(5);
plot(x1,y1,x2,y2,0,0,'o');
legend('Primera art.','Segunda art.','Base del robot')
title('Recorrido del brazo')
xlabel('Desplazamiento en x (m)')
ylabel('Desplazamiento en y (m)')
ylim([-.1,2.1]);
xlim([-2.1,2.1]);

saveas(5,'Mov_total','jpg')
cd ..




