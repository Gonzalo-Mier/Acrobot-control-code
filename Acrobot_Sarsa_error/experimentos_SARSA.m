% MAX_EXP = 10;
% clear
% 
% MAX_EXP = 3;


% alpha = 0.2/48;
% lambda = 0.9;
% gamma = 1;
% eps = 0;
% MaxIt=600;
% error_m = 0;0.05;

% alpha = 0.2/48;
% gamma = 0.9;
% eps = 0;
% MaxIt=200;
% error_m = 0;0.05;
% t_total = 4000;




for experimento = 1:MAX_EXP
    experimento
    Sarsa_Acrobot
    it_S_exp(experimento,:) = it_t;
end

% close all

for experimento = 1:MAX_EXP
    for i=1:MaxIt-9
        it_S_exp_median(experimento,i) = sum(it_S_exp(experimento,i:i+9))/10;
    end
end

% figure(1);
% for experimento = 1:MAX_EXP
%     plot(it_S_exp(experimento,:));
%     title('Pasos por iteración')
%     xlabel('Iteración')
%     ylabel('Pasos necesitados para completar la iteración')
%     hold on;
% end
% hold off;
% 
% figure(2);
% 
% for experimento = 1:MAX_EXP
%     plot(it_S_exp_median(experimento,:));
%     title('Pasos por iteración suavizado en 10 pasos')
%     xlabel('Iteración')
%     ylabel('Pasos necesitados para completar la iteración')
%     hold on;
% end
% hold off;
% 
% 
% figure(3);
% plot(sum(it_S_exp)/MAX_EXP);
% title('Media de pasos por iteración en 10 iteraciones')
% xlabel('Iteración')
% ylabel('Pasos necesitados para completar la iteración')
% 
% figure(4);
% plot(sum(it_S_exp_median)/MAX_EXP);
% title('Media de pasos por iteración en 10 iteraciones con suavizado en 10 pasos')
% xlabel('Iteración')
% ylabel('Pasos necesitados para completar la iteración')
% 
% 
% figure(5);
% plot(1:MaxIt,sum(it_S_exp)/MAX_EXP,1:MaxIt,max(it_S_exp),1:MaxIt,min(it_S_exp));
% title('Media, máximo y mínimo de pasos por iteración en 10 iteraciones')
% legend('Pasos medio','Pasos máximos','Pasos mínimo')
% xlabel('Iteración')
% ylabel('Pasos necesitados para completar la iteración')


% 
% for i = 1:5
%     saveas(i,strcat('Figuras/SARSA_',char(string(i))),'jpeg')
% end
% 
% close all
% 
% [A,B]=min(min(it_S_exp))



