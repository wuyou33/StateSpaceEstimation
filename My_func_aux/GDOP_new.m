function [K] = GDOP_new(in_struct,dT,visual)
%Функция выполняет расчет геометрического фактора
%Вход - структура out
%out.vis - признак видимости
%out.cH  - косинусы
%Выход - К - геометрический фактор

[nsmax] = size(in_struct,1);%количество спутников
N_t_max=size(in_struct(1).vis,2);%размерность по времени

vision=zeros(nsmax,N_t_max);
% cH=zeros(nsmax*3,N_t_max);

for ns = 1:nsmax
    vision(ns,:) = in_struct(ns).vis(:);
%     cH( (ns-1)*3+1:ns*3, : )     = in_struct(ns).cH;   
end

% sum_vis = zeros(1,N_t_max);     %количество видимых спутников. для задания размерности массива
% vision = vision';               %транспонирование матрицы


%алгоритм вычисления кол-ва видимых спутников и их номеров в определнный
%момент времени

sum_vis=sum(vision); %количество видимых спутников. для задания размерности массива
K=zeros(1,N_t_max); %пустой массив геометрического фактора

%алгоритм вычисления геометрического фактора

for i=1:N_t_max
    num_vis=find( vision(:,i) );
    H=ones(length(num_vis),4);
    
    for j=1:length(num_vis)
        H(j,1:3)=-in_struct( num_vis(j) ).cH(:,i);
%         -cH( ( num_vis(j)-1 )*3+1:num_vis(j)*3 ,i);
    end
    
    K(i) = sqrt( trace( inv(H'*H) ) );
end;


%{
% for i = 1:maximum
%     for j = 1:sum_vis(i)
%         for k = 1:4 
%             
%             switch k
%                 case 1
%                     H(j,k) = -cH(num_vis(i,j),i).cos_a;
%                 case 2
%                     H(j,k) = -cH(num_vis(i,j),i).cos_b;
%                 case 3
%                     H(j,k) = -cH(num_vis(i,j),i).cos_y;
%                 case 4
%                     H(j,k) = 1;
%             end
%         end
%     end
%     
%     K(i) = sqrt(trace((H'*H)^-1));
%     H=[];
% end
%}

if visual
  figure; 
  subplot(2,1,1), plot((1:N_t_max).*dT/60,sum_vis);  
  xlabel('Время моделирования, t [мин]');    
  ylabel('кол-во видемых спутников');
  subplot(2,1,2), plot((1:N_t_max).*dT/60, (K) );
  xlabel('Время моделирования, t [мин]');    
  ylabel('GDOP');
end;
