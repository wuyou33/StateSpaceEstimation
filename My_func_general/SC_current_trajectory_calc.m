function [Satpos_xyz_Rec_current]=SC_current_trajectory_calc(T_Ttotal_eval,time_refresh_data,T_dT_sc,T_Tstart,T_Tend,...
    Sun,Moon,y0_sc_current)

N_15min_inter_max=ceil(T_Ttotal_eval/time_refresh_data);%сколько раз по 15 мин укладываетс€ во всЄм времени проектировани€

N_max=T_Ttotal_eval/T_dT_sc+1;% полное число точек траектории полЄта
Ti_current=T_Tstart;

Coord_Inert_sc_current=zeros(3,N_max);
Vel_Inert_sc_current=zeros(3,N_max);

%цикл по времени- интервалы по 15 мин
for j=1:N_15min_inter_max
    
    K=time_refresh_data/T_dT_sc+1;%сколько точек траектории заполнитс€ на этой итерации
   
    
    if (Ti_current+time_refresh_data)<=T_Tend
        T=Ti_current:T_dT_sc:(Ti_current+time_refresh_data);%массив времени, дл€ которого будет производитс€ решение системы ƒ”.
        k=K;
        start=(j-1)*K+2-j;%начальный интевал отсчЄта по времени 1=T_dT_sc
        stop=(j-1)*K+1+k-j;%конечное значение пееменной итрации по времени
    elseif Ti_current==T_Tend
        continue;
    else
        T=Ti_current:T_dT_sc:T_Tend;%массив времени, дл€ которого будет производитс€ решение системы ƒ”.
        k=size(T,2);
        start=(j-1)*K+2-j;
        stop=start+k-1;
    end
     
    [~,Y]=ode113( @(t,y) d_State_Space_dT(t,y,Sun(:,j) ,Moon(:,j) ), T,y0_sc_current );%непосредственно решаем систему ƒ” дл€ текущей орбиты (с ошибкой)
 %ode45
%образуем данные из массива координат и скоростей в структуру 
    Coord_Inert_sc_current( 1:3, start:stop  )=Y(1:end,1:3)';
    Vel_Inert_sc_current( 1:3, start:stop  )=Y(1:end,4:6)';    
%----------
    y0_sc_current=[Y(end,1:3)'; Y(end,4:6)'];% задаЄм начальные услови€ дл€ ƒ”  
%----------
    
    Ti_current=Ti_current+time_refresh_data;%увеличиваем врем€ на 15 мин
end;


%}


%текуща€ траектори€ движени€  ј
Satpos_xyz_Rec_current=struct('x',zeros(1,N_max),'y',zeros(1,N_max),'z',zeros(1,N_max),'vx',zeros(1,N_max),'vy',zeros(1,N_max),'vz',zeros(1,N_max));
Satpos_xyz_Rec_current.x=Coord_Inert_sc_current(1,:);
Satpos_xyz_Rec_current.y=Coord_Inert_sc_current(2,:);
Satpos_xyz_Rec_current.z=Coord_Inert_sc_current(3,:);

Satpos_xyz_Rec_current.vx=Vel_Inert_sc_current(1,:);
Satpos_xyz_Rec_current.vy=Vel_Inert_sc_current(2,:);
Satpos_xyz_Rec_current.vz=Vel_Inert_sc_current(3,:);