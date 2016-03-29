%функция для выбора 3х случайных НКА из всех видимых для проверки
%работоспособности алгоритма при НКА<4
%меняем(добавляем НКА) 1 раз в 5 сек (time_sat_change)
function Visible_satellites=NKA_selectioon(Visible_satellites,observ_quant,T_dT_sc,time_sat_change)


num_sat=find( Visible_satellites(:,1,1) );
n_vis=length(num_sat);
combinations=combnk(1:n_vis,4);
select_sat_num=ceil( ( size(combinations,1) )*rand(1) );
select_sat_num=num_sat( combinations(select_sat_num,:) );

select_sat_num=[8,20,22,23];
% select_sat_num=[8,9,22,23];

%{
select_sat_num=round( ( length(num_sat)-1 ) *rand(3,1) );
for i=1:2
    flag=1;
    while flag
        double_num=find( select_sat_num==select_sat_num(i) );
        ln=length( double_num );
        if ln>1
            select_sat_num(i)=round( ( length(num_sat)-1 ) *rand(1) );
        else
            flag=0;
        end;
        
    end;
    
end;

select_sat_num=num_sat( select_sat_num+1 );
%}

Visible_satellites(:,:,1)=0;
Visible_satellites(select_sat_num,:,1)=1;


%{
Visible_satellites(:,1,1)=0;
Visible_satellites(select_sat_num,1,1)=1;


for i=2:observ_quant
    num_sat2=find( Visible_satellites(:,i,1) );
    
    if ( length(num_sat)~=length(num_sat2) )  % меняем НКА, если добавился/отвалился НКА или сменилась секунда
        
        for j=1:3 %проверяем, не отвалился ли 1 из выбраных НКА
           n_vis=find( num_sat2==select_sat_num(j) );
           if isempty( n_vis )
              flag=1; 
              while flag
                new_NKA=ceil( length(num_sat2)*rand(1) ); 
                new_NKA=num_sat2(new_NKA);
                if isempty( find( select_sat_num==new_NKA ) ) 
                    flag=0;
                end;
                
              end;
              
              select_sat_num(j)=new_NKA;
              
           end;
           
        end;      
    end;
    
    Visible_satellites(:,i,1)=0;
    Visible_satellites(select_sat_num,i,1)=1;
    
    if ( mod( i,(time_sat_change/T_dT_sc) )==1 )
        %меняем НКА, если никакой не пропал из виду
        n_vis=ceil( 3*rand(1) );%номер НКА, который будем менять
        flag=1; 
        while flag
            new_NKA=ceil( length(num_sat2)*rand(1) ); 
            new_NKA=num_sat2(new_NKA);
            if isempty( find( select_sat_num==new_NKA ) ) 
                flag=0;
            end; 
        end;
        select_sat_num(n_vis)=new_NKA;         
        
        Visible_satellites(select_sat_num,i,1)=1;
    end;
    
 
        
    num_sat=num_sat2;
end;
%}

