function [G] = G_angle(DN, angle)
%Функция G_angle
%Вычисляет значение коэффициента усиления антенны в направлении angle
%Входные данные: 
%Массив ДН антенны КА
%Угол между вектором напрвленным на землю и приемником/передатчиком
%Выходные данные:
%КУ антенны КА в угле angle.
G=zeros( 1,size(angle,2) );
for i = 1:size(DN,1)-1 
    A1=find(angle>DN(i));
    if not( isempty( A1 ) )
        A2=find(angle(A1)<DN(i+1));
       
        if not( isempty( A2 ) )
            A2=A1(A2);
            Koef_a     = DN(i+1) - DN(i);      %разница большего угла и меньшего
            delta_a    = angle(A2) - DN(i);                %разница между меньшим углом и а2                                   
            Koef_P     = DN(i+1,2) - DN(i,2);       %коэффициент мощности
            Koef       = Koef_P/Koef_a;               %относительный коэффициент Мощность/угол
            G(A2)          = delta_a*Koef+DN(i,2);        %КУ антенны при угле а2            
        end;
        delta_a=[];
       
    end;    
    A1=[];
    A2=[];
end;




