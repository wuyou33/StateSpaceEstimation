%% Функция перехода координат из Абсолютной Инерциальной СК в Гринвечскую СК ПЗ-90.11 
function [Coord_Out,Vel_Out]=CS_Transform_Inert2PZ(Coord_In,Vel_In,S)
global GL_W_rot_Earth

N=size(S,1);
W=[0 GL_W_rot_Earth 0;-GL_W_rot_Earth 0 0;0 0 0];%матрица угловых скоростей

Coord_Out=zeros(3,N);
Vel_Out=zeros(3,N);

for i=1:N
%     A_trans=[cos( S(i) ) sin( S(i) )  0;...
%                    -sin( S(i) ) cos( S(i) ) 0;...
%                           0          0            1];% матрица поворота
%                
% Coord_Out(:,i)=A_trans*Coord_In(:,i);
% Vel_Out(:,i)=A_trans*Vel_In(:,i)+W*Coord_Out(:,i);
    
Coord_Out(1,i)=Coord_In(1,i)*cos( S(i) )+Coord_In(2,i)*sin( S(i) );
Coord_Out(2,i)=-Coord_In(1,i)*sin( S(i) )+Coord_In(2,i)*cos( S(i) );
Coord_Out(3,i)=Coord_In(3,i);
Vel_Out(1,i)=Vel_In(1,i)*cos( S(i) )+Vel_In(2,i)*sin( S(i) )+GL_W_rot_Earth*Coord_Out(2,i);
Vel_Out(2,i)=-Vel_In(1,i)*sin( S(i) )+Vel_In(2,i)*cos( S(i) )-GL_W_rot_Earth*Coord_Out(1,i);
Vel_Out(3,i)=Vel_In(3,i);
    
end;
