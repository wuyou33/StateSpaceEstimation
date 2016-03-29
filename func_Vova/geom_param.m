function [H,a1,a2,D,cH] = geom_param(rec, trans)
%функция вычисляет номера видимых спутников на данный момент времени.
%входные данные:
%rec   - координаты КА
%trans - координаты НКА
%выходные данные:
%resut - массив данных о поазывающий наличие видимости (0 - невидно, 1 - видно)
%H     - длинна перпендикуляра опущенного из центра земли на вектор визирования
%a1    - угол между вектором визирования и вектора направленного на землю для КА
%а2    - угол между вектором визирования и вектора направленного на землю для НКА
%D     - длинна вектора визирования

N_total=length(rec.x);

% result =ones(1,N_total); %признак видимости
D = zeros(1,N_total);      %обнуление
X12=zeros(3,N_total);
X21=zeros(3,N_total);

SP1    = zeros(1,N_total);
SP2    = zeros(1,N_total);
% c      = 299792.458;

KA    = [rec.x;rec.y;rec.z];    
NKA   = [trans.x;trans.y;trans.z];

% vKA   = [rec.vx, rec.vy, rec.vz];
% vNKA  = [trans.vx, trans.vy, trans.vz];

%вычисление суммарного вектора

for i = 1:3
    DD = KA(i,:) - NKA(i,:);
    D  = D + DD.*DD;
    X12(i,:) = DD;            %вектор с направлением от НКА к КА
    X21(i,:) = -DD;           %вектор с направлением от КА к НКА
end;

D = sqrt(D);                %расстояние между объектами

R1 = sqrt(KA(1,:).^2+KA(2,:).^2+KA(3,:).^2);       %расстояние от 0 до КА
R2 = sqrt(NKA(1,:).^2+NKA(2,:).^2+NKA(3,:).^2);    %расстояние от 0 до НКА %???????


    
    for i = 1:3
        SP1 = SP1+X12(i,:).*KA(i,:);   %скалярное произведение вектора визирования и вектора 0 - приемник
    end
    a1 = acos(SP1./(D.*R1));%угол м/у линией визирования и расстоянием от КА до центра земли
    
    for i = 1:3
        SP2 = SP2+X21(i,:).*NKA(i,:);  %скалярное произведение вектора визирования и вектора 0 - передатчик
    end
    a2 = acos(SP2./(D.*R2));%угол м/у линией визирования и расстоянием от НКА до центра земли
    
    clear SP2 SP1 DD X12 X21 i R2
%-----------
%?????????????????????????????????????????????????
    H=zeros(1,N_total);
    
    pos_num_a=find(a1>pi/2);
    if not( isempty(pos_num_a) )
       H(pos_num_a) = R1(pos_num_a);       
    end;
      
    pos_num_a=find( a1<pi/2); 
    if not( isempty(pos_num_a) )
        H(pos_num_a) = sin( a1(pos_num_a) ).*R1(pos_num_a);
    end;
    
    %{
    pos_num_a=find((a1+a2)>pi/2);
    if not( isempty(pos_num_a) )
       pos_num_R=find( R1(pos_num_a)>=R2(pos_num_a) );
       
       if not( isempty(pos_num_R) )
           pos_num_R=pos_num_a(pos_num_R);
           H(pos_num_R) = R2(pos_num_R);
       end; 
       pos_num_R=[];
       
       pos_num_R=find( R1(pos_num_a)<R2(pos_num_a) ); 
       
       if not( isempty(pos_num_R) )
           pos_num_R=pos_num_a(pos_num_R);
           H(pos_num_R) = R1(pos_num_R);
       end; 
       
    end;
    
    
    pos_num_a=[];
    pos_num_a=find((a1+a2)<pi/2);
    
    if not( isempty(pos_num_a) )
        H(pos_num_a) = sin( a1(pos_num_a) ).*R1(pos_num_a);
    end;
   %}
    
    clear pos_num_a N_total R1
%-----------------
    
 %{   
    if ((a1+a2)>pi/2)  %90градусов в радианах %???????????1.570796326794897
        
      if (R1>R2)
          H = R2;
      else
          H = R1;
      end;
      
    else
        H = sin(a1)*R1;
    end;
%}    


%расчет косинусов для вычисления геометрического фактора
        
cH(1,:) = ( NKA(1,:) - KA(1,:) )./D;
cH(2,:) = ( NKA(2,:) - KA(2,:) )./D;
cH(3,:) = ( NKA(3,:) - KA(3,:) )./D;

%{
        cH.cos_a = ( NKA(1,:) - KA(1,:) )./D;
        cH.cos_b = ( NKA(2,:) - KA(2,:) )./D;
        cH.cos_y = ( NKA(3,:) - KA(3,:) )./D;
%}   
