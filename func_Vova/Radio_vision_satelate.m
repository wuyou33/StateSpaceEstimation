function [vision,P_input, cH] = Radio_vision_satelate(satpos_xyz_Rec, satpos_xyz_gln)
%???????????????????
%Функцмя выполняет расчет относительно радиовидемости КА и НКА а также
%вычисляет направляющие косинусы
%Вход:
%satpos_xyz_Rec - координаты КА
%satpos_xyz_gln - координаты НКА
%Выход:
%vision - признак видимости
%cH - направляющие косинусы

global GL_Threshold %Пороговая мощность
global GL_Temp;%[K] - температура приёмника
global GL_Ae %[км] - экваториальный радиус Земли
global GL_DN_GLN %[дБ] диаграмма направленности НКА Глонасс
global GL_DN_ARN %[дБ] диаграмма направленнсти КА
global GL_DN_ARN_LO %[dB] диаграмма направленнсти КА для НО
global GL_Hmin_ion %[км] высота ниже которой считается,что КА не виден 


P_GLN = 16.627578316815743; %Мощность передатчика ГЛОНАСС[dBВт]

k=1.380648813131313e-23;%[Дж/К] - постоянная больцмана
c=2.99792458e8; %[м/c] - скорость света
F0=1602e6;  %частота 1 диапазона
  
Lam = c/F0; %длинна волны

h_min=GL_Hmin_ion;%1000;%450; - минимальная высота над поверхностью земли (из-за ионосферы)

% N0 =-201.5;
N0 =10*log10(k*GL_Temp); %[dB-Hz ??? ] - Спектральная мощность шума

Lp = 0;%4  [dB] - поляризационные потери и прочие потери

%{
%ДН антенн в [dB]
% DN_GLN = [0 10.7; 9 11.763636;10 11.881818; 11 12; 12 12; 15 12; 19 9;...
%           20 7; 21 5.5; 22 4; 23 2.7; 24 1.5; 25 0; 30 -5; 35 -4.375;...
%           40 -5.875; 45 -7.875; 50 -8.25; 55 -10.875; 60 -12.75];
      
% DN_GPS = [0 13; 5 13.8; 6 14; 8 14.7; 9 14.9; 10 15; 11 14.9; 12 14.7;...
%           14 13.7; 18 6; 20 0; 22 -10; 23 -8; 24 -5; 26 -3; 30 -2; 34 -3;...
%           36 -5; 38 -10; 40 -15];
      
% DN_ARN = [0 3; 2 3; 4 3; 6 3; 8 3; 10 3; 12 2.96; 14 2.94; 16 2.92; 18 2.935;...
%           20 2.925; 22,2.9235; 24 2.914; 29 2.9; 60 1.5; 90 -2];
%}
%{
DN_GLN = [0.0    11.0; 4.0    11.5; 8.0    13.0; 10.0   13.6; 12.0   14.0; 16.0   13.2; 20.0   10.8;... %new
    24.0    6.0; 28.0   -1.0; 30.0   -3.0; 36.0   -6.0; 44.0   -8.0; 50.0  -10.0; 52.0  -12.0; 54.0  -15.0; 60.0  -11.0; ...
    64.0   -7.0; 70.0   -6.0; 80.0   -7.0; 90.0   -10.0];

DN_ARN = [0	14.8; 5	14.5; 10	14.0; 15	13; 20	10.5; 25	8; 30	3; 35	-6; 40	-6; 45	-2; 50	0; 55	-3;...%new
    60	-6; 65	-18; 70	-12;75	-10; 80	-12; 85	-18; 90	-30];
%}


%Вычисление визуальной видимости
      
        [H,a1,a2,D,cH] = geom_param(satpos_xyz_Rec, satpos_xyz_gln);
        
        
        maxNKA_angel = max(GL_DN_GLN);
        deg_a1 = 180.*a1./pi;
        deg_a2 = 180.*a2./pi;

        vision=zeros( 1,length(satpos_xyz_Rec.x) );
        P_input=zeros( 1,length(satpos_xyz_Rec.x) );

        G_GLN=-100*ones(1,length(satpos_xyz_Rec.x) );
        G_ARN=-100*ones(1,length(satpos_xyz_Rec.x) );
        Lf=-200*ones(1,length(satpos_xyz_Rec.x) );
        D  = D.* 1000; 
            
        r_sc_current=sqrt( satpos_xyz_Rec.x(:).^2+satpos_xyz_Rec.y(:).^2+satpos_xyz_Rec.z(:).^2 )-GL_Ae;
         
        R_arn_hi=find( r_sc_current>=5000 );
        if not( isempty( R_arn_hi ) )
            DN_ARN=GL_DN_ARN;
            
            maxKA_angel = max(DN_ARN);
            minKA_angel = min(DN_ARN);
           

            pos_num_a1=find( abs( deg_a1(R_arn_hi) )<maxKA_angel(1) );
            pos_num_a1=R_arn_hi(pos_num_a1);
            a=find( abs( deg_a1(pos_num_a1) )>minKA_angel(1) );%!!!!
            pos_num_a1=pos_num_a1(a); %!!!

            if not( isempty(pos_num_a1) )
                pos_num_a1_a2=find( abs( deg_a2(pos_num_a1) )<maxNKA_angel(1) );

                if not( isempty(pos_num_a1_a2) )
                    pos_num_a1_a2=pos_num_a1(pos_num_a1_a2);
                    G_GLN(pos_num_a1_a2) = G_angle(GL_DN_GLN, deg_a2(pos_num_a1_a2) );    %Вычисление КУ антенны в угле а2
                    G_ARN(pos_num_a1_a2) = G_angle(DN_ARN, deg_a1(pos_num_a1_a2) );  %Вычисление КУ антенны в угле а1

                    
                    Lf(pos_num_a1_a2) = 20.*log10( Lam./( 4.*pi.*D(pos_num_a1_a2) ) ); % [dB] - потери при РРВ

%                     P_priem = P_GLN + G_GLN + G_ARN + Lf - Lp;  %мощность на входе приёмника 
%                     P_por   = GL_Threshold + N0; %пороговое отношение дБ-Гц

%                     pos_num_P_H1=find(P_priem>P_por);
%                     if not( isempty( pos_num_P_H1 ) )
%                         pos_num_P_H2=find( ( H( pos_num_P_H1 ) - GL_Ae )>h_min );%450
% 
%                         if not( isempty( pos_num_P_H2 ) )
%                             pos_num_P_H2=pos_num_P_H1(pos_num_P_H2);
%                             vision(pos_num_P_H2)=1;
%                             P_input(pos_num_P_H2)=P_priem(pos_num_P_H2)-N0;
%                         end;
% 
%                     end;

                end;

            end;
            
        end;
    %
        R_arn_lo=find( r_sc_current<5000 );
        if not( isempty( R_arn_lo ) )
            DN_ARN=GL_DN_ARN_LO;
            
            maxKA_angel = max(DN_ARN);
            minKA_angel = min(DN_ARN);
                       
            pos_num_a1=find( abs( deg_a1(R_arn_lo) )<maxKA_angel(1) );
            pos_num_a1=R_arn_lo(pos_num_a1);
            a=find( abs( deg_a1(pos_num_a1) )>minKA_angel(1) );%!!!!
            pos_num_a1=pos_num_a1(a); %!!!

            if not( isempty(pos_num_a1) )
                pos_num_a1_a2=find( abs( deg_a2(pos_num_a1) )<maxNKA_angel(1) );

                if not( isempty(pos_num_a1_a2) )
                    pos_num_a1_a2=pos_num_a1(pos_num_a1_a2);
                    G_GLN(pos_num_a1_a2) = G_angle(GL_DN_GLN, deg_a2(pos_num_a1_a2) );    %Вычисление КУ антенны в угле а2
                    G_ARN(pos_num_a1_a2) = G_angle(DN_ARN, deg_a1(pos_num_a1_a2) );  %Вычисление КУ антенны в угле а1

%                     D  = D.* 1000; 
                    Lf(pos_num_a1_a2) = 20.*log10( Lam./( 4.*pi.*D(pos_num_a1_a2) ) ); % [dB] - потери при РРВ
                   

                end;

            end;
            
        end;
        %}
        
        
         P_priem = P_GLN + G_GLN + G_ARN + Lf - Lp;  %мощность на входе приёмника 
         P_por   = GL_Threshold + N0; %пороговое отношение дБ-Гц

         pos_num_P_H1=find(P_priem>P_por);
            if not( isempty( pos_num_P_H1 ) )
                pos_num_P_H2=find( ( H( pos_num_P_H1 ) - GL_Ae )>h_min );%450

                if not( isempty( pos_num_P_H2 ) )
                    pos_num_P_H2=pos_num_P_H1(pos_num_P_H2);
                    vision(pos_num_P_H2)=1;
                    P_input(pos_num_P_H2)=P_priem(pos_num_P_H2)-N0;
                end;

            end;
        
        
        
        
        
        
%{      
        if ( ( abs(deg_a1)<maxKA_angel(1) ) && ( abs(deg_a2) <maxNKA_angel(1) ) ) %???????? добавил модули
        
            G_GLN = G_angle(DN_GLN, deg_a2);    %Вычисление КУ антенны в угле а2
            G_ARN = G_angle(DN_ARN, deg_a1);   %Вычисление КУ антенны в угле а1
            
            D  = D * 1000; 
            Lf = 20*log10(Lam/(4*pi*D)); % [dB] - потери при РРВ

            P_priem = P_GLN + G_GLN + G_ARN + Lf - Lp;   
            P_por   = GL_Threshold + N0; 

            if (( P_priem > P_por) && (H - GL_Ae > 500))  %????????????
                vision = 1;
            else
                vision = 0;
            end;
        
        else
            vision = 0;
        end;
%}
        
        