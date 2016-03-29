function [result,H,a1,a2,D,tz,fdop] = geom_param(rec, trans, fs);
%������� ��������� ������ ������� ��������� �� ������ ������ �������.
%������� ������:
%rec   - ���������� ��
%trans - ���������� ���
%fs    - ������� ������� ���
%�������� ������:
%resut - ������ ������ � ����������� ������� ��������� (0 - �������, 1 - �����)
%H     - ������ �������������� ���������� �� ������ ����� �� ������ �����������
%a1    - ���� ����� �������� ����������� � ������� ������������� �� ����� ��� ��
%�2    - ���� ����� �������� ����������� � ������� ������������� �� ����� ��� ���
%D     - ������ ������� �����������
%tz    - ����� ��������
%fdop  - ������������ ����� �������

result = 1; %������� ���������
D      = 0;      %���������
SP1    = 0;
SP2    = 0;
c      = 299792.458;

KA    = [rec.x,rec.y,rec.z];    
NKA   = [trans.x,trans.y,trans.z];

vKA   = [rec.vx, rec.vy, rec.vz];
vNKA  = [trans.vx, trans.vy, trans.vz];

%���������� ���������� �������

for i = 1:3
    DD = KA(i) - NKA(i);
    D  = D + DD*DD;
    X12(i) = DD;            %������ � ������������ �� ��� � ��
    X21(i) = -DD;           %������ � ������������ �� �� � ���
end;

D = sqrt(D);

R1 = sqrt(KA(1)^2+KA(2)^2+KA(3)^2);    %���������� �� 0 �� ��
R2 = sqrt(NKA(1)^2+NKA(2)^2+NKA(3)^2);    %���������� �� 0 �� ���

if (R2>0)
    
    for i = 1:3
        SP1 = SP1+X12(i)*KA(i);  %��������� ������������ ������� ����������� � ������� 0 - ��������
    end
    a1 = acos(SP1/(D*R1));
    
    for i = 1:3
        SP2 = SP2+X21(i)*NKA(i);  %��������� ������������ ������� ����������� � ������� 0 - ����������
    end
    a2 = acos(SP2/(D*R2));
    
    if ((a1+a2)>1.570796326794897)  %90�������� � ��������
        
      if (R1>R2)
          H = R2;
      else
          H = R1;
      end;
      
    else
        H = sin(a1)*R1;
    end;
    
else
    result = 0;
    H = 0;
end;

%====================���������� �������� �������===========================

tz = D/c;

%===================���������� ������������� ������========================
%==========================================================================
%===========================��������� ���������============================

    S1 = 0;
    S2 = 0;
    
%==========================================================================    
    
    summ1 = sqrt(vKA(1)^2 + vKA(2)^2 + vKA(3)^2);     %������ ������� �������� ��
    summ2 = sqrt(vNKA(1)^2 + vNKA(2)^2 + vNKA(3)^2);  %������ ������� �������� ���
    
    for i = 1:3
        S1 = S1+X12(i)*vKA(i); %��������� ������������ ������� ����������� � ������� �������� ���������
    end;
    arg1 = acos(S1/(D*summ1));
    v1   = cos(arg1)*summ1;
    
    for i = 1:3
        S2 = S2+X21(i)*vNKA(i); %��������� ������������ ������� ����������� � ������� �������� �����������
    end;
    arg2 = acos(S2/(D*summ2));  %���� � ��������
    v2   = cos(arg2)*summ2;

% for i =1:3        %2� - ������ - ���� �������
%     dV(i) = vKA(i) - vNKA(i);
% end
% 
% for i =1:3
%     X12(i) = X12(i) / D;
%     X21(i) = X21(i) / D;
% end
% 
% S=0;
% for i =1:3
%     S=S+dV(i)*X12(i);
% end

fr   = 1600*10^6; 

fdop = fr*(v2 + v1)/c;