function [alm_test]=almanahReciver(PathToTheFileAlmanach)
%������� �������� ������ ��������� �������� ����� - OSCAR - 13
%���� - ����� ��������
%����� - ��������� ������ ���������� �������� ��������
%��� �������������, �������� ���������, ������ �����
%     alm_test.Lam    = 2.322003490145776;           %������� ����������� ����
%     alm_test.I      = 63.4;%0.00381591;           %���������� [�������]
%     alm_test.omegan = 4.712388980384690;%0.5233521;            %�������� ������� [�������]
%     alm_test.E      = 0.740969;%0.7005292892;%0.5005292892;         %��������������
%     alm_test.a      = 26553.4;%39800;                %������� �������
%     alm_test.tc     = 0;%34070.875;            %����� ����������� ������� ����, �
    
    
%     HEO2 - ARN 40000
%     alm_test.Lam    =83.041 ;%103.041 %[����] - ������� ����������� ����
%     alm_test.I      = 63.4;%[�������] ���������� 
%     alm_test.omegan = 270;%[�������] �������� ������� 
%     alm_test.E      = 0.740969;%  ��������������
%     alm_test.a      = 26553.4;%39800;                %������� �������
%     alm_test.tc     = 41425.8;%34070.875;            %����� ����������� ������� ����, �
%     
    
% 41425.8 41425.0

%     HEO2 - ARN 80000
%     alm_test.Lam    = 83.041;%[����] - ������� ����������� ���� 63.041
%     alm_test.I      = 63.4;%[�������] ���������� 
%     alm_test.omegan = 270;%[�������] �������� ������� 
%     alm_test.E      =0.85249;%  ��������������
%     alm_test.a      =46628.1;%                %������� �������
%     alm_test.tc     = 98538.1;%34070.875;            %����� ����������� ������� ����, �

    
%     GEO_dreyf
%     alm_test.Lam    = 70.2176;%[����] - ������� ����������� ���� 63.041
%     alm_test.I      = 0;%[�������] ���������� 
%     alm_test.omegan = 0;%[�������] �������� ������� 
%     alm_test.E      =0;%  ��������������
%     alm_test.a      =41456.5+6378.136;%                %������� �������
%     alm_test.tc     = 0;            %����� ����������� ������� ����, �
    
%     LO
    alm_test.Lam    = 66.5567;%[����] - ������� ����������� ���� 63.041
    alm_test.I      = 10;%[�������] ���������� 
    alm_test.omegan = 0;%[�������] �������� ������� 
    alm_test.E      =0;%  ��������������
    alm_test.a      =7378.14;%                %������� �������
    alm_test.tc     = 0;            %����� ����������� ������� ����, �
%     