function[koef_id, koef_abd, koef_ik, koef_abk] = koef 
% ��� �������:koef  
% ������� ������������� ��� �����  ������������� ��� ������� ������� 
koef_id  = [ 0, 0, 0, 0, 1 ; % 1 
 0, 0, 0, 0, 2 ; % 2 
-2, 0, 2, 0, 1 ; % 3 
 2, 0,-2, 0, 0 ; % 4 
-2, 0, 2, 0, 2 ; % 5 
 1,-1, 0,-1, 0 ; % 6 
 0,-2, 2,-2, 1 ; % 7 
 2, 0,-2, 0, 1 ; % 8 
 0, 0, 2,-2, 2 ; % 9 
 0, -1, 0, 0, 0 ; % 10  % ���������� 
 0, 1, 2,-2, 2 ; % 11 
 0,-1, 2,-2, 2 ; % 12 
 0, 0, 2,-2, 1 ; % 13 
 -2, 0, 0,2, 0 ; % 14  % ���������� 
 0, 0, 2,-2, 0 ; % 15 
 0, 2, 0, 0, 0 ; % 16 
 0, 1, 0, 0, 1 ; % 17 
 0, 2, 2,-2, 2 ; % 18 
 0,-1, 0, 0, 1 ; % 19 
-2, 0, 0, 2, 1 ; % 20 
 0,-1, 2,-2, 1 ; % 21 
 2, 0, 0,-2, 1 ; % 22 
 0, 1, 2,-2, 1 ; % 23 
 1, 0, 0,-1, 0 ; % 24 
 2, 1, 0,-2, 0 ; % 25 
 0, 0,-2, 2, 1 ; % 26 
 0, 1,-2, 2, 0 ; % 27 
 0, 1, 0, 0, 2 ; % 28 
-1, 0, 0, 1, 1 ; % 29 
 0, 1, 2,-2, 0 ];% 30 
  
koef_abd  = [ -171996.0,-174.2, 92025.0, 8.9; % 1 
2062.0,   0.2,  -895.0, 0.5; % 2 
  46.0,   0.0,   -24.0, 0.0; % 3 
  11.0,   0.0,     0.0, 0.0; % 4 
  -3.0,   0.0,     1.0, 0.0; % 5 
  -3.0,   0.0,     0.0, 0.0; % 6 
  -2.0,   0.0,     1.0, 0.0; % 7 
   1.0,   0.0,     0.0, 0.0; % 8 
-13187.0,  -1.6,  5736.0,-3.1; % 9 
 -1426.0,  3.4,    54.0,-0.1; % 10   % ���������� 
  -517.0,   1.2,   224.0,-0.6; % 11 
   217.0,  -0.5,   -95.0, 0.3; % 12 
   129.0,   0.1,   -70.0, 0.0; % 13 
   -48.0,   0.0,     1.0, 0.0; % 14 % ���������� 
   -22.0,   0.0,     0.0, 0.0; % 15 
    17.0,  -0.1,     0.0, 0.0; % 16 
    -15.0,   0.0,     9.0, 0.0; % 17 
  -16.0,   0.1,     7.0, 0.0; % 18 
  -12.0,   0.0,     6.0, 0.0; % 19 
     -6.0,   0.0,     3.0, 0.0; % 20 
   -5.0,   0.0,     3.0, 0.0; % 21 
    4.0,   0.0,    -2.0, 0.0; % 22 
    4.0,   0.0,    -2.0, 0.0; % 23 
     -4.0,   0.0,     0.0, 0.0; % 24 
    1.0,   0.0,     0.0, 0.0; % 25 
    1.0,   0.0,     0.0, 0.0; % 26 
   -1.0,   0.0,     0.0, 0.0; % 27 
    1.0,   0.0,     0.0, 0.0; % 28 
      1.0,   0.0,     0.0, 0.0; % 29 
   -1.0,   0.0,     0.0, 0.0];  % 30 
koef_ik = [ 0, 0, 2, 0, 2;  % 31 
 1, 0, 0, 0, 0;  % 32 
 0, 0, 2, 0, 1;  % 33 
 1, 0, 2, 0, 2;  % 34 
 1, 0, 0,-2, 0;  % 35 
-1, 0, 2, 0, 2;  % 36 
 0, 0, 0, 2, 0;  % 37 
 1, 0, 0, 0, 1;  % 38 
-1, 0, 0, 0, 1;  % 39 
-1, 0, 2, 2, 2;  % 40 
 1, 0, 2, 0, 1;  % 41 
 0, 0, 2, 2, 2;  % 42 
 2, 0, 0, 0, 0;  % 43 
 1, 0, 2,-2, 2;  % 44 
 2, 0, 2, 0, 2;  % 45 
 0, 0, 2, 0, 0;  % 46 
-1, 0, 2, 0, 1;  % 47 
-1, 0, 0, 2, 1;  % 48 
 1, 0, 0,-2, 1;  % 49 
-1, 0, 2, 2, 1;  % 50 
 1, 1, 0,-2, 0;  % 51 
 0, 1, 2, 0, 2;  % 52 
 0,-1, 2, 0, 2;  % 53 
 1, 0, 2, 2, 2;  % 54 
 1, 0, 0, 2, 0;  % 55 
 2, 0, 2,-2, 2;  % 56 
 0, 0, 0, 2, 1;  % 57 
 0, 0, 2, 2, 1;  % 58 
 1, 0, 2,-2, 1;  % 59 
 0, 0, 0,-2, 1;  % 60 
 1,-1, 0, 0, 0;  % 61 
 2, 0, 2, 0, 1;  % 62 
 0, 1, 0,-2, 0;  % 63 
 1, 0,-2, 0, 0;  % 64 
 0, 0, 0, 1, 0;  % 65 
 1, 1, 0, 0, 0;  % 66 
 1, 0, 2, 0, 0;  % 67 
 1,-1, 2, 0, 2;  % 68 
 -1,-1, 2, 2, 2;  % 69 
 -2, 0, 0, 0, 1;  % 70 
  3, 0, 2, 0, 2;  % 71 
  0,-1, 2, 2, 2;  % 72 
  1, 1, 2, 0, 2;  % 73 
 -1, 0, 2,-2, 1;  % 74 
  2, 0, 0, 0, 1;  % 75 
  1, 0, 0, 0, 2;  % 76 
  3, 0, 0, 0, 0;  % 77 
  0, 0, 2, 1, 2;  % 78 
 -1, 0, 0, 0, 2;  % 79 
  1, 0, 0,-4, 0;  % 80 
  -2, 0, 2, 2, 2;  % 81 
  -1, 0, 2, 4, 2;  % 82 
   2, 0, 0,-4, 0;  % 83 
  1, 1, 2,-2, 2;  % 84 
  1, 0, 2, 2, 1;  % 85 
 -2, 0, 2, 4, 2;  % 86 
 -1, 0, 4, 0, 2;  % 87 
  1,-1, 0,-2, 0;  % 88 
  2, 0, 2,-2, 1;  % 89 
  2, 0, 2, 2, 2;  % 90 
  1, 0, 0, 2, 1;  % 91 
  0, 0, 4,-2, 2;  % 92 
  3, 0, 2,-2, 2;  % 93 
  1, 0, 2,-2, 0;  % 94 
  0, 1, 2, 0, 1;  % 95 
 -1,-1, 0, 2, 1;  % 96 
  0, 0,-2, 0, 1;  % 97 
  0, 0, 2,-1, 2;  % 98 
  0, 1, 0, 2, 0;  % 99 
  1, 0,-2,-2, 0;  % 100 
  0,-1, 2, 0, 1;  % 101 
  1, 1, 0,-2, 1;  % 102 
  1, 0,-2, 2, 0;  % 103 
   2, 0, 0, 2, 0;  % 104 
   0, 0, 2, 4, 2;  % 105 
   0, 1, 0, 1, 0];   % 106 
koef_abk = [-2274.0, -0.2, 977.0, -0.5;  % 31 
    712.0,  0.1,  -7.0,  0.0;  % 32 
   -386.0, -0.4, 200.0,  0.0;  % 33 
   -301.0,  0.0, 129.0, -0.1;  % 34 
   -158.0,  0.0,  -1.0,  0.0;  % 35 
    123.0,  0.0, -53.0,  0.0;  % 36 
     63.0,  0.0,  -2.0,  0.0;  % 37 
     63.0,  0.1, -33.0,  0.0;  % 38 
    -58.0, -0.1,  32.0,  0.0;  % 39 
    -59.0,  0.0,  26.0,  0.0;  % 40 
    -51.0,  0.0,  27.0,  0.0;  % 41 
    -38.0,  0.0,  16.0,  0.0;  % 42 
     29.0,  0.0,  -1.0,  0.0;  % 43 
     29.0,  0.0, -12.0,  0.0;  % 44 
    -31.0,  0.0,  13.0,  0.0;  % 45 
     26.0,  0.0,  -1.0,  0.0;  % 46 
     21.0,  0.0, -10.0,  0.0;  % 47 
     16.0,  0.0,  -8.0,  0.0;  % 48 
    -13.0,  0.0,   7.0,  0.0;  % 49 
    -10.0,  0.0,   5.0,  0.0;  % 50 
     -7.0,  0.0,   0.0,  0.0;  % 51 
      7.0,  0.0,  -3.0,  0.0;  % 52 
     -7.0,  0.0,   3.0,  0.0;  % 53 
     -8.0,  0.0,   3.0,  0.0;  % 54 
      6.0,  0.0,   0.0,  0.0;  % 55 
      6.0,  0.0,  -3.0,  0.0;  % 56 
     -6.0,  0.0,   3.0,  0.0;  % 57 
     -7.0,  0.0,   3.0,  0.0;  % 58 
      6.0,  0.0,  -3.0,  0.0;  % 59 
     -5.0,  0.0,   3.0,  0.0;  % 60 
      5.0,  0.0,   0.0,  0.0;  % 61 
     -5.0,  0.0,   3.0,  0.0;  % 62 
     -4.0,  0.0,   0.0,  0.0;  % 63 
      4.0,  0.0,   0.0,  0.0;  % 64 
     -4.0,  0.0,   0.0,  0.0;  % 65 
     -3.0,  0.0,   0.0,  0.0;  % 66 
      3.0,  0.0,   0.0,  0.0;  % 67 
     -3.0,  0.0,   1.0,  0.0;  % 68 
     -3.0,  0.0,   1.0,  0.0;  % 69 
     -2.0,  0.0,   1.0,  0.0;  % 70 
     -3.0,  0.0,   1.0,  0.0;  % 71 
     -3.0,  0.0,   1.0,  0.0;  % 72 
      2.0,  0.0,  -1.0,  0.0;  % 73 
     -2.0,  0.0,   1.0,  0.0;  % 74 
      2.0,  0.0,  -1.0,  0.0;  % 75 
     -2.0,  0.0,   1.0,  0.0;  % 76 
        2.0,  0.0,   0.0,  0.0;  % 77 
        2.0,  0.0,  -1.0,  0.0;  % 78 
        1.0,  0.0,  -1.0,  0.0;  % 79 
       -1.0,  0.0,   0.0,  0.0;  % 80 
        1.0,  0.0,  -1.0,  0.0;  % 81 
       -2.0,  0.0,   1.0,  0.0;  % 82 
       -1.0,  0.0,   0.0,  0.0;  % 83 
        1.0,  0.0,  -1.0,  0.0;  % 84 
       -1.0,  0.0,   1.0,  0.0;  % 85 
       -1.0,  0.0,   1.0,  0.0;  % 86 
        1.0,  0.0,   0.0,  0.0;  % 87 
        1.0,  0.0,   0.0,  0.0;  % 88 
        1.0,  0.0,  -1.0,  0.0;  % 89 
       -1.0,  0.0,   0.0,  0.0;  % 90 
       -1.0,  0.0,   0.0,  0.0;  % 91 
        1.0,  0.0,   0.0,  0.0;  % 92 
        1.0,  0.0,   0.0,  0.0;  % 93 
       -1.0,  0.0,   0.0,  0.0;  % 94 
        1.0,  0.0,   0.0,  0.0;  % 95 
        1.0,  0.0,   0.0,  0.0;  % 96 
       -1.0,  0.0,   0.0,  0.0;  % 97 
       -1.0,  0.0,   0.0,  0.0;  % 98 
       -1.0,  0.0,   0.0,  0.0;  % 99 
       -1.0,  0.0,   0.0,  0.0;  % 100 
       -1.0,  0.0,   0.0,  0.0;  % 101 
       -1.0,  0.0,   0.0,  0.0;  % 102 
       -1.0,  0.0,   0.0,  0.0;  % 103 
        1.0,  0.0,   0.0,  0.0;  % 104 
       -1.0,  0.0,   0.0,  0.0;  % 105 
        1.0,  0.0,   0.0,  0.0];   % 106 