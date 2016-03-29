function [s_true]=Sideral_time_calc(Date_time_mass)

%

    r = 1296000.0;
    %longitude_h=2.5086667;?????

    iyear =Date_time_mass(1);
    im = Date_time_mass(2);
    iday = Date_time_mass(3);
    l_hour = Date_time_mass(4)-3;
    l_min = Date_time_mass(5);
    l_sec = Date_time_mass(6);
    longitude_h=0;

    %//20-ja stroka
    month = [31,28,31,30,31,30,31,31,30,31,30,31];

    %//double
 %   //	iday,iy,Dweek : word;
%     int i;
   % //integer
%     double t,sm,p,e,q,d,f,m,l,pl,ps,s0;%//,ihs,ims,sec;
  %  // : real;


%     //Begin
%     //{*--------------------------------------------------------------------*
%     //*                        Calculate iday and t                        *
%     //*--------------------------------------------------------------------*}
%     //35-ja stroka
%     //--	if (id <=0) or (im <= 0) or (iyear <= 0) then
%     //--	GetDate(iYear,iM,iD,Dweek);
%     //--	iday := id;

    if (im ~= 1) 
        %   //--	i = (iyear div 4)*4;
        tmp=fix(iyear/4);
        i=4*tmp;
        if (iyear == i) 
            month(2) = 29;
        end;
            
    %    //[1] - v JS s 0 numeracija massiva!
        for i=1:im-1        
            iday = iday+month(i);
        end;
    end;
  %  //50-ja stroka

    iy = iyear - 1900;
  %  //--- iday = (iday-1)+(iy-1) div 4;

  %  //document.form.test.value=iy+', d='+iday+"; hms="+l_hour+":"+l_min+":"+l_sec;
    iday = (iday-1)+ fix( (iy-1)/4 );%!!!!!!!!!!!!!!!!!!!!!!!!!! округление до целой части

   % //document.form.test.value=iy+', d='+iday+"; hms="+l_hour+":"+l_min+":"+l_sec;


   % //--t := Longint(iday) + Longint(iy)*365.0;
    t=iday + iy*365.0;
    t = (t+0.5)/36525.0;   %// {   ! 00.01.1900 12h UT}
    t = t - 1; %   // {   ! 01.01.2000 12h UT1}
%     //60-ja stroka
%     //{*--------------------------------------------------------------------*
%     //*                    Calculate mean sidereal time                    *
%     //*--------------------------------------------------------------------*}
    sm = 24110.548410 + 8640184.8128660*t + 0.093104*t*t- 0.00000620*t*t*t;

    while (sm <= 0) 
        sm = sm + 86400.0;
    end;
    while (sm > 86400)
        sm = sm - 86400.0;
    end
%     //{*--------------------------------------------------------------------*
%     //*             Calculate long and short periodic nutation             *
%     //*--------------------------------------------------------------------*}
%     //70-ja stroka
    p = pi/180.0/3600.0;
%     //{-------------------}
    e = p*(84381.448 - 46.8150*t - 0.00059*t*t + 0.0018130*t*t*t);
%     //{-------------------}
    q = p*( 450160.280 -   5.0*r*t - 482890.539*t+ 7.455*t*t + 0.0080*t*t*t);
    d = p*(1072261.3070 + 1236.0*r*t + 1105601.328*t - 6.891*t*t+ 0.0190*t*t*t);
    f = p*( 335778.8770 + 1342.0*r*t + 295263.1370*t - 13.2570*t*t+ 0.0110*t*t*t);
    m = p*(1287099.804 +  99.0*r*t+1292581.2240*t -  0.5770*t*t - 0.0120*t*t*t);
    l = p*( 485866.7330+1325.0*r*t + 715922.633*t + 31.3100*t*t+ 0.0640*t*t*t);
%     //80-ja stroka{*-------------------}


    pl =  -(17.19960 + 0.017420*t)*sin(q);

%     //pl=sin(q);
    pl = pl + (0.20620 + 0.000020*t)*sin(2.0*q);
    pl = pl +   0.00460            *sin(q+2.0*f-2.0*l);
    pl = pl +   0.00110            *sin(2.0*(l-f));
    pl = pl -   0.00030            *sin(2.0*(q+f-l));
    pl = pl-   0.00030            * sin(l-m-d);
    pl = pl-   0.00020            * sin(q-2.0*d+2.0*f-2.0*m);
    pl = pl+   0.00010            * sin(q-2.0*f+2.0*l);
    pl = pl-( 1.31870+ 0.000160*t)* sin(2.0*(q-d+f));
    pl = pl+(  0.14260-0.000340*t)* sin(m);
    pl = pl-(  0.05170-0.000120*t)* sin(2.0*q-2.0*d+2.0*f+m);
    pl = pl+(  0.02170-0.000050*t)* sin(2.0*q-2.0*d+2.0*f-m);
    pl = pl+(  0.01290+0.000010*t)* sin(q-2.0*d+2.0*f);
    pl = pl+   0.00480            * sin(2.0*(l-d));
    pl = pl-   0.00220            * sin(2.0*(f-d));
    pl = pl+(  0.00170-0.000010*t)* sin(2.0*m);
    pl = pl-   0.00150            * sin(q+m);
    pl = pl-(  0.00160-0.000010*t)* sin(2.0*(q-d+f+m));
    pl = pl-   0.00120            * sin(q-m);
    pl = pl-   0.00060            * sin(q+2.0*d-2.0*l);
    pl = pl-   0.00050            * sin(q-2.0*d+2.0*f-m);
    pl = pl+   0.00040            * sin(q-2.0*d+2.0*l);
    pl = pl+   0.00040            * sin(q-2.0*d+2.0*f+m);
    pl = pl-   0.00040            * sin(l-d);
    pl = pl+   0.00010            * sin(2.0*l+m-2.0*d);
    pl = pl+   0.00010            * sin(q+2.0*d-2.0*f);
    pl = pl-   0.00010            * sin(2.0*d-2.0*f+m);
    pl = pl+   0.00010            * sin(2.0*q+m);
    pl = pl+   0.00010            * sin(q+d-l);
    pl = pl-   0.00010            * sin(m+2.0*f-2.0*d);
%     //111-ja stroka{*------------------- }
    ps =   -(  0.22740+0.000020*t)* sin(2.0*(q+f));
    ps = ps+(  0.07120+0.000010*t)* sin(l);
    ps = ps-(  0.03860+0.000040*t)* sin(q+2.0*f);
    ps = ps-   0.03010            * sin(2.0*q+2.0*f+l);
    ps = ps-   0.01580            * sin(l-2.0*d);
    ps = ps+   0.01230            * sin(2.0*q+2.0*f-l);
    ps = ps+   0.00630            * sin(2.0*d);
    ps = ps+(  0.00630+0.000010*t)* sin(q+l);
    ps = ps-(  0.00580+0.000010*t)* sin(q-l);
    ps = ps-   0.00590            * sin(2.0*q+2.0*d+2.0*f-l);
    ps = ps-   0.00510            * sin(q+2.0*f+l);
    ps = ps-   0.00380            * sin(2.0*(q+d+f));
    ps = ps+   0.00290            * sin(2.0*l);
    ps = ps+   0.00290            * sin(2.0*q-2.0*d+2.0*f+l);
    ps = ps-   0.00310            * sin(2.0*(q+f+l));
    ps = ps+   0.00260            * sin(2.0*f);
    ps = ps+   0.00210            * sin(q+2.0*f-l);
    ps = ps+   0.00160            * sin(q+2.0*d-l);
    ps = ps-   0.00130            * sin(q-2.0*d+l);
    ps = ps-   0.00100            * sin(q+2.0*d+2.0*f-l);
    ps = ps-   0.00070            * sin(l+m-2.0*d);
    ps = ps+   0.00070            * sin(2.0*q+2.0*f+m);
    ps = ps-   0.00070            * sin(2.0*q+2.0*f-m);
    ps = ps-   0.00080            * sin(2.0*q+2.0*d+2.0*f+l);
    ps = ps+   0.00060            * sin(2.0*d+l);
    ps = ps+   0.00060            * sin(2.0*(q-d+f+l));
    ps = ps-   0.00060            * sin(q+2.0*d);
    ps = ps-   0.00070            * sin(q+2.0*d+2.0*f);
    ps = ps+   0.00060            * sin(q-2.0*d+2.0*f+l);
    ps = ps-   0.00050            * sin(q-2.0*d);
    ps = ps+   0.00050            * sin(l-m);
    ps = ps-   0.00050            * sin(q+2.0*f+2.0*l);
    ps = ps-   0.00040            * sin(m-2.0*d);
    ps = ps+   0.00040            * sin(l-2.0*f);
    ps = ps-   0.00040            * sin(d);
    ps = ps-   0.00030            * sin(l+m);
    ps = ps+   0.00030            * sin(l+2.0*f);
    ps = ps-   0.00030            * sin(2.0*q+2.0*f-m+l);
    ps = ps-   0.00030            * sin(2.0*q+2.0*d+2.0*f-m-l);
    ps = ps-   0.00020            * sin(q-2.0*l);
    ps = ps-   0.00030            * sin(2.0*q+2.0*f+3.0*l);
    ps = ps-   0.00030            * sin(2.0*q+2.0*d+2.0*f-m);
    ps = ps+   0.00020            * sin(2.0*q+2.0*f+m+l);
    ps = ps-   0.00020            * sin(q-2.0*d+2.0*f-l);
    ps = ps+   0.00020            * sin(q+2.0*l);
    ps = ps-   0.00020            * sin(2.0*q+l);
    ps = ps+   0.00020            * sin(3.0*l);
    ps = ps+   0.00020            * sin(2.0*q+d+2.0*f);
    ps = ps+   0.00010            * sin(2.0*q-l);
    ps = ps-   0.00010            * sin(l-4.0*d);
    ps = ps+   0.00010            * sin(2.0*(q+d+f-l));
    ps = ps-   0.00020            * sin(2.0*q+4.0*d+2.0*f-l);
    ps = ps-   0.00010            * sin(2.0*l-4.0*d);
    ps = ps+   0.00010            * sin(2.0*q-2.0*d+2.0*f+m+l);
    ps = ps-   0.00010            * sin(q+2.0*d+2.0*f+l);
    ps = ps-   0.00010            * sin(2.0*q+4.0*d+2.0*f-2.0*l);
    ps = ps+   0.00010            * sin(2.0*q+4.0*f-l);
    ps = ps+   0.00010            * sin(l-m-2.0*d);
    ps = ps+   0.00010            * sin(q-2.0*d+2.0*f+2.0*l);
    ps = ps-   0.00010            * sin(2.0*(q+d+f+l));
    ps = ps-   0.00010            * sin(q+2.0*d+l);
    ps = ps+   0.00010            * sin(2.0*q-2.0*d+4.0*f);
    ps = ps+   0.00010            * sin(2.0*q-2.0*d+2.0*f+3.0*l);
    ps = ps-   0.00010            * sin(l+2.0*f-2.0*d);
    ps = ps+   0.00010            * sin(q+2.0*f+m);
    ps = ps+   0.00010            * sin(q+2.0*d-m-l);
    ps = ps-   0.00010            * sin(q-2.0*f);
    ps = ps-   0.00010            * sin(2.0*q-d+2.0*f);
    ps = ps-   0.00010            * sin(2.0*d+m);
    ps = ps-   0.00010            * sin(l-2.0*f-2.0*d);
    ps = ps-   0.00010            * sin(q+2.0*f-m);
    ps = ps-   0.00010            * sin(q-2.0*d+m+l);
    ps = ps-   0.00010            * sin(l-2.0*f+2.0*d);
    ps = ps+   0.00010            * sin(2.0*(l+d));
    ps = ps-   0.00010            * sin(2.0*q+4.0*d+2.0*f);
    ps = ps+   0.00010            * sin(d+m);

%     //{*--------------------------------------------------------------------*
%     //*                    Calculate true sidereal time                    *
%     //*--------------------------------------------------------------------*}
    s0 = sm+(pl+ps)/15.0*cos(e);
%     //{*--------------------------------------------------------------------*


    s0=s0/3600.0; %//v chasah
%     //    double t_1=int(l_hour)+l_min/60+l_sec/3600;
     t_=fix(l_hour)+l_min/60.0+l_sec/3600.0;
     s_true=s0+longitude_h+(t_)*1.002737909350795;

    if (s_true<0) 
        s_true=s_true+24;
    end;
    if (s_true>=24) 
        s_true=s_true-24;
    end;

    
  s_true=  s_true*15*pi/180;

%{
    

    double iyear =UTC.YYYY;
    double im = UTC.MM;
    double iday = UTC.DD;
    double l_hour = UTC.hh;
    double l_min = UTC.mm;
    double l_sec = UTC.ss;
    longitude_h=0;

    //20-ja stroka
    int month[] = {31,28,31,30,31,30,31,31,30,31,30,31};

    //double
    //	iday,iy,Dweek : word;
    int i;
    //integer
    double t,sm,p,e,q,d,f,m,l,pl,ps,s0;//,ihs,ims,sec;
    // : real;


    //Begin
    //{*--------------------------------------------------------------------*
    //*                        Calculate iday and t                        *
    //*--------------------------------------------------------------------*}
    //35-ja stroka
    //--	if (id <=0) or (im <= 0) or (iyear <= 0) then
    //--	GetDate(iYear,iM,iD,Dweek);
    //--	iday := id;

    if (im != 1) {
        //--	i = (iyear div 4)*4;
        int tmp_=int(iyear/4);
        i=4*tmp_;
        if (iyear == i)
        {
            month[1] = 29;};
        //[1] - v JS s 0 numeracija massiva!
        for (i=1; i<=im-1; i++)
        {
            iday = iday+month[i-1];
        }
    }
    //50-ja stroka

    int iy = iyear - 1900;
    //--- iday = (iday-1)+(iy-1) div 4;

    //document.form.test.value=iy+', d='+iday+"; hms="+l_hour+":"+l_min+":"+l_sec;
    iday = (iday-1)+(iy-1)/4;

    //document.form.test.value=iy+', d='+iday+"; hms="+l_hour+":"+l_min+":"+l_sec;


    //--t := Longint(iday) + Longint(iy)*365.0;
    t=iday + iy*365.0;
    t = (t+0.5)/36525.0;   // {   ! 00.01.1900 12h UT}
    t = t - 1;    // {   ! 01.01.2000 12h UT1}
    //60-ja stroka
    //{*--------------------------------------------------------------------*
    //*                    Calculate mean sidereal time                    *
    //*--------------------------------------------------------------------*}
    sm = 24110.548410 + 8640184.8128660*t + 0.093104*t*t- 0.00000620*t*t*t;

    while (sm <= 0) {sm = sm + 86400.0;}
    while (sm > 86400) {sm = sm - 86400.0;}
    //{*--------------------------------------------------------------------*
    //*             Calculate long and short periodic nutation             *
    //*--------------------------------------------------------------------*}
    //70-ja stroka
    p = M_PI/180.0/3600.0;
    //{-------------------}
    e = p*(84381.448 - 46.8150*t - 0.00059*t*t + 0.0018130*t*t*t);
    //{-------------------}
    q = p*( 450160.280 -   5.0*r*t - 482890.539*t+ 7.455*t*t + 0.0080*t*t*t);
    d = p*(1072261.3070 + 1236.0*r*t + 1105601.328*t - 6.891*t*t+ 0.0190*t*t*t);
    f = p*( 335778.8770 + 1342.0*r*t + 295263.1370*t - 13.2570*t*t+ 0.0110*t*t*t);
    m = p*(1287099.804 +  99.0*r*t+1292581.2240*t -  0.5770*t*t - 0.0120*t*t*t);
    l = p*( 485866.7330+1325.0*r*t + 715922.633*t + 31.3100*t*t+ 0.0640*t*t*t);
    //80-ja stroka{*-------------------}


    pl =  -(17.19960 + 0.017420*t)*sin(q);

    //pl=sin(q);
    pl = pl + (0.20620 + 0.000020*t)*sin(2.0*q);
    pl = pl +   0.00460            *sin(q+2.0*f-2.0*l);
    pl = pl +   0.00110            *sin(2.0*(l-f));
    pl = pl -   0.00030            *sin(2.0*(q+f-l));
    pl = pl-   0.00030            * sin(l-m-d);
    pl = pl-   0.00020            * sin(q-2.0*d+2.0*f-2.0*m);
    pl = pl+   0.00010            * sin(q-2.0*f+2.0*l);
    pl = pl-( 1.31870+ 0.000160*t)* sin(2.0*(q-d+f));
    pl = pl+(  0.14260-0.000340*t)* sin(m);
    pl = pl-(  0.05170-0.000120*t)* sin(2.0*q-2.0*d+2.0*f+m);
    pl = pl+(  0.02170-0.000050*t)* sin(2.0*q-2.0*d+2.0*f-m);
    pl = pl+(  0.01290+0.000010*t)* sin(q-2.0*d+2.0*f);
    pl = pl+   0.00480            * sin(2.0*(l-d));
    pl = pl-   0.00220            * sin(2.0*(f-d));
    pl = pl+(  0.00170-0.000010*t)* sin(2.0*m);
    pl = pl-   0.00150            * sin(q+m);
    pl = pl-(  0.00160-0.000010*t)* sin(2.0*(q-d+f+m));
    pl = pl-   0.00120            * sin(q-m);
    pl = pl-   0.00060            * sin(q+2.0*d-2.0*l);
    pl = pl-   0.00050            * sin(q-2.0*d+2.0*f-m);
    pl = pl+   0.00040            * sin(q-2.0*d+2.0*l);
    pl = pl+   0.00040            * sin(q-2.0*d+2.0*f+m);
    pl = pl-   0.00040            * sin(l-d);
    pl = pl+   0.00010            * sin(2.0*l+m-2.0*d);
    pl = pl+   0.00010            * sin(q+2.0*d-2.0*f);
    pl = pl-   0.00010            * sin(2.0*d-2.0*f+m);
    pl = pl+   0.00010            * sin(2.0*q+m);
    pl = pl+   0.00010            * sin(q+d-l);
    pl = pl-   0.00010            * sin(m+2.0*f-2.0*d);
    //111-ja stroka{*------------------- }
    ps =   -(  0.22740+0.000020*t)* sin(2.0*(q+f));
    ps = ps+(  0.07120+0.000010*t)* sin(l);
    ps = ps-(  0.03860+0.000040*t)* sin(q+2.0*f);
    ps = ps-   0.03010            * sin(2.0*q+2.0*f+l);
    ps = ps-   0.01580            * sin(l-2.0*d);
    ps = ps+   0.01230            * sin(2.0*q+2.0*f-l);
    ps = ps+   0.00630            * sin(2.0*d);
    ps = ps+(  0.00630+0.000010*t)* sin(q+l);
    ps = ps-(  0.00580+0.000010*t)* sin(q-l);
    ps = ps-   0.00590            * sin(2.0*q+2.0*d+2.0*f-l);
    ps = ps-   0.00510            * sin(q+2.0*f+l);
    ps = ps-   0.00380            * sin(2.0*(q+d+f));
    ps = ps+   0.00290            * sin(2.0*l);
    ps = ps+   0.00290            * sin(2.0*q-2.0*d+2.0*f+l);
    ps = ps-   0.00310            * sin(2.0*(q+f+l));
    ps = ps+   0.00260            * sin(2.0*f);
    ps = ps+   0.00210            * sin(q+2.0*f-l);
    ps = ps+   0.00160            * sin(q+2.0*d-l);
    ps = ps-   0.00130            * sin(q-2.0*d+l);
    ps = ps-   0.00100            * sin(q+2.0*d+2.0*f-l);
    ps = ps-   0.00070            * sin(l+m-2.0*d);
    ps = ps+   0.00070            * sin(2.0*q+2.0*f+m);
    ps = ps-   0.00070            * sin(2.0*q+2.0*f-m);
    ps = ps-   0.00080            * sin(2.0*q+2.0*d+2.0*f+l);
    ps = ps+   0.00060            * sin(2.0*d+l);
    ps = ps+   0.00060            * sin(2.0*(q-d+f+l));
    ps = ps-   0.00060            * sin(q+2.0*d);
    ps = ps-   0.00070            * sin(q+2.0*d+2.0*f);
    ps = ps+   0.00060            * sin(q-2.0*d+2.0*f+l);
    ps = ps-   0.00050            * sin(q-2.0*d);
    ps = ps+   0.00050            * sin(l-m);
    ps = ps-   0.00050            * sin(q+2.0*f+2.0*l);
    ps = ps-   0.00040            * sin(m-2.0*d);
    ps = ps+   0.00040            * sin(l-2.0*f);
    ps = ps-   0.00040            * sin(d);
    ps = ps-   0.00030            * sin(l+m);
    ps = ps+   0.00030            * sin(l+2.0*f);
    ps = ps-   0.00030            * sin(2.0*q+2.0*f-m+l);
    ps = ps-   0.00030            * sin(2.0*q+2.0*d+2.0*f-m-l);
    ps = ps-   0.00020            * sin(q-2.0*l);
    ps = ps-   0.00030            * sin(2.0*q+2.0*f+3.0*l);
    ps = ps-   0.00030            * sin(2.0*q+2.0*d+2.0*f-m);
    ps = ps+   0.00020            * sin(2.0*q+2.0*f+m+l);
    ps = ps-   0.00020            * sin(q-2.0*d+2.0*f-l);
    ps = ps+   0.00020            * sin(q+2.0*l);
    ps = ps-   0.00020            * sin(2.0*q+l);
    ps = ps+   0.00020            * sin(3.0*l);
    ps = ps+   0.00020            * sin(2.0*q+d+2.0*f);
    ps = ps+   0.00010            * sin(2.0*q-l);
    ps = ps-   0.00010            * sin(l-4.0*d);
    ps = ps+   0.00010            * sin(2.0*(q+d+f-l));
    ps = ps-   0.00020            * sin(2.0*q+4.0*d+2.0*f-l);
    ps = ps-   0.00010            * sin(2.0*l-4.0*d);
    ps = ps+   0.00010            * sin(2.0*q-2.0*d+2.0*f+m+l);
    ps = ps-   0.00010            * sin(q+2.0*d+2.0*f+l);
    ps = ps-   0.00010            * sin(2.0*q+4.0*d+2.0*f-2.0*l);
    ps = ps+   0.00010            * sin(2.0*q+4.0*f-l);
    ps = ps+   0.00010            * sin(l-m-2.0*d);
    ps = ps+   0.00010            * sin(q-2.0*d+2.0*f+2.0*l);
    ps = ps-   0.00010            * sin(2.0*(q+d+f+l));
    ps = ps-   0.00010            * sin(q+2.0*d+l);
    ps = ps+   0.00010            * sin(2.0*q-2.0*d+4.0*f);
    ps = ps+   0.00010            * sin(2.0*q-2.0*d+2.0*f+3.0*l);
    ps = ps-   0.00010            * sin(l+2.0*f-2.0*d);
    ps = ps+   0.00010            * sin(q+2.0*f+m);
    ps = ps+   0.00010            * sin(q+2.0*d-m-l);
    ps = ps-   0.00010            * sin(q-2.0*f);
    ps = ps-   0.00010            * sin(2.0*q-d+2.0*f);
    ps = ps-   0.00010            * sin(2.0*d+m);
    ps = ps-   0.00010            * sin(l-2.0*f-2.0*d);
    ps = ps-   0.00010            * sin(q+2.0*f-m);
    ps = ps-   0.00010            * sin(q-2.0*d+m+l);
    ps = ps-   0.00010            * sin(l-2.0*f+2.0*d);
    ps = ps+   0.00010            * sin(2.0*(l+d));
    ps = ps-   0.00010            * sin(2.0*q+4.0*d+2.0*f);
    ps = ps+   0.00010            * sin(d+m);

    //{*--------------------------------------------------------------------*
    //*                    Calculate true sidereal time                    *
    //*--------------------------------------------------------------------*}
    s0 = sm+(pl+ps)/15.0*cos(e);
    //{*--------------------------------------------------------------------*


    s0=s0/3600.0; //v chasah
    //    double t_1=int(l_hour)+l_min/60+l_sec/3600;
    double t_=int(l_hour)+l_min/60.0+l_sec/3600.0;
    double s_=s0+longitude_h+(t_)*1.002737909350795;

    if (s_<0) s_+=24;
    if (s_>=24) s_-=24;

    return s_*15*M_PI/180.0;
}
  %}
