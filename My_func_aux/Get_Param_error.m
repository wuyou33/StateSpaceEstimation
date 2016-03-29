function [E_T,E_Fd]=Get_Param_error(Out_param_current,Out_param_mean,visible_satellites)

Nsmax=size(visible_satellites,1);
j=0;
E_T=zeros(Nsmax,1);
E_Fd=zeros(Nsmax,1);

for i=1:Nsmax
    if visible_satellites(i,1)==1
        j=j+1;
        E_T(j)=Out_param_mean(i).T_navig_resid(1)-Out_param_current(i).T_navig_resid(1);
        E_Fd(j)=Out_param_mean(i).Fdop_1(1)-Out_param_current(i).Fdop_1(1);
    end
    
end;
E_T(j+1:end)=[];
E_Fd(j+1:end)=[];
