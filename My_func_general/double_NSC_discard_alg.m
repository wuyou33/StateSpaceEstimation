%% алгоритм удаления НКА с одинаковыми литерами
%оставляем НКА с большей энергетикой
function [Visible_satellites,Out_param_current]=double_NSC_discard_alg(nsmax,alm_gln,Visible_satellites,Out_param_current)

for ns=1:nsmax
        lit=alm_gln.Nn(ns);
        aa=find(alm_gln.Nn==lit);
        bb= aa==ns;
        aa(bb)=[];
        
        Vis_bouth=Visible_satellites(ns,:)+Visible_satellites(aa,:);
        
        bb=find(Vis_bouth==2, 1);
        
        if not( isempty( bb ) )
            left_greater_right=Out_param_current(ns).Power( : )>=Out_param_current(aa).Power( : );
            pos_LGR=find(left_greater_right);
            Visible_satellites(aa,pos_LGR)=0;
            Out_param_current(aa).vis( pos_LGR )=0;
            
            left_smaller_right=Out_param_current(ns).Power( : )<Out_param_current(aa).Power( : );
            pos_LSR=find(left_smaller_right);
            Visible_satellites(ns,pos_LSR)=0;
            Out_param_current(ns).vis( pos_LSR )=0;
        end;

end;