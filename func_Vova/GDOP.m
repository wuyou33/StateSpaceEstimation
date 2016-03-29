function [K] = GDOP(out)
%‘ункци€ выполн€ет расчет геометрического фактора
%¬ход - структура out
%out.vis - признак видимости
%out.cH  - косинусы
%¬ыход -   - геометрический фактор

[nsmax,maximum] = size(out);

for ns = 1:nsmax
    for i = 1:maximum
        vision(ns,i) = out(ns,i).vis;
        cH(ns,i)     = out(ns,i).cH;
    end
end

sum_vis = zeros(maximum,1);     %количество видимых спутников. дл€ задани€ размерности массива

vision = vision';               %транспонирование матрицы

%алгоритм вычислени€ кол-ва видимых спутников и их номеров в определнный
%момент времени


for i = 1:maximum
    j = 1;
    for ns = 1:nsmax
        if (vision(i,ns) == 1)
            sum_vis(i) = sum_vis(i)+1;
            num_vis(i,j) = ns;
            j = j+1;
        end
    end   
end


%алгоритм вычислени€ геометрического фактора

for i = 1:maximum
    for j = 1:sum_vis(i)
        for k = 1:4 
            
            switch k
                case 1
                    H(j,k) = -cH(num_vis(i,j),i).cos_a;
                case 2
                    H(j,k) = -cH(num_vis(i,j),i).cos_b;
                case 3
                    H(j,k) = -cH(num_vis(i,j),i).cos_y;
                case 4
                    H(j,k) = 1;
            end
        end
    end
    
    K(i) = sqrt(trace((H'*H)^-1));
    H=[];
end



    fig1 = figure;
    subplot(2,1,1), plot(sum_vis);
    xlabel('¬рем€ моделировани€, t')    
    ylabel('кол-во видемых спутников');
    subplot(2,1,2), plot(abs(K));
    xlabel('¬рем€ моделировани€, t')    
    ylabel('GDOP');
