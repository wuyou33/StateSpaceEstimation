function noise = sampleGamma(noiseDataSet, N)
%% Generate N samples of a noise source specified by the NoiseDS data structure
%%
    alpha = noiseDataSet.alpha;
    beta = noiseDataSet.beta;

    if (alpha==1)
        noise = -log(1-rand(1,N))*beta;
        return
    end

    flag=0;

    if (alpha<1)
        flag=1;
        alpha=alpha+1;
    end

    gamma=alpha-1;
    eta=sqrt(2.0*alpha-1.0);
    c=.5-atan(gamma/eta)/pi;

    y(N)=0;

    for k=1:N,
        aux=-.5;
        while(aux<0)
            y(k)=-.5;
            while(y(k)<=0)
                u=rand(1,1);
                y(k) = gamma + eta * tan(pi*(u-c)+c-.5);
            end
            v=-log(rand(1,1));
            aux=v+log(1.0+((y(k)-gamma)/eta)^2)+gamma*log(y(k)/gamma)-y(k)+gamma;
        end
    end

    if (flag==1)
        noise = y.*beta.*(rand(1,N)).^(1.0/(alpha-1));
    else
        noise = y.*beta;
    end
end