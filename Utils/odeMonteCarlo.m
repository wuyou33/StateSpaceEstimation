function [ output_args ] = odeMonteCarlo( input_args )
    %ODEMONTECARLO function to calculate Monte Carlo integration of an any function odeFunc between x0 and xn interval.
    rng('state', 0); %generating random number with the sequence reproducible
    
    ssum=0;
    %calculate the function of a total of 1000 times
    %printing an estimate of the integral after every 10 iterations
    for i=1:100
        for j=1:10
            x=rand;
            ssum=ssum+exp(x);
        end
        N=i*10;
        integ=ssum/N;
        %calculate the relative error; from the known value
        fprintf('N=%6.0f of,integration=%15.12f,error=%15.12f\n',N,integ,err)
    end
    
end
