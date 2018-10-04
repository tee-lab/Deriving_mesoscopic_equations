%Author: Jitesh Jhawar
%Date:25/09/2018
%This code simulates the numerical integration of a stochastic differential
%equation of pairwise interaction model

repetitions = 10;   %Number of repetitions
store = zeros(repetitions,1);
for rep = 1:repetitions
    n = 100000; %Number of simulation steps
    x = zeros(n,1);
    
    x(1,1) = 0.5;
    del_t = 0.01;
    i = 1;
    % simulation begins
    while i ~= n-1
        
        
        % Coefficients:
        
        Nc = 1; N = 50; e = 1/250;
        diff = (Nc/N)*(1+e-x(i)^2);
        %ensuring that diffusion is not negative by redoing the previous time step
        if (diff < 0)
            i = i - 1;
        end
        %stochastic function
        diff = (Nc/N)*(1+e-x(i)^2);
        
        %deterministic function
        drift = -1*(x(i));
        %random number
        r = randn(1,1);
        %calculating value at next time step using Euler method and ito's
        %interpretation
        x(i+1) = x(i) + (drift)*del_t + r*(diff)^0.5*(del_t^0.5);
        i = i + 1;
        
    end
    store(rep,1) = x(n,1);
end

%Plotting
figure,
plot(1:n,(x))
xlim([1,10000])
% xlabel('time','fontweight','bold')
% ylabel('OP from SDE','fontweight','bold')
% ylim([0,1])
figure,
hist(store,100)
xlabel('\rho','fontweight','bold','FontSize',22)
ylabel('Counts','fontweight','bold','FontSize',18)
% xlim([0,1])

