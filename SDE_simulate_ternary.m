%Author: Jitesh Jhawar
%Date:25/09/2018
%This code simulates the numerical integration of a stochastic differential
%equation of ternary interaction model


n = 50000001;
x = zeros(n,1);
time = zeros(n,1);
x(1,1) = 0.1;
del_t = 0.01;
i = 1;
N = 400; % system size
s = 0.05; % spontaneous rate
c = 0.005; % pairwise interaction rate
h = 0.21; % ternary interaction rate
% simulation begins
while i ~= n-1
    
    diff = 4/N*(s+(2*c+h)*(1-x(i,1)^2)/4);
    %ensuring that diffusion is not negative by redoing the previous time step
    if (diff < 0)
        i = i - 1;
    end
    % stochastic function
    diff = 4/N*(s+(2*c+h)*(1-x(i,1)^2)/4);
    % determinisic function
    drift = -2*s*(x(i,1))+(x(i,1)*(1-x(i,1)^2)*h/2);
    %random number
    r = randn(1,1);
    %calculating value at next time step using Euler method and ito's
    %interpretation
    x(i+1) = x(i) + (drift)*del_t + r*(diff)^0.5*(del_t^0.5);
    time(i) = (i*del_t);
    i = i + 1;
end

%Plotting
figure,
plot(time,(x))
xlabel('time','fontweight','bold')
ylabel('OP from SDE','fontweight','bold')
ylim([-1,1])
figure,
hist((x),100)
xlabel('\rho','fontweight','bold','FontSize',22)
ylabel('Counts','fontweight','bold','FontSize',18)
% xlim([0,1])

