%Author: Jitesh Jhawar
%Date:25/09/2018
%This code runs stochastic simulation algorithm for ternary interaction
%individual based model

%initialize species, reaction propensities
mu  = 6;  % total number of reactions
N = 800;  % system size
C = zeros(mu,1);
A = C;
s = 0.05;  %spontaneous reaction rate
r1 = 0.005;  % pairwise interaction rate
r2 = 0.21;  % ternary interaction rate
C(1:2) = s;
C(3:4) = r1;
C(5:6) = r2;
T = 0;
steps = 50000000;   %number of simulation steps 
Tint = 1/(s);
Tend = steps*Tint;
rel = 1;    %number of repetitions
%Number of individuals in a given state storing matrices
S = zeros(steps,rel);
S1 = S; S2 = S1; sum_all = S1; tSample = S1;

% Loop for doing multiple repetions with different initial condtions
for iter = 1:rel
    iter
    X1 = ceil(rand(1,1)*N); X2 = N-X1;
    n = 0;
    T = 0;
    loop = 0;
    Tprint = 0.01;
    % Stochastic Simulation begins
    while (T < Tend)
        %         loop = loop + 1
        %Reaction propensities
        A(1) = C(1)*X2/N;  %1
        A(2) = C(2)*X1/N;  %2
        A(3) = C(3)*X1*X2/(N^2);   %3
        A(4) = C(4)*X1*X2/(N^2);   %4
        A(5) = C(5)*(X1^2)*X2/(N^3);   %5
        A(6) = C(6)*X1*(X2^2)/(N^3);   %6
        
        A0 = sum(A);
        A0 = A0(1,1);
        
        %Generating random numbers
        R(1) = rand(1,1);
        
        T = T + (log(1/R(1)))/A0;
        % Store system state at regular intervals (Tint)
        if T > Tprint
            n = n  + 1;
            S1(n,iter) = X1; S2(n,iter) = X2; sum_all(n,iter) = (X1+X2)/N;
            S(n,iter) = (X1 - X2)/N;
            tSample(n,iter) = T;
            Tprint = Tprint + Tint;
        end
        R(2) = rand(1,1);
        R2A0 = R(2)*A0;
        Sum = 0;
        nu = 1;
        %selecting next reaction that will occur
        while Sum <= R2A0 %&& nu <= 24
            mu = nu;
            Sum = Sum + A(nu);
            nu = nu + 1;
        end
        % Update states according to the current reaction
        if mu == 1 || mu == 4 || mu == 5
            X2 = X2 - 1;
            X1 = X1 + 1;
        elseif mu == 2 || mu == 3 || mu == 6
            X1 = X1 - 1;
            X2 = X2 + 1;
        end
    end
end

% Rescaling time to compare with SDE simulation
tSample = tSample/N;
% Plotting outputs
figure,
plot(tSample(:,iter),S(:,iter))
xlabel('time')
ylabel('Polarization')
% ylim([0,1])

figure,
nbins=25;
hist(S(:),nbins)
% xlim([0,1])

