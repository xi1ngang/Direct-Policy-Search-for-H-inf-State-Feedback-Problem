clear;
clc;
% This script implements the interpolated normalized gradient descent (INGD) method
% This script implement the gradient sampling
% Written by:
% -- 
% Xingang Guo            2022-10-01
% 
% email: xingang2@illinois.edu
% 
% Please send comments and especially bug reports to the
% above email address.

% Problem matrix (16)
n = 3;
p = 1;
A = [1 0 -5;-1 1 0;0 0 1];
B = [1;0;-1];
Q = [2 -1 0;-1 2 -1;0 -1 2];
R = 1;

% Initial condition 
K0 = [0.493136227935763	-0.136800899469019	-2.26542205030330];

% Algorithm parameters
delta0 = 1e-2;
epsilon0 = 1e-4;
delta = delta0;
epsilon = epsilon0;
K = K0;
JJ = comput_Hinf(A,B,K,Q,R,n,p);
Num = 100;
mu_delta = 0.7;
FF = [];
for N = 1:Num
    F = mininorm(K,delta,epsilon,A,B,Q,R,n,p);
    FF = [FF norm(F,'fro')];
    if norm(F,'fro') < epsilon
        delta = delta*mu_delta;
    else 
        K = K - delta*F./norm(F,"fro");
        J = comput_Hinf(A,B,K,Q,R,n,p);
    end
    JJ = [JJ J];
    formatSpec = 'Iteration number %d, cost function value %f, delta value %f \n';
    fprintf(formatSpec,N,J,delta);
end

% Plot result
Jopt = 7.3475;
figure()
plot(1:length(JJ),(JJ-Jopt)/Jopt,'-','LineWidth', 3)
xlabel('Iterations','FontSize',20,'Interpreter','latex')
ylabel('$(J(K)-J^*)J^*$','FontSize',20,'Interpreter','latex')
grid on;
set(gca,'FontSize',20)


