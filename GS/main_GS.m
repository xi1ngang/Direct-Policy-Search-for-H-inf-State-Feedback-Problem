clear;
clc;
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
K0 = [0.493136227935763	-0.136800899469019	-2.26542205030330];


% Uncomment here for the higher MIMO example 
% High dimension MIMO system 
% n = 4;
% p = 2;
% A = [1.7865 0.3912 0.8758 0.5996; 0.2756 1.3175 0.7692 0.4848; ...
%      0.4764 0.9786 1.0618 0.7591; 0.4489 0.7918 0.6014 1.7520]; 
% B = [0.1303 0.0312; 0.1309 0.0528; 0.7452 0.6727; 0.2460 0.0743];
% Q = 1.0613*eye(n);
% R = 1.1315*eye(p);
% K0 = [2.4364 12.1213 2.2337 -4.6823 2.4867 2.1718 1.5551 -2.5906];

Jm = [];
Km = [];

% Algorithm parameters
delta0 = 0.01;
epsilon0 = 100;
mu_delta = 0.5;
mu_epsilon = 0.5;
m = n*p+1;
beta = 0.5;
theta = 0.9;
% Initial points
delta = delta0;
epsilon = epsilon0;
K = K0;
Num = 100;
KK = K;
JJJ = comput_Hinf(A,B,K,Q,R,n,p);

% Algorithm loop
flag = 2; %flag = 1, using estimated gradient; flag = 2, using gradient formula
for N = 1:Num
    K_samples = randsphere(m,n*p,delta)+reshape(K,1,n*p).*ones(m,1);
    if flag == 2
        P = est_gradient(K_samples,A,B,Q,R,n,p);
    else
        P = exact_gradient(K_samples,A,B,Q,R,n,p);
    end
    F = mininorm_convexhull(P')';
    if norm(F,"fro") < epsilon
        epsilon = epsilon*mu_epsilon;
        delta = delta*mu_delta;
        tn = 0;
    else
        t = delta;
        Count = 0;
        Fnorm = norm(F,"fro");
        F_hat = F;
        J1 = comput_Hinf(A,B,K,Q,R,n,p);
        while Count == 0
            J2 = comput_Hinf(A,B,K-t*F_hat,Q,R,n,p);
            if J2 <= J1-beta*t*Fnorm
                tn = t;
                Count = 1;
            else
                t = t*theta;
            end
        end
        K = K - tn*F_hat;
    end
    J = comput_Hinf(A,B,K,Q,R,n,p);
    JJJ = [JJJ;J];
    KK = [KK;K];
    formatSpec = 'Iteration number %d, cost function value %f, delta,epsilon value (%f %f) \n';
    fprintf(formatSpec,N,J,delta,epsilon);
end


% Plot the iterations 
J_opt = 7.3475;
% J_opt = 43.26; % Optimal function value for high dimensional case
figure()
plot(1:length(JJJ),(JJJ-J_opt)/J_opt,'-','LineWidth', 3)
xlabel('Iterations','FontSize',20,'Interpreter','latex')
ylabel('$(J(K)-J^*)/J^*$','FontSize',20,'Interpreter','latex')
grid on;
set(gca,'FontSize',20)
title('Gradient Sampling Method');


