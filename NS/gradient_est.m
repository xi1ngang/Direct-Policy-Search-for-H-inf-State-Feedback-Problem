function Chi = gradient_est(K_samples,A,B,Q,R,alpha,n,p)
% This function compute the gradinet based on the Gupal estimation
[m,N] = size(K_samples);
xi = -0.5+rand(m,N);
Chi = zeros(m,N);
for k = 1:m
    for kk = 1:N
        K1 = K_samples(k,:)+alpha*xi(k,:);
        K1(kk) = K_samples(k,kk)+0.5*alpha;
        K2 = K_samples(k,:)+alpha*xi(k,:);
        K2(kk) = K_samples(k,kk)-0.5*alpha;
        JK1 = comput_Hinf(A,B,K1,Q,R,n,p);
        JK2 = comput_Hinf(A,B,K2,Q,R,n,p);
        Chi(k,kk) = 1/alpha*(JK1 - JK2);
    end
end
end