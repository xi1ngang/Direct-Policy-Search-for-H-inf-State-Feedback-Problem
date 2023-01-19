function Chi = gradient_est(K_samples,A,B,Q,R,alpha,length,n,p)

[m,N] = size(K_samples);
xi = -0.5+rand(m,N);
Chi = zeros(m,N);
for k = 1:m
    for kk = 1:N
        K1 = K_samples(k,:)+alpha*xi(k,:);
        K1(kk) = K_samples(k,kk)+0.5*alpha;
        K2 = K_samples(k,:)+alpha*xi(k,:);
        K2(kk) = K_samples(k,kk)-0.5*alpha;
        JK1 = comput_Hinf_data(A,B,K1,Q,R,n,p,length);
        JK2 = comput_Hinf_data(A,B,K2,Q,R,n,p,length);
        Chi(k,kk) = 1/alpha*(JK1 - JK2);
    end
end
end