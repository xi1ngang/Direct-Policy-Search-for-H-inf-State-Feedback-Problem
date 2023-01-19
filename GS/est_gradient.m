function P = est_gradient(K_samples,A,B,Q,R,n,p)

h = 0.00001;
[m,N] = size(K_samples);

E = eye(N);
P = zeros(m,N);

for k = 1:m
    for kk = 1:N
        K1 = K_samples(k,:)+h*E(kk,:);
        JK1 = comput_Hinf(A,B,K1,Q,R,n,p);
        K2 = K_samples(k,:);
        JK2 = comput_Hinf(A,B,K2,Q,R,n,p);
        P(k,kk) = (JK1-JK2)/h;
    end   
end


end