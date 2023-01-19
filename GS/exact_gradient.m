function P = exact_gradient(K_samples,A,B,Q,R,n,p)

h = 0.0001;
[m,N] = size(K_samples);
P = [];
for mk = 1:m
    grad = one_K_gradient(A,B,Q,R,K_samples(mk,:),n,p,h);
    P = [P;grad];
end

end