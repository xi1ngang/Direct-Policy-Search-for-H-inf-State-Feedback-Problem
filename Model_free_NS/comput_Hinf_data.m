function gamma = comput_Hinf_data(A,B,K,Q,R,n,p,N)
B1 = eye(n);
K = reshape(K,p,n);
Ccl = sqrtm(Q+K'*R*K);
D11 = zeros(n,n);
max_iterations = 300;
tolerance = 1e-10;
gamma = estimate_Hinf(A, B, B1, Ccl, D11, K, N, max_iterations, tolerance);

end