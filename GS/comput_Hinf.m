function JK = comput_Hinf(A,B,K,Q,R,n,p)
K = reshape(K,p,n);
Acl = A-B*K;
Bcl = eye(n);
Ccl = sqrtm(Q+K'*R*K);
sys = ss(Acl,Bcl,Ccl,[],-1);
JK = hinfnorm(sys);
end