function grad = one_K_gradient(A,B,Q,R,K,n,p,h)
% This function compute the gradient based on the analytical formula
K = reshape(K,p,n);
Acl = A-B*K;
Bcl = eye(n);
Ccl = sqrtm(Q+K'*R*K);
sys = ss(Acl,Bcl,Ccl,[],-1);
[ninf,fpeak] = hinfnorm(sys);

% Transfer function matrix
H1 = sqrtm(Q+K'*R*K);
H2 = inv(complex(cos(fpeak),sin(fpeak))*eye(n)-Acl);
T = H1*H2;
[U,S,V] = svd(T);
u1 = U(:,1);
v1 = V(:,1);
fun = @(x) expm(-H1.*x)*H2*v1*(u1)'*expm(-H1.*x);
I = integral(fun,0,Inf,'ArrayValued',true);
grad = real(R*K*(I+I')-(H2*v1*(u1)'*T*B)');
end