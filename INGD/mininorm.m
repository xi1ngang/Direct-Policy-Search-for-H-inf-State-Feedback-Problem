function F = mininorm(K,delta,epsilon,A,B,Q,R,n,p)
L = 100;
Y = randsphere(1,n*p,delta)+reshape(K,1,n*p);
% F = onegradient(Y,A,B,Q,R,n,p);
F = one_K_gradient(A,B,Q,R,Y,n,p,0.0001);
Fnorm = norm(F,"fro");
JK = comput_Hinf(A,B,K,Q,R,n,p);
JK_prim = comput_Hinf(A,B,K-delta*F/Fnorm,Q,R,n,p);

while (Fnorm > epsilon) && (0.01*delta/4*Fnorm >= JK - JK_prim)
    C = Fnorm*sqrt(1-(1-Fnorm^2/(128*L^2))^2);
    r = C*rand(1);
    H = randsphere(1,n*p,r)+F;
    K_prim = K - delta*H/norm(H,'fro');
    range = abs(K - K_prim);
    Y = min(K,K_prim) + range.*rand(1);
    gra_Y = one_K_gradient(A,B,Q,R,Y,n,p,0.0001);
%     gra_Y = onegradient(Y,A,B,Q,R,n,p);
    F = min_element(F,gra_Y)';
    Fnorm = norm(F,"fro");
    JK_prim = comput_Hinf(A,B,K-delta*F/Fnorm,Q,R,n,p);
end