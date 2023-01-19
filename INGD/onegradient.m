function F = onegradient(Y,A,B,Q,R,n,p)
% Finite difference gradient 
h = 0.00001;
N = size(Y,2);
E = eye(n*p);
F = zeros(1,n);
for k = 1:N
    Yprim = Y + h*E(k,:);
    JK = comput_Hinf(A,B,Y,Q,R,n,p);
    JK_prim = comput_Hinf(A,B,Yprim,Q,R,n,p);
    F(k) = (JK_prim-JK)/h;
end

end