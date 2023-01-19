function x = min_element(F1,F2)

n = size(F1,2);

lb = min(F1,F2);
ub = max(F1,F2);
options = optimoptions('quadprog','Display','none');
x = quadprog(eye(n),[],[],[],[],[],lb,ub,[],options);

end