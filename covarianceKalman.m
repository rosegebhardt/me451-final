function [x_est,x_pred,K_est,K_pred] = covarianceKalman(x_prev,K_prev,y,u,A,B,C,Qv,Qw)

% x_prev = x(n|n-1), K_prev = K(n,n-1)
% x_est = x(n|n), K_est = K(n,n)
% x_pred = x(n+1|n), K_pred = K(n+1|n)

G = K_prev*C'*pinv(C*K_prev*C' + Qw);
alpha = y - C*x_prev;
x_est = x_prev + G*alpha;
K_est = K_prev - G*C*K_prev;

if (min(eig(K_est)) < 0)
    warning("WARNING: Condition K(n,n) > 0 failed.")
end

K_pred = A*K_est*A' + Qv;
x_pred = A*x_est + B*u;

end