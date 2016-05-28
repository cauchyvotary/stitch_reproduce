function [x] = bundle_lm_solver(calFun, Cp_inv,  x0)

epsilon1 = 1e-8; 
epsilon2 = 1e-8;
tau = 1e-3;
kmax = 2000;
k = 0; 
nu = 2; 
x = x0; 
[F, r, J] = calFun(x0);
F
A = J' * J;
g = J' * r;
found = norm(g, 'inf') <= epsilon1;
mu = tau * max(diag(A));
while ~found  && k < kmax;
    k = k + 1;  
    
    h_lm = (A + mu * Cp_inv) \ (-g);
    %h_lm = (A + mu * inv( diag(diag(A)) ) ) \ (-g);
    
    if norm(h_lm) <= epsilon2 * (norm(x) + epsilon2)
        found = true;
    else
        x_new = x + h_lm;
        [F_new, r_new, J_new] = calFun(x_new);
        varrho = (F - F_new) / (0.5*h_lm' * (mu * h_lm - g));
        if varrho > 0
            x = x_new;
            F = F_new;
            J = J_new;
            r = r_new;
            A = J_new' * J_new;
            g = J_new' * r_new;
            found = norm(g, 'inf') < epsilon1;
            mu = mu * max(1/3, 1 - (2*varrho - 1)^3);
            nu = 2;
        else
            mu = mu * nu; nu = 2*nu;
        end
    end
    
end
F
k

