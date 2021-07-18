function  [teta , P] = RLS(phi_k_1 , u_k_1, e_k, teta, P, sigma)

Phi = phi_k_1;
P = P - P*Phi'*(1 + Phi * P * Phi')^(-1)*Phi*P;
teta = teta + P*Phi'/(1 + Phi * P * Phi')*(e_k + sigma*(Phi*teta-u_k_1)); 

end