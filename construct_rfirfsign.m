function [ YYirf_new, valid ] = construct_rfirfsign(PHI, SIGMA, nirf, signrestrictionsindx);

nvar        = size(PHI,2);
nlags       = floor(size(PHI,1)/nvar);

% construt IRFs to orthogonalized but not structural shocks
[ YYirf_hat ]  = construct_rfirf(PHI, SIGMA, nirf);

 YYirf_restricted_hat = select_irf(YYirf_hat,signrestrictionsindx); 
[ psi_hat_restricted, index_phi ] = construct_psi_restricted(YYirf_restricted_hat); 


%sample unit length vector q
temp_q = randn(nvar,1);
q = temp_q/norm(temp_q);

%check whether the restrictions are satisfied
mu_dim    = size(signrestrictionsindx(:,3),1);
q_dim     = size(q,1);
Sq        = kron(eye(mu_dim),q');
Sq        = Sq(:,index_phi);

csi_hat   = (Sq*psi_hat_restricted).*signrestrictionsindx(:,3);
n_pos   = sum(csi_hat >= 0 );

 if n_pos == size(csi_hat,1)
    valid = 1;
    psi_hat_all = reshape(YYirf_hat',nvar*nvar*nirf,1);  
    Sq_all   = kron(eye(nirf*nvar),q');
    YYirf_temp = Sq_all*psi_hat_all; %responses are [var1h1;var2h1;var3h1;var4h1;...;var1hnirf;var2hnirf;var3hnirf;var4hnirf]
    
    %reshape the IRF to make it var1h1;var1h2;var1h3;...;var1hnirf; var2h1;
    %var2h2;....;var4h1;...var4hnirf]
    for i = 1:nvar
        YYirf_new((i-1)*nirf+1:i*nirf,1) = YYirf_temp([i:nvar:(nirf-1)*nvar+i],1);
    end
    
 else
     valid = 0;
     YYirf_new = [];
    
 end

end

